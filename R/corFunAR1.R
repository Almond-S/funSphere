

#' Functional AR-1 covariance structure based on covariance smoothing
#'
#' @param value smoothing parameter for symmetric tensor-product smoother (\code{gam} smoother class \code{symm.smooth})
#' @param form a one sided formula of the form ~ t, or ~ t | g,
#' specifying a time covariate t and, optionally, a grouping factor g.
#' When a grouping factor is present in form, the correlation structure is assumed
#' to apply only to observations within the same grouping level;
#' observations with different grouping levels are assumed to be uncorrelated.
#' @param fixed an optional logical value indicating whether the coefficients
#' should be allowed to vary in the optimization, or kept fixed at their initial value.
#' Defaults to \code{FALSE}, in which case the coefficients are allowed to vary.
#'
#' @return an object of class \code{corFunAR1}, representing an covariance smoother autocorrelation structure.
#' @import mgcv nlme
#' @export
#'
corFunAR1 <- function(value = c(0,0), form = ~1, fixed = FALSE,
                      working_correlation = corExp,
                      working_control = list(),
                      s_xt = list(), s_k = -1, s_m = NA, verbose = FALSE) {
  if (any(value < 0)) {
    stop("penalty parameter for covariance smoothing must be non-negative")
  }
  value <- notLog2(value)

  attr(value, "s_args") <- list(xt = s_xt, k = s_k, m = s_m)
  attr(value, "verbose") <- verbose
  class(value) <- c("corFunAR1", "corStruct")

  corDynamic(value, working_correlation = working_correlation,
             working_control = working_control,
             form = form, fixed = fixed)
}


#' @export
#' @import nlme
#' @rdname nlme::coef.corStruct
#'
coef.corFunAR1 <- function (object, unconstrained = TRUE, ...) {
  if (unconstrained) {
    if (attr(object, "fixed")) {
      return(numeric(0))
    }
    else {
      return(as.vector(object))
    }
  }
  aux <- notExp2(as.vector(object))
  aux
}

#' @export
#' @import nlme
#'
Initialize.corFunAR1 <- function(object, data, ...) {
  object <- NextMethod()
  covar <- getCovariate(object)[[1]]
  which_int <- names(covar)[sapply(as.list(covar), is.integer)]
  if(length(which_int) != 1)
    stop("Exactly one covariate, corresponding to the discrete time of the AR process,
         has to be provided as integer.")
  attr(object, "ARtime") <- which_int
  object
}



#' @export
#' @import mgcv nlme Matrix
#' @rdname nlme::corMatrix.corStruct
#'
corMatrix.corFunAR1 <- function(object, covariate = getCovariate(object),
                                covariance = TRUE, ...) {

  # get residuals
  Residuals <- attr(object, "residuals") # assigned by update.corDynamic_init / update.corSmooth
  if(is.null(Residuals))
    Residuals <- c(attr(object, "get_residuals")())

  # check whether models needs to be refit
  refit <- is.null(attr(object, "model"))
  if(!refit) {
    oldpars <- as.vector(attr(object, "model")$smooth[[1]]$sp)
    refit <- !(all.equal(oldpars, coef(object, unconstrained = FALSE)) == TRUE)
  }

  # name of the time variable of the AR process
  ARtime <- attr(object, "ARtime")
  covnames <- setdiff(all.vars(getCovariateFormula(object)), ARtime)
  covnames_ <- paste0(covnames, "_")

  if(refit) {
    grps <- getGroups(object)
    if(!is.list(Residuals)) {
      Residuals <- split(Residuals, grps)
    }

    ## first estimate lag 0 covariance as in corSmooth -------------------------

    # build covariance data (assuming vector covariate for now)
    make_lag0_data <- function(x, r) {
      if(nrow(x) < 2)
        return(NULL)
      idx <- combn(seq_len(nrow(x)), 2)
      d <- x[idx[1,], , drop = FALSE]
      d[covnames_] <- x[idx[2, ], , drop = FALSE]
      d$residuals2 <- combn(r, 2, prod)
      d$diagonal <- 0
      x[covnames_] <- x
      x$residuals2 <- r^2
      x$diagonal <- 1
      rbind(d, x)
    }

    covariate_comb <- Map(function(x, r) {
      art <- x[[ARtime]]
      x <- split(x[covnames], art)
      r <- split(r, art)
      covco <- Map(make_lag0_data, x, r)
      covco <- do.call(rbind, covco)
    }, covariate, Residuals)

    covariate_comb <- do.call(rbind, covariate_comb)

    # fit lag 0 covariance model
    args <- as.list(attr(object, "s_args"))
    args$xt <- as.list(args$xt)
    args$xt$absorb.cons <- FALSE
    kform <- as.formula(paste("residuals2 ~ 0 + s(",
                              paste(c(covnames, covnames_), collapse = ","),
                              ", bs = 'symm', xt = args$xt, k = args$k, m = args$m) + diagonal"))
    k <- gam(kform, data = covariate_comb,
             sp = head(coef(object, unconstrained = FALSE), length(covnames)))
    if(attr(object, "verbose"))
      plot(k, asp = 1, main = paste("Lag 0 Smoother penalty:", k$smooth[[1]]$sp))
    k <- force_non_negative(k)
    if(attr(object, "verbose"))
      cat("Eigenvalues:", attr(k, "eigen(coefMat)")$values)

    if(coef(k)["diagonal"] < 0)
      k$coefficients["diagonal"] <- 0 # ensure non-negative error variance
    if(attr(object, "verbose"))
      cat(" --- Noise variance:", k$coefficients["diagonal"], "\n")

    # store model object
    f <- attr(object, "lme_env")
    if(!is.null(f$lmeSt)) {
      attr(f$lmeSt$corStruct, "model") <- k
    }
  }

  # obtain fitted values --------------

  get_lag0_fit <- function(x) {
    idx <- seq_len(nrow(x))
    idx <- expand.grid(V1 = idx, V2 = idx)
    d <- cbind(x[idx$V1, , drop = FALSE],
               structure(x[idx$V2, , drop = FALSE], names = covnames_))
    d$diagonal <- as.numeric(idx$V1 == idx$V2)
    matrix(predict(
      if(refit) k else attr(object, "model"),
      d), nrow = nrow(x))
  }

  val0 <- lapply(covariate, function(x) {
    x <- split(x[covnames], x[[ARtime]])
    # return inner list of covariance matrices
    lapply(x, get_lag0_fit)
  })

  ## then estimate lag 1 covariance analogously -------------------------------

  make_lag1_data <- function(x1, x2, r1, r2) {
    dims <- c(nrow(x1), nrow(x2))
    idx <- list(seq_len(dims[1]), seq_len(dims[2]))
    idx <- expand.grid(V1 = idx[[1]],
                     V2 = idx[[2]])
    d <- cbind(x1[idx$V1, , drop = FALSE],
               structure(x2[idx$V2, , drop = FALSE], names = paste0(names(x2), "_")))

    d$residuals2 <- kronecker(r2, r1)
    attr(d, "dims") <- dims
    d
  }

  grps <- split(grps, grps) # for reordering

  covariate_comb <- Map(function(x, r, gr) {
    art <- x[[ARtime]]
    x$grps <- gr
    x <- split(x, art)
    r <- split(r, art)
    covco <- Map(make_lag1_data,
                 x1 = x[-length(x)], x2 = x[-1],
                 r1 = r[-length(x)], r2 = r[-1])
    dims <- lapply(covco, attr, "dims")
    covco <- do.call(rbind, covco)
    attr(covco, "dims") <- dims
    covco
  }, covariate, Residuals, grps)
  cov_dims <- lapply(covariate_comb, attr, "dims")
  covariate_comb <- do.call(rbind, covariate_comb)

  # fit lag 0 covariance model
  kform <- as.formula(paste("residuals2 ~ 0 + ti(",
                            paste(c(covnames, covnames_), collapse = ","),
                            ", bs = args$xt$bsmargin,
                            k = args$k, m = args$m, mc = c(FALSE, FALSE))"))
  k1 <- gam(kform, data = covariate_comb,
           sp = rep(tail(coef(object, unconstrained = FALSE), length(covnames)), 2))
  if(attr(object, "verbose"))
    plot(k1, asp = 1, main = paste("Lag 1 Smoother penalty:", k1$smooth[[1]]$sp))

  # store model object
  f <- attr(object, "lme_env")
  if(!is.null(f$lmeSt)) {
    attr(f$lmeSt$corStruct, "model_lag1") <- k1
  }


  # obtain fitted values ---------

  browser()

  covariate_comb$fitted.values <- k1$fitted.values
  val1 <- split(covariate_comb,
                paste(covariate_comb$grps, covariate_comb$grps_))
  val1 <- Map(function(x, dims) {
    x <- split(x, paste(x[[ARtime]], x[[paste(ARtime, "_")]]))
    Map(function(x, dims) matrix(x$residuals2, nrow = dims[1], ncol = dims[2]),
        x, dims)
  }, val1, cov_dims)


  if(!covariance) {
    for(i in seq_along(val)) {
      for(j in seq_along(val[[i]])) {
        val[[i]][[j]] <- cov2cor(val[[i]][[j]])
      }
    }
  }

  # compute factor
  e <- lapply(val, lapply, eigen, symmetric = TRUE)
  fac <- lapply(e, function(x) {
    facl <- lapply(x, function(x) 1/sqrt(x$values) * t(x$vectors))
    as.matrix(bdiag(facl))
  })
  attr(fac, "logDet") <- sum(log(unlist(lapply(e, lapply, `[[`, "values"))))
  attr(val, "factor") <- fac

  # partly extend covmats
  val <- lapply(val, function(x) as.matrix(bdiag(x)))

  val
}

