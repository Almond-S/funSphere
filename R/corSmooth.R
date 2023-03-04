


#' Correlation structure based on covariance smoothing
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
#' @param G an object of class \code{gam.prefit} as returned by \code{gam(..., fit = FALSE)}
#' containing the mean model structure.
#'
#' @return an object of class \code{corSmooth}, representing an covariance smoother autocorrelation structure.
#' @import mgcv nlme Matrix
#' @export
#'
corSmooth <- function(value = 0, form = ~1, fixed = FALSE,
                      working_correlation = corCAR1,
                      working_control = list(), verbose = FALSE) {
  if (any(value < 0)) {
    stop("penalty parameter for covariance smoothing must be non-negative")
  }
  value <- notLog2(value)

  attr(value, "verbose") <- verbose
  class(value) <- c("corSmooth", "corStruct")

  # extract suitable formula
  gf <- getGroupsFormula(form)
  af <- getCovariateFormula(form)
  cf <- interpret.gam(af); cf <- cf$pred.formula
  attr(value, "gam_formula") <- af
  # update form
  i <- if(length(form) == 2) 2 else 3
  form[[i]][[2]] <- cf[[2]]

  corDynamic(value, working_correlation = working_correlation,
             working_control = working_control,
             form = form, fixed = fixed)
}


#' @export
#' @import nlme
# #' @rdname nlme::coef.corStruct
#'
coef.corSmooth <- function (object, unconstrained = TRUE, ...) {
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


row_tensor_square <- function(x) x[, rep(1:ncol(x), each = ncol(x)), drop = FALSE] *
  x[, rep(1:ncol(x), ncol(x)), drop = FALSE]

#' @export
#' @import mgcv nlme MASS
# #' @rdname nlme::corMatrix.corStruct
#'
corMatrix.corSmooth <- function(object, covariate = getCovariate(object), corr = TRUE,
                                covariance = TRUE, ...) {

  # get residuals
  Residuals <- attr(object, "residuals") # assigned by update.corDynamic_init / update.corSmooth
  if(is.null(Residuals))
    Residuals <- c(attr(object, "get_residuals")())

  # check whether models needs to be refit
  refit <- TRUE # TODO: need to update refit condition
  if(!refit) {
    oldpars <- as.vector(attr(object, "model")$smooth[[1]]$sp)
    refit <- !(all.equal(oldpars, coef(object, unconstrained = FALSE)) == TRUE)
  }

  if(refit) {
    grps <- getGroups(object)
    if(!is.list(Residuals)) {
      Residuals <- split(Residuals, grps)
    }
    idx <- split(seq_along(grps), grps)

    # fit covariance via linear array model (Currie et al. 2006)
    X <- lapply(idx, function(id) attr(object, "smooth")$X[id, , drop = FALSE])
    S <- attr(object, "penalty")[[1]]
    # TODO: diagonals still need to be subtracted
    XxX_YxY <- matrix(0, nrow = nrow(XxX_XxX))
    # TODO: exclude design matrices with only one row when removing diagonals?
    for(i in seq_along(X)) {
      X_Y <- crossprod(X[[i]], as.matrix(Residuals[[i]]))
      XxX_YxY <- XxX_YxY + kronecker(X_Y, X_Y)
    }

    # compute cofficients
    coefs <- solve( attr(object, "XxX_XxX") + coef(object, unconstrained = FALSE) * S, c(XxX_YxY))

    browser()

    # build covariance data (assuming vector covariate for now)
    covariate_comb <- Map(function(x, r) {
      if(nrow(x) < 2)
        return(NULL)
      idx <- combn(seq_len(nrow(x)), 2)
      d <- x[idx[1,], , drop = FALSE]
      d[paste0(names(x), "_")] <- x[idx[2, ], , drop = FALSE]
      d$residuals2 <- combn(r, 2, prod)
      d$diagonal <- 0
      x[paste0(names(x), "_")] <- x
      x$residuals2 <- r^2
      x$diagonal <- 1
      rbind(d, x)
    }, covariate, Residuals)
    covariate_comb <- do.call(rbind, covariate_comb)

    # fit covariate model
    args <- as.list(attr(object, "s_args"))
    args$xt <- as.list(args$xt)
    args$xt$absorb.cons <- FALSE
    kform <- as.formula(paste("residuals2 ~ 0 + s(",
                              paste(names(covariate[[1]]), collapse = ","), ",",
                              paste(paste0(names(covariate[[1]]), sep = "_"), collapse = ","),
                              ", bs = 'symm', xt = args$xt, k = args$k, m = args$m) + diagonal"))
    if(is.null(attr(object, "model_prefit"))) {
      k_pre <- gam(kform, data = covariate_comb, fit = FALSE)
      # store prefit object
      f <- attr(object, "lme_env")
      if(!is.null(f$lmeSt))
        attr(f$lmeSt$corStruct, "model_prefit") <- k_pre
    } else {
      k_pre <- attr(object, "model_prefit")
      # update response
      k_pre$y <- covariate_comb$residuals2
    }

    # fit model
    k <- gam(G = k_pre, sp = coef(object, unconstrained = FALSE) )

    if(attr(object, "verbose"))
      plot(k, asp = 1, main = paste("Smoother penalty:", k$full.sp))
    k <- force_non_negative(k)
    if(attr(object, "verbose"))
      cat("Eigenvalues:", attr(k, "eigen(coefMat)")$values)

    if(coef(k)["diagonal"] < 0) {
      # do Manuel Pfeuffer's positivity trick
      # to ensure non-negative error variance
      thisdiag <- which(names(coef(k)) == "diagonal")
      sddiag <- sqrt(k$Vp[thisdiag, thisdiag])
      # set to mean of normal truncated at 0
      k$coefficients["diagonal"] <- k$coefficients["diagonal"] + 2*dnorm(0)*sddiag
    }

    if(attr(object, "verbose"))
      cat(" --- Noise variance:", k$coefficients["diagonal"], "\n")

    # store model object
    f <- attr(object, "lme_env")
    if(!is.null(f$lmeSt)) {
      attr(f$lmeSt$corStruct, "model") <- k
    }
  }

  # get appropriate marginal bases
  if(is.null(attr(object, "marginalDesign"))) {
    marginalDesign <- lapply(covariate, Predict.matrix, object = k$smooth[[1]]$margin[[1]])
    # store marginal design matrix
    f <- attr(object, "lme_env")
    if(!is.null(f$lmeSt))
      attr(f$lmeSt$corStruct, "marginalDesign") <- marginalDesign
  } else {
    marginalDesign <- attr(object, "marginalDesign")
  }

  val <- lapply(marginalDesign, function(X) {
    # extract design and coefficient matrices
    coefs <- if(refit) k$coefficients else attr(object, "model")$coefficients
    Z <- if(refit) k$smooth[[1]]$Z else attr(object, "model")$smooth[[1]]$Z
    coefs <- list(smooth = Z%*%coefs[-1], nugget = coefs[1])

    # return prediction
    pred <- X %*%
      tcrossprod( matrix(coefs$smooth, ncol = sqrt(length(coefs$smooth))),
                  X)
    diag(pred) <- diag(pred) + coefs$nugget
    pred
  })

  if(!covariance) {
    for(i in seq_along(val)) {
      val[[i]] <- cov2cor(val[[i]])
    }
  }

  if(corr)
    return(val)

  # otherwise compute factor:
  e <- lapply(val, eigen, symmetric = TRUE)
  # fac <- unlist(lapply(e, function(x) c(1/sqrt(x$values) * t(x$vectors))))
  fac <- lapply(e, function(x) 1/sqrt(x$values) * t(x$vectors))
  # log determinant of factor
  lD <- -1/2*sum(log(unlist(lapply(e, `[[`, "values"))))
  attr(fac, "logDet") <- lD

  fac
}



#' @import nlme
#' @export
corFactor.corSmooth <- function(object, ...) {
  if(!is.null(aux <- attr(object, "factor"))) {
    return(aux)
  }
  corMatrix(object, ..., corr = FALSE)
}

#' @import nlme
#' @export
Initialize.corSmooth <- function(object, data, ...) {
  object <- NextMethod()

  # construct smoother
  frm <- attr(object, "gam_formula")
  sm <- eval(frm[[length(frm)]])
  sm <- smooth.construct(sm, data, attr(object, "knots"))

  attr(object, "smooth") <- sm
  pls <- list()
  attr(object, "penalty") <- S <- tensor.prod.penalties(rep(sm$S, 2))[[1]]

  grps <- getGroups(object)
  idx <- split(seq_along(grps), grps)

  # compute design matrix innter product using linear array model (Currie et al. 2006)
  X <- lapply(idx, function(id) sm$X[id, , drop = FALSE])
  # TODO: diagonals still need to be subtracted
  XxX_XxX <- array(0, dim = dim(S))
  # TODO: exclude design matrices with only one row when removing diagonals?
  for(i in seq_along(X)) {
    RT <- row_tensor_square(X[[i]])
    ONE <- matrix(1, nrow = nrow(X[[i]]), ncol = nrow(X[[i]]))
    # matrix obtained via array model has to be reorganized:
    XxX_XxX_ <- array(crossprod(RT, ONE %*% RT), dim = rep(ncol(X[[i]]), 4))
    XxX_XxX <- XxX_XxX + matrix(aperm(XxX_XxX_, c(1,3,2,4)), ncol = ncol(RT))
  }

  # use Demmler-Reinsch type form to speed up computation

  # first decomposition of design product
  eX <- eigen(XxX_XxX)
  # left cholesky-type factor
  L <- sweep(eX$vectors, 2, 1/sqrt(eX$values), `*`)
  # adjust penalty
  K <- crossprod( L, S ) %*% L

  # then decomposition of adjusted penalty
  eK <- eigen(K)
  # update side factor
  L <- L %*% eK$vectors

  attr(object, "demmler_reinsch") <- function(sp, Xy) {
    sweep(L, 2, 1 + sp*eK$values, `/`) %*% crossprod(L, Xy)
  }

browser()

  attr(object, "XxX_XxX") <- XxX_XxX


  object
}




