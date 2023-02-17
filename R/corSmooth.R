


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
#' @import mgcv nlme
#' @export
#'
corSmooth <- function(value = 0, form = ~1, fixed = FALSE, G, s_xt = list(), s_k = -1, s_m = NA, verbose = TRUE) {
  if (any(value < 0)) {
    stop("penalty parameter for covariance smoothing must be non-negative")
  }
  value <- log(value)

  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  attr(value, "G") <- G
  attr(value, "s_args") <- list(xt = s_xt, k = s_k, m = s_m)
  attr(value, "verbose") <- verbose
  class(value) <- c("corSmooth", "corStruct")
  value
}




#' @export
#' @import nlme
#' @rdname nlme::Initialize
#'
Initialize.corSmooth <- function(object, data, ...) {
  object <- NextMethod()
  object
}

#' @export
#' @import nlme
#' @rdname nlme::coef.corStruct
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
  aux <- exp(as.vector(object))
  len_sp_mean <- length(attr(object, "G")$sp)
  aux
}


#' Extract gam smoothing paramter estimates from gamm object
#'
#' Internal helper function for \code{corMatrix.corSmooth}
#'
#' @param var.param
#'
#' @return vector of smoothing parameters for gam object.
#'
get_sp <- function(var.param) {
  1/mgcv:::notExp2(var.param)
}


#' @export
#' @import mgcv nlme
#' @rdname nlme::corMatrix.corStruct
#'
corMatrix.corSmooth <- function(object, covariate = getCovariate(object), ...) {
  # as in e.g. corMatrix.AR1:
  # corD <- Dim(object, if (is.list(covariate)) {
  #   if (is.null(names(covariate)))
  #     names(covariate) <- seq_along(covariate)
  #   rep(names(covariate), lengths(covariate))
  # }
  # else rep(1, length(covariate)))

  # obtain residuals
  sp <- NULL
  for(i in 2:4) {
    f <- parent.frame(i)
    if(is.list(f$object)) {
      sp <- f$object$sp
      if(!is.null(sp))
        break
      if(!is.null(f$object$reStruct)) {
        sp <- get_sp(coef(f$object$reStruct))
        break
      }
    }
  }

  if(is.null(sp) || is.na(sp)) {
    warning("Current penalty parameter not found - going back to default specified in G.")
    m <- gam(G = attr(object, "G"))
  } else {
    m <- gam(G = attr(object, "G"), sp = sp)
  }
  if(attr(object, "verbose")) {
    plot(m, main = paste("Smoother penalty:", sp))
  }

  # build covariance data (assuming vector covariate for now)
  grps <- getGroups(object)
  res <- split(m$residuals, grps)
  covariate_comb <- Map(function(x, r) {
    if(length(x) < 2)
      return(NULL)
    d <- as.data.frame(t(combn(x, 2)))
    d$residuals2 <- combn(r, 2, prod)
    d$diagonal <- 0
    d <- rbind(d, data.frame(V1 = x, V2 = x, residuals2 = r^2, diagonal = 1))
    d
  }, covariate, res)
  covariate_comb <- do.call(rbind, covariate_comb)

  # fit covariate model
  args <- as.list(attr(object, "s_args"))
  args$xt <- as.list(args$xt)
  args$xt$absorb.cons <- FALSE
  k <- gam(residuals2 ~ 0 +
             s(V1, V2, bs = "symm", xt = args$xt, k = args$k, m = args$m) + diagonal,
           data = covariate_comb,
           sp = coef(object, unconstrained = FALSE))
  if(attr(object, "verbose"))
    plot(k, asp = 1, main = paste("Smoother penalty:", k$sp))
  k <- force_non_negative(k)
  if(attr(object, "verbose"))
    cat("Eigenvalues:", attr(k, "eigen(coefMat)")$values)

  if(coef(k)["diagonal"] < 0)
    k$coefficients["diagonal"] <- 0 # ensure non-negative error variance
  if(attr(object, "verbose"))
    cat("Noise variance:", k$coefficients["diagonal"])

  attr(val, "cov_model") <- k

  val <- lapply(covariate, function(x) {
    d <- expand.grid(V1 = x, V2 = x)
    d$diagonal <- as.numeric(d$V1 == d$V2)
    matrix(predict(k, d), nrow = length(x))
  })

  ## earlier version mostly without prediction on new data
  # # obtain covariance matrices
  # grps_ <- ordered(unlist(split(grps, grps)), levels = levels(grps))
  # val <- split(k$fitted.values, grps_)
  # val <- Map(function(r, L) {
  #   val <- matrix(0, nrow = L, ncol = L)
  #   idx <- combn(seq_len(L), 2, function(x) x[2] + L*(x[1]-1) )
  #   idx_ <- combn(seq_len(L), 2, function(x) x[1] + L*(x[2]-1) )
  #   val[idx] <- val[idx_] <- r
  # }, val, lengths(covariate))
  #
  # # obtain/estimate variances
  # vc <- unsplit(covariate, grps)
  # diag(val) <- predict(k, newdata = data.frame(V1 = vc, V2 = vc))

  ## no positive definiteness guarantee for the moment
  # # restrict to positive definite part
  # e <- eigen(val, symmetric = TRUE)
  # npos <- sum(e$values > 0)
  # val <- crossprod(sqrt(head(e$values, npos)) * t(e$vectors[, seq_len(npos)]))

  # assume constant error variance
  # nugget <- Map(function(r,v) r^2 - diag(v), res, val)
  # nugget <- unlist(nugget)
  # nugget <- mean(nugget[nugget>0])
  for(i in seq_along(val)) {
    # diag(val[[i]]) <- diag(val[[i]]) + nugget
    val[[i]] <- cov2cor(val[[i]])
  }
  # if(any(abs(sapply(val, determinant)) < 1e-15))
  #   browser()

  val
}
