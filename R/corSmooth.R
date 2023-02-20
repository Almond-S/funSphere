


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
corSmooth <- function(value = 0, working_correlation = NULL,
                      form = formula(working_correlation), fixed = FALSE,
                      s_xt = list(), s_k = -1, s_m = NA, verbose = TRUE) {
  if (any(value < 0)) {
    stop("penalty parameter for covariance smoothing must be non-negative")
  }
  value <- log(value)

  attr(value, "s_args") <- list(xt = s_xt, k = s_k, m = s_m)
  attr(value, "verbose") <- verbose
  class(value) <- c("corSmooth", "corStruct")

  corDynamic(value, working_correlation, form = form, fixed = fixed)
}


#' @export
#' @import nlme
#' @rdname nlme::coef.corStruct
#'
coef.corSmooth <- function (object, unconstrained = TRUE, ...) {
  i_bims_coef.corSmooth <- 1
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

get_grps <- function(object) {
  groups <- getGroupsFormula(object$reStruct)

}


#' @export
#' @import mgcv nlme
#' @rdname nlme::corMatrix.corStruct
#'
corMatrix.corSmooth <- function(object, covariate = getCovariate(object),
                                initial.corStruct = NULL,
                                return.model = FALSE,
                                parent.frame.modelStruct = !return.model, ...) {

  # get residuals
  Residuals <- c(attr(object, "residuals")())
  # grps <- with(attr(object, "lme_env"), grps[revOrder, , drop = FALSE])
  grps <- getGroups(object)

  res <- split(Residuals, grps)

  # build covariance data (assuming vector covariate for now)
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

  if(return.model)
    return(k)

  val <- lapply(covariate, function(x) {
    d <- expand.grid(V1 = x, V2 = x)
    d$diagonal <- as.numeric(d$V1 == d$V2)
    matrix(predict(k, d), nrow = length(x))
  })

  for(i in seq_along(val)) {
    # diag(val[[i]]) <- diag(val[[i]]) + nugget
    val[[i]] <- cov2cor(val[[i]])
  }

  # compute factor
  e <- lapply(val, eigen, symmetric = TRUE)
  fac <- lapply(e, function(x) 1/sqrt(x$values) * t(x$vectors))
  attr(fac, "logDet") <- sum(log(unlist(lapply(e, `[[`, "values"))))
  attr(val, "factor") <- fac

  val
}

corFactor.corSmooth <- function(object, ...) {
  # so far only copy of corFactor.corStruct
  # setting corMatrix(..., parent.frame.update = FALSE)
  if (!is.null(aux <- attr(object, "factor"))) {
    return(aux)
  }
  corD <- Dim(object)
  val <- .C(corStruct_factList, as.double(unlist(
    corMatrix(object, parent.frame.modelStruct = FALSE))),
            as.integer(unlist(corD)), factor = double(corD[["sumLenSq"]]),
            logDet = double(1))[c("factor", "logDet")]
  lD <- val[["logDet"]]
  val <- val[["factor"]]
  attr(val, "logDet") <- lD
  val
}

