


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

  coefs <- attr(object, "coefficients")
  # check whether models needs to be refit
  refit <- is.null(coefs)
  if(!refit) {
    oldpars <- attr(coefs, "sp")
    refit <- !(all.equal(oldpars, coef(object, unconstrained = FALSE)) == TRUE)
  }

  if(refit) {
    grps <- getGroups(object)
    if(!is.list(Residuals)) {
      Residuals <- split(Residuals, grps)
    }

    coefs <- attr(object, "solvePLS")(
      sp = coef(object, unconstrained = FALSE),
      y = Residuals
      )

    # get positive definite part (in form of its eigen decomposition)
    ecoefs <- eigen(coefs, symmetric = T)
    pos <- ecoefs$values > 0
    ecoefs$values <- ecoefs$values[pos]
    ecoefs$vectors <- ecoefs$vectors[, pos, drop = FALSE]
    attr(ecoefs, "sp") <- attr(coefs, "sp")

    # store eigen decomposition of coefficients
    f <- attr(object, "coefficients")
    if(!is.null(f$lmeSt))
        attr(f$lmeSt$corStruct, "coefficients") <- ecoefs
    }

    # if(attr(object, "verbose"))
    #   image(coefs, asp = 1, main = paste("Smoother penalty:", attr(coefs, "sp")))

    if(attr(object, "verbose"))
      cat("Eigenvalues:", ecoefs$values)

    # obtain predictions
    X <- environment(attr(object, "solvePLS"))$X
    val <- lapply(X, function(x) {
      x <- sweep( x %*% ecoefs$vectors, 2, ecoefs$values, `/`)
      tcrossprod(x)
    })
    res2 <- unlist(Map(function(pred, res) {
      res^2 - diag(pred)
    }, val, Residuals))
    sigma2 <- mean(res2)

    browser()

    # if(coef(k)["diagonal"] < 0) {
    #   # do Manuel Pfeuffer's positivity trick
    #   # to ensure non-negative error variance
    #   thisdiag <- which(names(coef(k)) == "diagonal")
    #   sddiag <- sqrt(k$Vp[thisdiag, thisdiag])
    #   # set to mean of normal truncated at 0
    #   k$coefficients["diagonal"] <- k$coefficients["diagonal"] + 2*dnorm(0)*sddiag
    # }

    if(attr(object, "verbose"))
      cat(" --- Noise variance:", k$coefficients["diagonal"], "\n")


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
  S <- tensor.prod.penalties(rep(sm$S, 2))[[1]]

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

  # prepare fitting function
  attr(object, "solvePLS") <- function(sp, y) {
    stopifnot(is.list(y) & length(y) == length(X))
    # compute "Xy"
    XxX_YxY <- matrix(0, nrow = nrow(XxX_XxX))
    for(i in seq_along(X)) {
      X_Y <- crossprod(X[[i]], as.matrix(y[[i]]))
      XxX_YxY <- XxX_YxY + kronecker(X_Y, X_Y)
    }

    # solve PLS to get coefficients
    structure(matrix(
        sweep(L, 2, 1 + sp*eK$values, `/`) %*% crossprod(L, XxX_YxY),
      ncol = sqrt(ncol(L))), sp = sp) # store smoothing parameter used for fitting
  }

  object
}




