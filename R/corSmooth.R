


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
    f <- attr(object, "lme_env")
    if(!is.null(f$lmeSt))
        attr(f$lmeSt$corStruct, "coefficients") <- ecoefs

    if(attr(object, "verbose"))
      cat("Eigenvalues:", ecoefs$values)
  }


  if(refit | !covariance | ((corr | is.null(attr(object, "fac"))) &
       is.null(attr(object, "covariance"))) ) {

    if(!refit)
      ecoefs <- attr(object, "coefficients")

    # obtain predictions
    X <- environment(attr(object, "solvePLS"))$X

    val <- lapply(X, function(x) {
      x <- sweep( x %*% ecoefs$vectors, 2, sqrt(ecoefs$values), `*`)
      tcrossprod(x)
    })
    res2 <- unlist(Map(function(pred, res) {
      res^2 - diag(pred)
    }, val, Residuals))
    # residual noise variance
    sigma2 <- max(mean(res2), min(ecoefs$values)/2)

    val <- lapply(val, function(x) {
      diag(x) <- diag(x) + sigma2
      x})
    attr(val, "sigma2noise") <- sigma2

    if(attr(object, "verbose"))
      cat(" --- Noise variance:", sigma2, "\n")

    # store eigen decomposition of coefficients
    f <- attr(object, "lme_env")
    if(!is.null(f$lmeSt))
      attr(f$lmeSt$corStruct, "covariance") <- val

    if(!covariance) {
      for(i in seq_along(val)) {
        val[[i]] <- cov2cor(val[[i]])
      }
    }

      if(corr)
        return(val)

    } else {
      if(!is.null(attr(object, "fac")))
        return(attr(object, "fac"))

      # only remaining case: fac wanted but only covariance available
      val <- attr(objec, "covariance")
    }

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


#' @param X marginal design matrix of \code{kronecker(X, X)}.
#' @param W matrix containing diagonal of diagonal weight matrix of the kronecker design.
#' Default: zero weights on diagonal and elsewhere one.
#'
#' @export
get_XxXtXxX <- function(X, W = 1- diag(nrow = nrow(X)), RX = row_tensor_square(X)) {
  # matrix obtained via array model has to be reorganized:
  XxXtXxX_ <- array(crossprod(RX, W %*% RX), dim = rep(ncol(X), 4))
  matrix(aperm(XxXtXxX_, c(1,3,2,4)), ncol = ncol(RX))
}


#' Title
#'
#' @param X marginal designmatrix in \code{kronecker{X,X}}.
#' @param Y marginal response in \code{kronecker{Y,Y}} minus diagonal.
#'
#' @export
#' @import nlme sparseFLMM
get_XxXtYxY_noDiag <- function(X, Y, RX = row_tensor_square(X)) {
  XtY <- crossprod(X, Y)
  kronecker(XtY, XtY) - crossprod(RX, Y^2) # subtract diagonal
}

#' Title
#'
#' @param XtX inner product matrix of design matrix columns
#' @param S penalty matrix
#' @param XtX_trafo a basis transformation matrix for restricting to a subspace
#' basis \code{X %*% XtX_trafo}. Note that the argument \code{XtY} of the
#' returned fitting function will stay untransformed and also the coefficients
#' of the untransformed basis will be returned by it.
#'
#' @export
get_demmlerreinsch_solver <- function(XtX, S, X_trafo = NULL) {
  # first decomposition of design product
  eX <- if(is.null(X_trafo)) eigen(XtX) else
    eigen(crossprod(X_trafo, XtX) %*% X_trafo)
  # left cholesky-type factor
  L <- sweep(eX$vectors, 2, 1/sqrt(eX$values), `*`)
  # adjust penalty
  K <- if(is.null(X_trafo)) crossprod( L, S ) %*% L else
    crossprod( L, crossprod(X_trafo, S) %*% X_trafo ) %*% L

  # then decomposition of adjusted penalty
  eK <- eigen(K)
  # update side factor
  L <- L_right <- L %*% eK$vectors
  # take trafo into prediction matrix
  if(!is.null(X_trafo)) {
    L_right <- X_trafo %*% L
  }

  # return fitting function
  function(sp, # smoothing parameter
           XtY) {
    # solve PLS to get coefficients
    coefs <- structure(
      sweep(L, 2, 1 + sp*eK$values, `/`) %*% crossprod(L_right, XtY),
      ncol = sqrt(ncol(L)), sp = sp) # store smoothing parameter used for fitting
    if(is.null(X_trafo))
      coefs else X_trafo %*% coefs
  }
}


#' @import nlme sparseFLMM
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
  nrowX <- sapply(X, nrow)
  ncolX <- sqrt(ncol(S))

  XxXtXxX <- array(0, dim = dim(S))
  RX <- list()
  # omit X with only one observation
  for(i in which(nrowX>1)) {
    # Note that diagonals are omitted in get get_XxXtXxX() for covariance estimation
    RX[[i]] <- row_tensor_square(X[[i]])
    XxXtXxX <- XxXtXxX + get_XxXtXxX(X[[i]], RX = RX[[i]])
  }

  # use symmetry
  Z <- make_summation_matrix(ncolX)
  # use Demmler-Reinsch type form to speed up computation
  DRsolve <- get_demmlerreinsch_solver(XxXtXxX, S, X_trafo = Z)

  # prepare fitting function
  attr(object, "solvePLS") <- function(sp, y) {
    stopifnot(is.list(y) & length(y) == length(X))
    # compute "Xy"
    XxXtYxY <- matrix(0, nrow = nrow(XxXtXxX))
    # again omit X with only one observation
    for(i in which(nrowX>1)) {
      XxXtYxY <- XxXtYxY + get_XxXtYxY_noDiag(X[[i]], y[[i]], RX = RX[[i]])
    }
    # solve PLS to get coefficient matrix
    matrix(DRsolve(sp, XxXtYxY), ncol = ncolX, nrow = ncolX)
  }
  object
}




