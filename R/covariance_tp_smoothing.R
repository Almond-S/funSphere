
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

#' @param X marginal designmatrix in \code{kronecker{X,X}}.
#' @param Y marginal response in \code{kronecker{Y,Y}} minus diagonal.
#'
#' @export
#' @import nlme sparseFLMM
get_XxXtYxY_noDiag <- function(X, Y, RX = row_tensor_square(X)) {
  XtY <- crossprod(X, Y)
  kronecker(XtY, XtY) - crossprod(RX, Y^2) # subtract diagonal
}

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
  eX <- if(is.null(X_trafo)) eigen(XtX, symmetric = TRUE) else
    eigen(crossprod(X_trafo, XtX) %*% X_trafo, symmetric = TRUE)
  # left cholesky-type factor
  L <- if(all(eX$values > 0))
          sweep(eX$vectors, 2, 1/sqrt(eX$values), `*`) else {
            k <- sum(eX$values>0)
            message(cat("Not all eigenvalues of XtX are positive.",
                        paste("Restricting to the k =", k, "basis functions where they are.")))
            sweep(eX$vectors[, 1:k], 2, 1/sqrt(eX$values[1:k]), `*`)
          }

  # adjust penalty
  K <- if(is.null(X_trafo)) crossprod( L, S ) %*% L else
    crossprod( L, crossprod(X_trafo, S) %*% X_trafo ) %*% L

  # then decomposition of adjusted penalty
  eK <- eigen(K, symmetric = TRUE)
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


#' @param X1,X2 marginal design matrices of \code{kronecker(X2, X1)}.
#' @param w a weight vector of length \code{nrow(X1)*nrow(X2)}.
#' @param RX1,RX2 row tensor products if design matrices
#' (only available as arguments to allow precomputation).
#'
#' @export
get_X2xX1tX2xX1_weighted <- function(X1, X2 = X1, w = rep(1, nrow(X1)*nrow(X2)),
                                     RX1 = row_tensor_square(X1),
                                     RX2 = row_tensor_square(X2)) {
  # matrix obtained via array model has to be reorganized:
  X2xX1tX2xX1_ <- array(crossprod(RX1, matrix(w, ncol = nrow(RX2)) %*% RX2),
                        dim = rep(c(ncol(X1), ncol(X2)), each = 2))
  matrix(aperm(X2xX1tX2xX1_, c(1,3,2,4)), ncol = ncol(X1)*ncol(X2))
}


#' @param X1,X2 marginal design matrices of \code{kronecker(X2, X1)}.
#' (only available as arguments to allow precomputation).
#'
#' @export
get_X2xX1tX2xX1 <- function(X1, X2 = X1) {
  kronecker(crossprod(X2, X2), crossprod(X1, X1))
}


#' @param X1,X2 marginal design matrices in \code{kronecker{X2,X1}}.
#' @param Y1,Y2 marginal responses in \code{kronecker{Y2,Y1}}.
#'
#' @export
#' @import nlme sparseFLMM
get_X2xX1tY2xY1 <- function(X1, X2, Y1, Y2) {
  X1tY1 <- crossprod(X1, Y1); X2tY2 <- crossprod(X2, Y2)
  kronecker(X2tY2, X1tY1) # subtract diagonal
}


#' Set up smooth covariance fitting
#'
#' @param X list of according marginal design matrices
#' @param S list with penalty matrix (or matrices) of marginal design matrices
#'
#' @return symmetric covariance fitting function returning coefficient matrix of
#' estimated tensor product smooth and taking the arguments
#' \code{y}: a list of zero-mean responses (typically model residuals) and
#' \code{sp}: a numeric (vector) specifying the smoothing parameter(s) for the penalty.
#' @import sparseFLMM
#' @export
cov_symm_setup <- function(X, S) {
  # get tensor product penalty matrix
  S <- tensor.prod.penalties(rep(S, 2))[[1]]

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

  # use Demmler-Reinsch type form to speed up computation in case of repeated evaluation
  DRsolve <- get_demmlerreinsch_solver(XxXtXxX, S, X_trafo = Z)

  # return fitting function
  function(y, sp, return.fun = FALSE) {
    stopifnot(is.list(y) & length(y) == length(X))
    # compute "Xy"
    XxXtYxY <- matrix(0, nrow = nrow(XxXtXxX))
    # again omit X with only one observation
    for(i in which(nrowX>1)) {
      XxXtYxY <- XxXtYxY + get_XxXtYxY_noDiag(X[[i]], y[[i]], RX = RX[[i]])
    }
    # solve PLS to get coefficient matrix
    fit <- function(sp) matrix(DRsolve(sp, XxXtYxY), ncol = ncolX, nrow = ncolX)
    if(return.fun)
      return(fit)
    fit(sp)
  }
}

#' @export
predict_square_smooths <- function(thetas, smooth_obj, newdata,
                                   decompose = FALSE,
                                   Gramian = NULL,
                                   Gramian_chol = NULL,
                                   Gramian_chol_inv = NULL) {
  if(!is.list(thetas))
    thetas <- list(pred = thetas)
  X <- Predict.matrix(smooth_obj, data = newdata)

  lapply(thetas, predict_square_smooth, X,
         decompose = decompose, Gramian = Gramian,
         Gramian_chol = Gramian_chol,
         Gramian_chol_inv = Gramian_chol_inv)
}


#' Predict a single-tensor-product smooth
#'
#' @param theta an object representing the estimated coefficients
#' @param smooth_obj a gam smoother object with a \code{Predict.matrix} method
#' @param newdata data for prediction
#' @param decompose logical, should prediction be returned in form of its SVD?
#' @param Gramian an optional Gramian providing the Gam-matrix of the splines for SVD.
#' For NULL, the default, the SVD is conducted with respect to an inner product in
#' which the basis functions are orthonormal.
#'
#' @export
predict_square_smooth <- function(theta, X, decompose = FALSE, Gramian = NULL,
                                  Gramian_chol = NULL,
                                  Gramian_chol_inv = NULL, ...) {
  UseMethod("predict_square_smooth")
}

#' @export
predict_square_smooth.matrix <- function(theta, X, decompose = FALSE, symmetric = FALSE,
                                         Gramian = NULL,
                                         Gramian_chol = NULL,
                                         Gramian_chol_inv = NULL) {

  if(!decompose) return(
    wtcrossprod(X, w = theta)
  )
  # alternatively make SVD first
  if(!is.null(Gramian)) {
    U <- if(is.null(Gramian_chol)) chol(Gramian) else Gramian_chol
    U_ <- if(is.null(Gramian_chol_inv)) solve(U) else Gramian_chol_inv
    theta <- wtcrossprod(U, w = theta)
    X <- X %*% U_
  }

  if(symmetric) {
    e <- as.list(eigen(theta, symmetric = TRUE))
    names(e) <- c("d", "u")
    e$u <- X %*% e$u
    return(e)
  }

  SVD <- svd(theta)
  SVD[c("u", "v")] <- lapply(SVD[c("u", "v")], function(y) X %*% y)

  SVD
}

#' @export
predict_square_smooth.eigen <- function(theta, X, decompose = FALSE,
                                        Gramian = NULL,
                                        Gramian_chol = NULL,
                                        Gramian_chol_inv = NULL) {

  if(!decompose) {
    X <- X %*% theta$vectors
    return( wtcrossprod(X, w = theta$values) )
    }
  # alternatively reduce to matrix case
  theta <- wtcrossprod(theta$vectors, w = theta$values)
  predict_square_smooth.matrix(theta, X, decompose = decompose,
                               Gramian = Gramian,
                               Gramian_chol = Gramian_chol,
                               Gramian_chol_inv = Gramian_chol_inv)
  # if(!is.null(Gramian)) {
  #   Gramian <- wcrossprod(theta$vectors, w = Gramian)
  #   U <- chol(Gramian)
  #   U_ <- solve(U)
  #   theta$values <- wtcrossprod(U, w = theta$values)
  #   X <- X %*% U_
  # }
  # SVD <- if(is.matrix(theta$values))
  #   svd(theta$values,
  #       nu = nrow(theta$values),
  #       nv = nrow(theta$values)) else list(d = theta$values,
  #                               u = theta$vectors, v = theta$vectors)
  # SVD[c("u", "v")] <- lapply(SVD[c("u", "v")], function(y) tcrossprod(X, y))
  #
  # SVD
}


#' Set up smooth covariance fitting for two cross-covariances
#'
#' @param X1,X2 lists of according marginal design matrices
#' @param S1,S2 lists with penalty matrix (or matrices) of marginal design matrices
#'
#' @return symmetric covariance fitting function returning coefficient matrix of
#' estimated tensor product smooth and taking the arguments
#' \code{y1,y2}: lists of zero-mean responses (typically model residuals) and
#' \code{sp1,sp2}: numeric (vectors) specifying the smoothing parameter(s) for the penalty.
#' @import sparseFLMM
#' @export
cov_cross_setup <- function(X1, X2 = X1, S1, S2 = S1) {
  stopifnot(length(X1) == length(X2))

  # get tensor product penalty matrices
  S <- tensor.prod.penalties(c(S1,S2))[[1]]

  nrowX1 <- sapply(X1, nrow); nrowX2 <- sapply(X2, nrow)
  ncolX1 <- ncol(S1[[1]]); ncolX2 <- ncol(S2[[1]])

  X2xX1tX2xX1 <- array(0, dim = dim(S))
  for(i in seq_along(X1)) {
    X2xX1tX2xX1 <- X2xX1tX2xX1 + get_X2xX1tX2xX1(X1[[i]], X2[[i]])
  }

  DRsolve2 <- get_demmlerreinsch_solver(X2xX1tX2xX1, S)

  # return fitting function
  function(y1, y2 = y1, sp, return.fun = FALSE) { # U is a basis trafo matrix
    stopifnot(is.list(y1) & length(y1) == length(X1))
    stopifnot(is.list(y2) & length(y2) == length(X2))
    # compute "Xy"
    X2xX1tY2xY1 <- matrix(0, nrow = nrow(X2xX1tX2xX1))
    # again omit X with only one observation
    for(i in seq_along(y1)) {
      X2xX1tY2xY1 <- X2xX1tY2xY1 + get_X2xX1tY2xY1(
            X1[[i]], X2[[i]], y1[[i]], y2[[i]])
      }
    # solve PLS to get coefficient matrix
    fit <- function(sp) matrix(DRsolve2(sp, X2xX1tY2xY1), ncol = ncolX2, nrow = ncolX1)
    if(return.fun)
      return(fit)
    fit(sp)
  }
}
