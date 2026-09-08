
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
#' @param ridge a numerical safeguard for ill-conditioned designs: this
#' multiple of the mean diagonal of the (transformed) \code{XtX} is added to
#' its diagonal before the decomposition, so that directions the design
#' barely identifies are damped rather than amplified by \code{1/sqrt(d)}.
#' The default \code{0} is the unregularised solver; \code{1e-8} is what
#' \code{eigenfun::eigenfun_covsmooth()} uses on its pencil and was needed
#' for a 118-function soap-film basis (condition number 1e7 of the marginal
#' design, 1e14 of the pair design), where the unregularised solver returned
#' a surface with a relative error of 0.9 at every smoothing parameter.
#'
#' @export
get_demmlerreinsch_solver <- function(XtX, S, X_trafo = NULL, ridge = 0) {
  # X_trafo may be sparse (a Matrix): the two reductions below are then cheap
  # and the returned solver keeps the transformation rather than the dense
  # X_trafo %*% L it used to precompute -- at k = 118 basis functions a
  # 13924 x 7021 matrix (780 MB) and a 1.4 TFlop product per call of this
  # function, i.e. per cross-validation fold
  reduce <- function(M) if(is.null(X_trafo)) M else
    as.matrix(crossprod(X_trafo, M) %*% X_trafo)
  # first decomposition of design product
  XtX <- reduce(XtX)
  if(ridge > 0) diag(XtX) <- diag(XtX) + ridge * mean(diag(XtX))
  eX <- eigen(XtX, symmetric = TRUE)
  rm(XtX)
  # left cholesky-type factor
  L <- if(all(eX$values > 0))
          sweep(eX$vectors, 2, 1/sqrt(eX$values), `*`) else {
            k <- sum(eX$values>0)
            message(cat("Not all eigenvalues of XtX are positive.",
                        paste("Restricting to the k =", k, "basis functions where they are.")))
            sweep(eX$vectors[, 1:k], 2, 1/sqrt(eX$values[1:k]), `*`)
          }
  rm(eX)

  # adjust penalty
  K <- crossprod( L, reduce(S) ) %*% L

  # then decomposition of adjusted penalty
  eK <- eigen(K, symmetric = TRUE)
  rm(K)
  # update side factor
  L <- L %*% eK$vectors
  d <- eK$values
  rm(eK)
  ncolL <- ncol(L)

  # return fitting function; its environment holds L, d and X_trafo only
  solver <- function(sp, # smoothing parameter
                     XtY) {
    # solve PLS to get coefficients: L diag(1 / (1 + sp d)) L' X_trafo' XtY,
    # mapped back to the untransformed basis
    if(!is.null(X_trafo)) XtY <- as.matrix(crossprod(X_trafo, XtY))
    coefs <- L %*% (crossprod(L, XtY) / (1 + sp * d))
    if(!is.null(X_trafo)) coefs <- as.matrix(X_trafo %*% coefs)
    structure(coefs, ncol = sqrt(ncolL), sp = sp) # store smoothing parameter used for fitting
  }
  environment(solver) <- list2env(list(L = L, d = d, X_trafo = X_trafo, ncolL = ncolL),
                                  parent = environment(get_demmlerreinsch_solver))
  solver
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
#' @param ridge passed on to \code{\link{get_demmlerreinsch_solver}}: a relative
#' ridge on the pair design product, \code{0} by default.
#'
#' @return symmetric covariance fitting function returning coefficient matrix of
#' estimated tensor product smooth and taking the arguments
#' \code{y}: a list of zero-mean responses (typically model residuals) and
#' \code{sp}: a numeric (vector) specifying the smoothing parameter(s) for the penalty.
#' @import sparseFLMM
#' @export
cov_symm_setup <- function(X, S, ridge = 0) {
  # get tensor product penalty matrix: S (x) I, the first of the pair
  # tensor.prod.penalties(rep(S, 2)) forms -- kept sparse, it is k^2 x k^2
  S <- kronecker(Matrix(S[[1]], sparse = TRUE), Diagonal(ncol(S[[1]])))

  nrowX <- sapply(X, nrow)
  ncolX <- ncol(X[[1]])
  # omit X with only one observation
  keep <- which(nrowX > 1)
  # Note that diagonals are omitted in get get_XxXtXxX() for covariance estimation
  RX <- vector("list", length(X))
  RX[keep] <- lapply(X[keep], row_tensor_square)

  # sum_i RX_i' W_i RX_i with W_i = 1 - I, assembled in ONE crossprod of the
  # stacked row tensor products: (W_i RX_i)[j, ] = colSums(RX_i) - RX_i[j, ].
  # The per-curve version formed a k^2 x k^2 matrix (1.5 GB at k = 118) and a
  # permuted copy of it for every curve; the arithmetic is the same.
  RXs <- do.call(rbind, RX[keep])
  gid <- rep(keep, nrowX[keep])
  WRX <- rowsum(RXs, gid, reorder = FALSE)[match(gid, keep), , drop = FALSE] - RXs
  XxXtXxX <- crossprod(RXs, WRX)
  rm(RXs, WRX, gid)
  # matrix obtained via array model has to be reorganized (in place)
  dim(XxXtXxX) <- rep(ncolX, 4)
  XxXtXxX <- aperm(XxXtXxX, c(1,3,2,4))
  dim(XxXtXxX) <- rep(ncolX^2, 2)

  # use symmetry
  Z <- make_summation_matrix(ncolX, sparse = TRUE)

  # use Demmler-Reinsch type form to speed up computation in case of repeated evaluation
  DRsolve <- get_demmlerreinsch_solver(XxXtXxX, S, X_trafo = Z, ridge = ridge)
  q0 <- nrow(XxXtXxX)
  rm(XxXtXxX, S, Z)

  # return fitting function
  function(y, sp, return.fun = FALSE) {
    stopifnot(is.list(y) & length(y) == length(X))
    # compute "Xy"
    XxXtYxY <- matrix(0, nrow = q0)
    # again omit X with only one observation
    for(i in keep) {
      XxXtYxY <- XxXtYxY + get_XxXtYxY_noDiag(X[[i]], y[[i]], RX = RX[[i]])
    }
    # solve PLS to get coefficient matrix. The function's environment holds
    # the solver and the right-hand side only, so that a cross-validation
    # keeping one per fold does not keep every fold's design products too
    fit <- function(sp) matrix(DRsolve(sp, XxXtYxY), ncol = ncolX, nrow = ncolX)
    environment(fit) <- list2env(list(DRsolve = DRsolve, XxXtYxY = XxXtYxY, ncolX = ncolX),
                                 parent = environment(cov_symm_setup))
    if(return.fun)
      return(fit)
    fit(sp)
  }
}

#' Choose the smoothing parameter of a symmetric covariance smooth by k-fold
#' cross-validation
#'
#' Splits the *curves* into \code{kfolds} folds, refits
#' \code{\link{cov_symm_setup}} on each training set, and picks the smoothing
#' parameter minimising the held-out squared prediction error of the
#' off-diagonal products \eqn{y_{ij} y_{ik}}, \eqn{j \neq k}. The diagonal is
#' excluded on both sides, matching the criterion \code{cov_symm_setup()} fits.
#'
#' Splitting by curve is what makes the criterion honest: the products within a
#' curve are exactly what the fit interpolates, so a split that leaves any of a
#' test curve's observations in the training set measures training error.
#'
#' Each curve contributes the *sum* over its held-out pairs, so curves observed
#' more often carry more weight -- they also carry more information. Curves with
#' fewer than two observations contribute no pair and are dropped from the test
#' folds, as they are from the fit.
#'
#' @param X list of per-curve marginal design matrices, as for
#' \code{\link{cov_symm_setup}}.
#' @param S list with the penalty matrix (or matrices) of the marginal basis.
#' @param y list of zero-mean per-curve responses, typically model residuals.
#' @param kfolds number of folds. Ignored when \code{fold} is given.
#' @param log_sp_range range of \code{log(sp)} searched by
#' \code{\link[stats]{optimize}}.
#' @param fold optional integer vector of length \code{length(X)} assigning each
#' curve to a fold, for reproducing or sharing a split.
#'
#' @param ridge passed on to \code{\link{cov_symm_setup}}, a relative ridge on
#' the pair design product for ill-conditioned bases; \code{0} by default.
#' @return a list with the selected \code{sp} and its \code{logsp}, the
#' \code{objective} there, the \code{fold} used, and \code{criterion}, the
#' cross-validation criterion as a function of \code{log(sp)} -- useful for
#' plotting the curve behind the choice.
#'
#' @seealso \code{\link{cov_symm_setup}}
#' @export
cov_symm_cv <- function(X, S, y, kfolds = 5L, log_sp_range = c(-5, 5),
                        fold = NULL, ridge = 0) {
  stopifnot(is.list(X), is.list(y), length(X) == length(y))
  n <- length(X)
  if(is.null(fold)) {
    kfolds <- min(as.integer(kfolds), n)
    if(kfolds < 2L)
      stop(sQuote("kfolds"), " must be at least 2 to leave curves out.")
    fold <- sample(rep_len(seq_len(kfolds), n))
  }
  fold <- as.integer(fold)
  if(length(fold) != n)
    stop(sQuote("fold"), " must assign one fold per curve (", n, ").")

  # everything that does not depend on sp is built once per fold, so the search
  # below costs one Demmler-Reinsch back-solve per candidate rather than a refit
  parts <- lapply(sort(unique(fold)), function(k) {
    train <- which(fold != k)
    test <- which(fold == k)
    test <- test[sapply(X[test], nrow) > 1]
    if(!length(test)) return(NULL)
    list(fit = cov_symm_setup(X[train], S, ridge = ridge)(y[train], return.fun = TRUE),
         X = X[test], yy = lapply(y[test], tcrossprod))
  })
  parts <- parts[!sapply(parts, is.null)]
  if(!length(parts))
    stop("no test fold contains a curve with at least two observations.")

  criterion <- function(logsp) {
    mean(vapply(parts, function(p) {
      theta <- p$fit(exp(logsp))
      pred <- lapply(p$X, wtcrossprod, w = theta)
      mean(unlist(Map(function(yy, pr) sum((yy - pr)^2) -
                        sum((diag(yy) - diag(pr))^2), p$yy, pred)))
    }, numeric(1)))
  }

  opt <- stats::optimize(criterion, log_sp_range)
  list(sp = exp(opt$minimum), logsp = opt$minimum, objective = opt$objective,
       fold = fold, criterion = criterion, log_sp_range = log_sp_range)
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
  # the L2 default (see predict_square_smooth()) only depends on the evaluation
  # grid, so it is built once here rather than once per theta -- and `newdata`
  # IS that grid, so the quadrature rule is the proper one here without the
  # caller having to say so
  if(decompose && is.null(Gramian) && is.null(Gramian_chol))
    Gramian <- L2_gramian(X, grid = newdata)

  lapply(thetas, predict_square_smooth, X,
         decompose = decompose, Gramian = Gramian,
         Gramian_chol = Gramian_chol,
         Gramian_chol_inv = Gramian_chol_inv)
}


#' Quadrature weights for an evaluation grid
#'
#' The weights of the inner product \eqn{\langle f, g\rangle = \sum_i w_i f_i
#' g_i}, normalized to \code{sum(w) == 1} so that the scale is the mean-square
#' one estimated eigenfunctions are usually reported in.
#'
#' \code{rule = "trapezoid"} gives interior points the mean of their two
#' neighbouring cell widths and the endpoints half of their single one.
#' \code{"uniform"} weighs every node alike -- the plain mean, which is what
#' \code{\link{L2_gramian}} used unconditionally before 2026-09-08. The
#' difference is not cosmetic: the plain mean gives both endpoints of the grid
#' full weight and is therefore only first-order accurate, and on a 100-point
#' grid over \eqn{[0,1]} it is off by 2.5\% on the Gram matrix of a cubic
#' spline basis with five functions and by 6.1\% with twenty, against 0.07\%
#' and 1.1\% for the trapezoid rule on the very same nodes.
#'
#' A **multivariate** grid is handled when it is a regular tensor product in
#' \code{\link{expand.grid}} order: the weights are then the outer product of
#' the marginal ones, and a coordinate held fixed (a slice through a cube)
#' carries weight one, so the rule integrates over the slice. A grid that is
#' not recognisably of that form -- an irregular subset of a grid, say, whose
#' cells all have the same area anyway -- falls back to uniform weights, which
#' the returned \code{"rule"} attribute records; check it rather than assuming.
#'
#' @param grid the evaluation grid: a numeric vector for a univariate index, or
#' a \code{data.frame}/matrix with one column per index variable.
#' @param rule \code{"trapezoid"} or \code{"uniform"}.
#'
#' @return a numeric vector of \code{NROW(grid)} positive weights summing to
#' one, with the rule actually applied in attribute \code{"rule"}.
#'
#' @examples
#' quadrature_weights(seq(0, 1, length.out = 5))   # 1/8, 1/4, 1/4, 1/4, 1/8
#'
#' @export
quadrature_weights <- function(grid, rule = c("trapezoid", "uniform")) {
  rule <- match.arg(rule)
  n <- NROW(grid)
  unif <- function() structure(rep(1 / n, n), rule = "uniform")
  if(rule == "uniform" || n < 2L) return(unif())

  # the 1-D rule, for a sorted vector of nodes. A coordinate held FIXED (one
  # unique value, as in a slice) carries no measure of its own and gets weight
  # one, so that the tensor rule integrates over the slice.
  trap1 <- function(x) {
    k <- length(x)
    if(k == 1L) return(1)
    if(is.unsorted(x)) return(NULL)
    w <- c(x[2L] - x[1L],
           if(k > 2L) x[3L:k] - x[1L:(k - 2L)],
           x[k] - x[k - 1L]) / 2
    w / sum(w)
  }

  if(is.null(dim(grid)) || NCOL(grid) == 1L) {
    w <- trap1(as.numeric(if(is.null(dim(grid))) grid else grid[[1L]]))
    return(if(is.null(w)) unif() else structure(w, rule = "trapezoid"))
  }

  # a tensor product in expand.grid() order: the first column varies fastest
  g <- as.data.frame(grid)
  u <- lapply(g, function(z) sort(unique(as.numeric(z))))
  len <- vapply(u, length, integer(1L))
  if(prod(len) != n) return(unif())
  before <- cumprod(c(1L, len))[seq_along(len)]
  after <- rev(cumprod(c(1L, rev(len))))[-1L]
  ok <- all(vapply(seq_along(u), function(j)
    isTRUE(all.equal(as.numeric(g[[j]]),
                     rep(rep(u[[j]], each = before[j]), times = after[j]))),
    logical(1L)))
  if(!ok) return(unif())
  wl <- lapply(u, trap1)
  if(any(vapply(wl, is.null, logical(1L)))) return(unif())
  w <- Reduce(function(a, b) as.vector(outer(a, b)), wl)
  structure(w / sum(w), rule = "trapezoid")
}

#' The L2 Gramian of a basis on an evaluation grid
#'
#' The matrix \eqn{\int b_i b_j} of inner products of the basis functions,
#' approximated by quadrature on the rows of \code{X}: \code{crossprod(X *
#' sqrt(w))} with weights \code{w} summing to one. It is the inner product
#' eigenfunctions of a covariance operator are conventionally defined and
#' normalized in.
#'
#' **The quadrature rule matters.** Given \code{grid} or \code{weights} the
#' trapezoid rule is used (see \code{\link{quadrature_weights}}); with neither,
#' nothing is known about the nodes and the plain mean \code{crossprod(X)/nrow(X)}
#' is all that is left -- which weighs the two endpoints of a grid like interior
#' points and is only first-order accurate. On a 100-point grid that is a 2.5 to
#' 6.1\% error on a cubic spline basis' Gram matrix, enough to move the
#' eigenfunctions of a fitted kernel by several degrees where the estimate is
#' otherwise accurate to one. **Pass the grid.** Until 2026-09-08 this function
#' took only \code{X} and always used the plain mean;
#' \code{\link{predict_square_smooths}} now passes its \code{newdata} on for
#' you.
#'
#' @param X a design matrix, typically \code{Predict.matrix(smooth_obj, newdata)}
#' evaluated on a grid.
#' @param grid the grid \code{X} was evaluated on, passed to
#' \code{\link{quadrature_weights}}; \code{NULL} leaves the rule unknown.
#' @param weights quadrature weights for the rows of \code{X}, overriding
#' \code{grid}; they are rescaled to sum to one.
#'
#' @return a symmetric \code{ncol(X)} by \code{ncol(X)} matrix, with the
#' quadrature rule used in attribute \code{"rule"}.
#'
#' @seealso \code{\link{quadrature_weights}}
#'
#' @export
L2_gramian <- function(X, grid = NULL, weights = NULL) {
  w <- if(!is.null(weights)) structure(weights / sum(weights), rule = "weights")
       else if(!is.null(grid)) quadrature_weights(grid)
       else structure(rep(1 / nrow(X), nrow(X)), rule = "uniform")
  structure(crossprod(X * sqrt(as.numeric(w))), rule = attr(w, "rule"))
}


#' Predict a single-tensor-product smooth
#'
#' @param theta an object representing the estimated coefficients
#' @param smooth_obj a gam smoother object with a \code{Predict.matrix} method
#' @param newdata data for prediction
#' @param decompose logical, should prediction be returned in form of its SVD?
#' @param Gramian the Gram matrix of the basis, defining the inner product the
#' decomposition is taken in. Defaults to the \eqn{L^2} Gramian
#' \code{\link{L2_gramian}(X, grid)} of the evaluation grid, so that the
#' returned functions are the \eqn{L^2} eigenfunctions of the fitted kernel --
#' the convention estimated eigenfunctions are normally reported in. Pass
#' \code{diag(ncol(X))} for the oldest default, a decomposition in the inner
#' product in which the basis functions themselves are orthonormal.
#' @param grid the evaluation grid \code{X} belongs to, used for the quadrature
#' rule of the default \code{Gramian} -- **pass it**: without it the rule falls
#' back to the plain mean over the rows of \code{X}, which weighs the endpoints
#' of a grid like interior points and is only first-order accurate. See
#' \code{\link{L2_gramian}}. \code{\link{predict_square_smooths}} passes its
#' \code{newdata} automatically.
#'
#' @export
predict_square_smooth <- function(theta, X, decompose = FALSE, Gramian = NULL,
                                  Gramian_chol = NULL,
                                  Gramian_chol_inv = NULL, grid = NULL, ...) {
  UseMethod("predict_square_smooth")
}

#' @export
predict_square_smooth.matrix <- function(theta, X, decompose = FALSE, symmetric = FALSE,
                                         Gramian = NULL,
                                         Gramian_chol = NULL,
                                         Gramian_chol_inv = NULL,
                                         grid = NULL) {

  if(!decompose) return(
    wtcrossprod(X, w = theta)
  )
  # alternatively make SVD first
  if(is.null(Gramian) && is.null(Gramian_chol))
    Gramian <- L2_gramian(X, grid = grid)
  U <- if(is.null(Gramian_chol)) chol(Gramian) else Gramian_chol
  U_ <- if(is.null(Gramian_chol_inv)) solve(U) else Gramian_chol_inv
  theta <- wtcrossprod(U, w = theta)
  X <- X %*% U_

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
                                        Gramian_chol_inv = NULL,
                                        grid = NULL) {

  if(!decompose) {
    X <- X %*% theta$vectors
    return( wtcrossprod(X, w = theta$values) )
    }
  # alternatively reduce to matrix case
  theta <- wtcrossprod(theta$vectors, w = theta$values)
  predict_square_smooth.matrix(theta, X, decompose = decompose,
                               Gramian = Gramian,
                               Gramian_chol = Gramian_chol,
                               Gramian_chol_inv = Gramian_chol_inv,
                               grid = grid)
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
