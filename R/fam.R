
#' Regression for sparse functional responses
#'
#' @param ...
#'
#' @return
#' @import mgcv
#' @export
#'
fam <- function(...) {
  # prepare arguments

  # fit model
  ret <- fam_fit(...)

  # post-process and return
  ret
    } # fam

wcrossprod <- function(x, y = NULL, w = NULL) {
  if(is.null(w))
    return(crossprod(x,y))
  if(is.null(y))
    y <- x
  if(is.matrix(w)) {
      return(crossprod(x, w %*% y))
  } else return(crossprod(w * x, y))
}

wtcrossprod <- function(x, y = NULL, w = NULL) {
  if(is.null(w))
    return(tcrossprod(x,y))
  if(is.null(y))
    y <- x
  if(is.matrix(w)) {
    return(tcrossprod(x %*% w, y))
  } else return(tcrossprod(x, sweep(y, 2, w, `*`)))
}

#' Fit functional regression model from within \code{fam}
#'
#' This returns OLS mean estimates underlying also
#' all estimates of covariance components.
#' @param gam_prefit
#' @param cov_smooth
#' @param grouping
#' @param ...
#'
#' @return
#' @import mgcv
#' @export
#'
fam_init <- function(gam_prefit, cov_smooth, id, cov_sp = 0,
                    quadrature_dat = NULL,
                    quadrature_weights = 1/length(quadrature_dat[[1]]),
                    truncate = TRUE, verbose = FALSE,
                    ...) {

  if(verbose) require(tictoc)
  ttoc <- function() {toc(); tic()}
  vcat <- function(x) {if(verbose) {cat(x); cat(" for which "); ttoc()}}
  if(verbose) tic()

  # checks:
  stopifnot(all(id == is.integer(id)) | is.ordered(id))
  if(length(cov_sp) == 1)
    cov_sp <- rep(cov_sp, 3) else
      stopifnot(length(cov_sp) == 3)

  # initial fits ------------------------------------------
  init <- list()
  # initial 1/4: fit mean model
  init$mean_model <- gam(G = gam_prefit, ...)
  vcat("Initial mean model fitted...")

  # prepare covariance fits
  res <- split(init$mean$residuals, id)
  idx <- split(seq_along(id), id)
  X <- lapply(idx, function(id) cov_smooth$X[id, , drop = FALSE])
  if(any(sapply(X, is.null)))
    stop("cov_smooth$X must not be NULL.")
  if(!is.null(quadrature_dat)) {
    X0 <- Predict.matrix(cov_smooth, quadrature_dat)
    init$Gramian <- wcrossprod(X0, w = quadrature_weights)
    init$Gramian_chol <- chol(init$Gramian)
    init$Gramian_chol_inv <- solve(init$Gramian_chol)
    vcat("Gramian computed...")
  }

  # initial 2/4: fit lag0 covariance
  cov_symm_fit <- cov_symm_setup(X, cov_smooth$S)
  vcat("lag 0 covariance estimation prepared...")
  # get estimated coefficients
  init$cov_lag0 <- cov_symm_fit(res, cov_sp[1])
  vcat("and done...")
  # if(is.null(quadrature_dat))
  eigen0 <- eigen(init$cov_lag0, symmetric = TRUE)
  # else {
  #   eigen0 <- eigen(wtcrossprod(init$Gramian_chol, w = init$cov_lag0), symmetric = TRUE)
  # }

  if(truncate) {
    which_pos <- which(eigen0$values > 0)
    eigen0$values <- eigen0$values[which_pos]
    eigen0$vectors <- eigen0$vectors[, which_pos, drop = FALSE]
  }

  # store coefficients in eigenbasis + eigenvectors of variance space
  init$cov_lag0 <- eigen0
  # if(!is.null(quadrature_dat))
  #   init$cov_lag0$vectors <- init$Gramian_chol_inv %*% eigen0$vectors

  # initial 3/4: fit lag1 covariance
  idgrid <- if(is.factor(id)) levels(id) else min(id):max(id)
  idgrid <- data.frame(a = head(idgrid, -1), b = tail(idgrid, -1))
  # filter id grid keeping only pairs contained int the data
  idcombs <- idgrid[idgrid$a %in% names(res) & idgrid$b %in% names(res), ]
  cov_cross_fit <- cov_cross_setup(X[idcombs$a], X[idcombs$b], cov_smooth$S)
  vcat("lag 1 covariance estimation prepared...")
  init$cov_lag1 <- cov_cross_fit(res[idcombs$a], res[idcombs$b], sp = cov_sp[2])
  vcat("and done...")
  # project into variance space
  # if(is.null(quadrature_dat))
    init$cov_lag1 <- list(
    values = wcrossprod(eigen0$vectors, w = init$cov_lag1),
    vectors = eigen0$vectors
    )
    # else {
    #   init$cov_lag1 <- list(
    #   # first map to orthonormal basis
    #   values = wcrossprod(eigen0$vectors,
    #                       w = wtcrossprod(init$Gramian_chol, init$cov_lag1)),
    #   vectors = init$cov_lag0$vectors
    #   )
    #   }
  class(init$cov_lag1) <- "eigen"

  # initial 4/4: estimate autocorrelation operator
  # S <- cov_smooth$S[[1]]
  # init$autocor <- tcrossprod(
  #   solve(crossprod(init$cov_lag0) + cov_sp * diag(nrow = nrow(init$cov_lag1)), #S,
  #         init$cov_lag0), init$cov_lag1)
  # ... in variance eigen space
  init$autocor <- list(
    values = t(init$cov_lag0$values / (init$cov_lag0$values^2 + cov_sp[3]) * init$cov_lag1$values),
    vectors = init$cov_lag0$vectors
    )
  class(init$autocor) <- "eigen"
  vcat("autocorrelation operator estimated...")

  # return predict function
  init$predict <- function(what = c("mean", "cov_lag0", "cov_lag1", "autocor"),
                           newdata = NULL, decompose = TRUE, ...) {
    ret <- list()
    if("mean" %in% what)
      ret$mean <- predict(init$mean, newdata = newdata, ...)
    whatelse <- setdiff(what, "mean")
    if(length(whatelse) > 0)
      ret[whatelse] <- predict_square_smooths(init[whatelse], cov_smooth,
                                             newdata = newdata,
                                             decompose = decompose,
                                             Gramian = init$Gramian,
                                             Gramian_chol = init$Gramian_chol,
                                             Gramian_chol_inv = init$Gramian_chol_inv)
    ret[what]
  }

  # post process and return -------------------------------
  init
} # fam_fit



#' Run EM algorithm for estimation of timeseries FAM
#'
#' @param gam_prefit
#' @param cov_smooth
#' @param id
#' @param cov_sp
#' @param truncate
#' @param verbose
#' @param ...
#'
#' @export
fam_EM <- function(init, gam_prefit, cov_smooth, id, cov_sp = 0,
                   truncate = TRUE, verbose = FALSE,
                   ...) {
  # compute conditional means and variances
  mu <- predict(init$mean)
  # ar <-

  # obtain updated means and covariance operators

}


