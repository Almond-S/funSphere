
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
  } else return(tcrossprod(w, sweep(y, 2, w, `*`)))
}

#' Fit functional regression model from within \code{fam}
#'
#' @param gam_prefit
#' @param cov_smooth
#' @param grouping
#' @param ...
#'
#' @return
#' @import mgcv
#' @export
#'
fam_fit <- function(gam_prefit, cov_smooth, id, cov_sp = 0,
                    quadrature_dat = NULL, quadrature_weights = NULL, ...) {

  # checks:
  stopifnot(all(id == is.integer(id)) | is.ordered(id))

  # initial fits ------------------------------------------
  init <- list()
  # initial 1/4: fit mean model
  init$mean_model <- gam(G = gam_prefit, ...)

  # prepare covariance fits
  res <- split(init$mean$residuals, id)
  idx <- split(seq_along(id), id)
  X <- lapply(idx, function(id) cov_smooth$X[id, , drop = FALSE])
  if(!is.null(quadrature_dat)) {
    X0 <- Predict.matrix(cov_smooth, quadrature_dat)
    G <- wcrossprod(X0, w = quadrature_weights)
    U <- chol(G)
    U_ <- solve(U)
  }

  # initial 2/4: fit lag0 covariance
  cov_symm_fit <- cov_symm_setup(X, cov_smooth$S)
  # get estimated coefficients
  init$cov_lag0 <- cov_symm_fit(res, cov_sp)

  # initial 3/4: fit lag1 covariance
  idgrid <- if(is.factor(id)) levels(id) else min(id):max(id)
  idgrid <- data.frame(a = head(idgrid, -1), b = tail(idgrid, -1))
  # filter id grid keeping only pairs contained int the data
  idcombs <- idgrid[idgrid$a %in% names(res) & idgrid$b %in% names(res), ]
  cov_cross_fit <- cov_cross_setup(X[idcombs$a], X[idcombs$b], cov_smooth$S)
  init$cov_lag1 <- cov_cross_fit(res[idcombs$a], res[idcombs$b], sp = cov_sp)

  # initial 4/4: estimate autocorrelation operator
  S <- cov_smooth$S[[1]]
  init$autocor <- tcrossprod(
    solve(crossprod(init$cov_lag0) + cov_sp * S, init$cov_lag0), init$cov_lag1)

  # return predict function
  init$predict <- function(what = c("mean", "cov_lag0", "cov_lag1", "autocor"),
                           newdata = NULL, decompose = FALSE, ...) {
    ret <- list()
    if("mean" %in% what)
      ret$mean <- predict(init$mean, newdata = newdata, ...)
    whatelse <- setdiff(what, "mean")
    if(length(whatelse) > 0)
      ret[whatelse] <- predict_square_smooth(init[whatelse], cov_smooth,
                                             newdata = newdata,
                                             decompose = decompose)
    ret[what]
  }

  # run EM algorithm --------------------------------------


  # post process and return -------------------------------
  init
} # fam_fit
