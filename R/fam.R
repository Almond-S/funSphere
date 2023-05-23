
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

#' Fit initial functional regression model from within \code{fam}
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

  # # estimate (iid) residual variance
  # var_lag0 <- (unlist(lapply(X, apply, 1, wcrossprod, w = init$cov_lag0)))
  # init$sigma2 <- pmax(init$mean_model$sig2)

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

  # initial 5/5: estimate increment variance


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

  init$res <- res
  init$X <- X

  # post process and return -------------------------------
  init
} # fam_init



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
fam_EM <- function(init, gam_prefit, id, cov_sp = 0, maxIter = 1,
                   truncate = TRUE, verbose = FALSE,
                   ...) {
  mean_model <- init$mean_model
  autocor <- wtcrossprod(init$autocor$vectors, w = init$autocor$values)
  cov_lag0 <- wtcrossprod(init$cov_lag0$vectors, w = init$cov_lag0$values)
  sigma2 <- mean_model$sig2

  y <- split(init$mean_model$y, id)

  X <- init$X
  res <- init$res
  n <- length(X)
  XX <- lapply(X, crossprod)

  ### compute conditional means and variances
  # following Shumway and Stoffer (1982) 'An Approach to Time Series Smoothing
  # and Forecasting using the EM Algorithm'

  # forward recursion to compute means/variances conditional on previous history
  forward <- function(x_prev, P_prev, M, y) {
    # mean and covariance conditional on previous history
    x_ <- autocor %*% x_prev
    P_ <- wtcrossprod(autocor, w = P_prev) + cov_lag0
    # data update
    P_M <- tcrossprod(P_, M)
    K <- P_M %*% solve(M %*% P_M + diag(sigma2, nrow = nrow(M)))
    list(
      x = x_ + K %*% (y - M %*% x_),
      P_ = P_,
      P = P_ - tcrossprod(K, P_M)
         )
  }

  # backward recursion to compute means/variances conditional on all data
  backward <- function(full_x_prev, full_P_prev, x, P, P_, P2prev_prev, J_prev, P_prev) {
    ret <- list()
    # ret$J <- solve(P, tcrossprod(P, autocor))
    ret$J <- tcrossprod(P, autocor) %*% ginv(P) # TODO: check what to do about this
    ret$x <- x + ret$J %*% (full_x_prev - autocor %*% x)
    ret$P <- P + ret$J %*% tcrossprod(full_P_prev - P_, ret$J)
    ret$P2prev <- crossprod(P_prev, ret$J) + J_prev %*% crossprod(P2prev_prev - autocor %*% P_prev, ret$J)
    ret
  }
  # where P2prev is the covariance between the previous (t) and the current time point (t-1) given all data
  # i.e. P^n_{t,t-1} in the notation of Shumway and Stoffer (shifted by 1 wrt their Eq. A11)

  # start EM algo -----------------------------------------------------------

  is <- 1:n
  EMiter <- 1
  while(EMiter <= maxIter) {

    ### compute conditional means and variances

    # run forward recursion
    fw <- list(list(x = numeric(ncol(autocor)), P = cov_lag0))
    for(i in 1:n) {
      fw[[i+1]] <- forward(fw[[i]]$x, fw[[i]]$P, X[[i]], res[[i]])
    }

    # run backward recursion
    bw <- list()
    P_M <- tcrossprod(fw[[n+1]]$P_, X[[n]])
    K <- P_M %*% solve(X[[n]] %*% P_M + diag(sigma2, nrow = nrow(X[[n]])))
    bw[[n]] <- fw[[n+1]]
    bw[[n]]$P2prev <- (diag(ncol(autocor)) - K %*% X[[n]]) %*% autocor %*% fw[[n]]$P
    bw[[n]]$J <- with(bw[[n]],  tcrossprod(P, autocor) %*% ginv(P)) # TODO: check what to do about this
    for(i in (n-1):1) {
      bw[[i]] <- backward(bw[[i+1]]$x, bw[[i+1]]$P, fw[[i+1]]$x, fw[[i+1]]$P,
                          fw[[i+1]]$P_, bw[[i+1]]$P2prev, bw[[i+1]]$J, fw[[i+1]]$P)
    }
    # bw <- bw[-1]

    ### compute
    P <- array(sapply(bw, `[[`, "P"), dim = c(dim(cov_lag0), length(bw)))
    P2prev <- array(sapply(bw, `[[`, "P2prev"), dim = c(dim(cov_lag0), length(bw)))
    x <- array(sapply(bw, `[[`, "x"), dim = c(ncol(autocor), length(bw)))

    # response minus AR predictions
    y_tilde <- Map(function(y, M, i) y - M %*% x[, i], y, X, is)
    y_tilde <- unsplit(y_tilde, id)

    # update mean model
    gam_prefit$y <- y_tilde
    mean_model <- gam(G = gam_prefit, ...)

    # update AR components
    xx <- tcrossprod(x)
    A <- rowSums(P[,,-length(bw)], dims = 2) + xx - tcrossprod(x[,length(bw)])
    B <- rowSums(P, dims = 2) + tcrossprod(x[,-1], x[,-length(bw)])
    C <- rowSums(P[,,-1], dims = 2) + xx - tcrossprod(x[,1])

    autocor <- ginv(A, B) # TODO: should probably be regularized!?!?!
    cov_lag0 <- 1/n*(C - tcrossprod(autocor, B))
    sigma2 <- mean_model$sig2 + mean(
      mapply(function(i,MM) mean(diag(P[,,i]%*%MM)), is, XX)
    )

    EMiter <- EMiter + 1
  }

  init$autocor <- autocor
  init$cov_lag0 <- cov_lag0
  init$sigma2 <- sigma2
  init$mean_model <- mean_model
  init$x_pred <- x

  init
}


