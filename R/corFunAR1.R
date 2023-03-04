

#' Functional AR-1 covariance structure based on covariance smoothing
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
#'
#' @return an object of class \code{corFunAR1}, representing an covariance smoother autocorrelation structure.
#' @import mgcv nlme Matrix MASS
#' @export
#'
corFunAR1 <- function(value = 0, form = ~1, fixed = FALSE, # first: version with same smoothing parameter for auto- and cross-covariance
                      working_correlation = corCAR1,
                      working_control = list(),
                      s_xt = list(), s_k = -1, s_m = NA, verbose = FALSE) {
  if (any(value < 0)) {
    stop("penalty parameter for covariance smoothing must be non-negative")
  }
  value <- notLog2(value)

  attr(value, "s_args") <- list(xt = s_xt, k = s_k, m = s_m)
  attr(value, "verbose") <- verbose
  class(value) <- c("corFunAR1", "corStruct")

  corDynamic(value, working_correlation = working_correlation,
             working_control = working_control,
             form = form, fixed = fixed)
}


#' @export
#' @import nlme
# #' @rdname nlme::coef.corStruct
#'
coef.corFunAR1 <- function (object, unconstrained = TRUE, ...) {
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

#' @export
#' @import nlme
#'
Initialize.corFunAR1 <- function(object, data, ...) {
  object <- NextMethod()
  covar <- getCovariate(object)[[1]]
  which_int <- names(covar)[sapply(as.list(covar), is.integer)]
  if(length(which_int) != 1)
    stop("Exactly one covariate, corresponding to the discrete time of the AR process,
         has to be provided as integer.")
  attr(object, "ARtime") <- which_int
  object
}


row_tensor_square <- function(x) x[, rep(1:ncol(x), each = ncol(x))] *
  x[, rep(1:ncol(x), ncol(x))]


#' @export
#' @import mgcv nlme Matrix
# #' @rdname nlme::corMatrix.corStruct
#'
corMatrix.corFunAR1 <- function(object, covariate = getCovariate(object),
                                corr = TRUE, # named to be consistent with other corMatrix methods
                                # -> if corr = FALSE, cholesky factor of precision is computed
                                covariance = TRUE,
                                orthodat = NULL, ...) {

  if(!covariance)
    stop("Currently only covariance matrices and no correlation
                       matrices are computed.")

  # get residuals
  Residuals <- attr(object, "residuals") # assigned by update.corDynamic_init / update.corSmooth
  if(is.null(Residuals))
    Residuals <- c(attr(object, "get_residuals")())

  # check whether models needs to be refit
  refit <- is.null(attr(object, "model"))
  if(!refit) {
    oldpars <- as.vector(attr(object, "model")$smooth[[1]]$sp)
    refit <- !(all.equal(oldpars, coef(object, unconstrained = FALSE)) == TRUE)
  }

  # name of the time variable of the AR process
  ARtime <- attr(object, "ARtime")
  covnames <- setdiff(all.vars(getCovariateFormula(object)), ARtime)
  covnames_ <- paste0(covnames, "_")

  grps <- getGroups(object)

  if(refit) {

        if(!is.list(Residuals)) {
          Residuals <- split(Residuals, grps)
        }

        ## first estimate lag 0 covariance as in corSmooth -------------------------

        # build covariance data (assuming vector covariate for now)
        make_lag0_data <- function(x, r) {
          if(nrow(x) < 2)
            return(NULL)
          idx <- combn(seq_len(nrow(x)), 2)
          d <- x[idx[1,], , drop = FALSE]
          d[covnames_] <- x[idx[2, ], , drop = FALSE]
          d$residuals2 <- combn(r, 2, prod)
          d$diagonal <- 0
          x[covnames_] <- x
          x$residuals2 <- r^2
          x$diagonal <- 1
          rbind(d, x)
        }

        covariate_comb <- Map(function(x, r) {
          art <- x[[ARtime]]
          x <- split(x[covnames], art)
          r <- split(r, art)
          covco <- Map(make_lag0_data, x, r)
          covco <- do.call(rbind, covco)
        }, covariate, Residuals)

        covariate_comb <- do.call(rbind, covariate_comb)

        # fit lag 0 covariance model
        args <- as.list(attr(object, "s_args"))
        args$xt <- as.list(args$xt)
        args$xt$absorb.cons <- FALSE
        kform <- as.formula(paste("residuals2 ~ 0 + s(",
                                  paste(c(covnames, covnames_), collapse = ","),
                                  ", bs = 'symm', xt = args$xt, k = args$k, m = args$m) + diagonal"))

        if(is.null(attr(object, "model_prefit_lag0"))) {
          k_pre <- gam(kform, data = covariate_comb, fit = FALSE)
          # store prefit object
          f <- attr(object, "lme_env")
          if(!is.null(f$lmeSt))
            attr(f$lmeSt$corStruct, "model_prefit_lag0") <- k_pre
        } else {
          k_pre <- attr(object, "model_prefit_lag0")
          # update response
          k_pre$y <- covariate_comb$residuals2
        }

        # fit model
        k <- gam(G = k_pre, sp = head(coef(object, unconstrained = FALSE), length(covnames)) )

        if(attr(object, "verbose")) {
          opar <- par(mfrow = c(1,2))
          plot(k, asp = 1, main = paste("Lag 0 Smoother penalty:", k$full.sp))
        }
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
          attr(f$lmeSt$corStruct, "model_lag0") <- k
        }

      # get appropriate marginal design matrices needed for several purposes later
        if(is.null(attr(object, "marginalDesign"))) {
          marginalDesign <- lapply(covariate, function(covs) {
            covs <- split(covs, covs[[ARtime]])
            lapply(covs, Predict.matrix, object = k$smooth[[1]]$margin[[1]])
          })

          # basis orthogonalization
          if(is.null(orthodat)) {
            X <- do.call(rbind, lapply(marginalDesign, do.call, what = rbind))
          } else {
            X <- predict(k$smooth[[1]]$margin[[1]], newdata = orthodat)
          }
          R <- qr.R(qr(X))
          attr(marginalDesign, "orthogonalizeDesign") <- solve(R)
          attr(marginalDesign, "orthogonalizeDesign_inv") <- R

          # store marginal design matrix
          f <- attr(object, "lme_env")
          if(!is.null(f$lmeSt))
            attr(f$lmeSt$corStruct, "marginalDesign") <- marginalDesign
        } else {
          marginalDesign <- attr(object, "marginalDesign")
        }

        ## then estimate lag 1 covariance analogously -------------------------------
browser()
        ### manually fit lag 1 covariance model using the basis of the lag 0 model k
        ## => use linear array model (Currie et al, 2006)
        # get marginal design matrices of positive definite subspace of k
        D <- attr(k, "eigen(coefMat)")$vectors
        # Xm <- lapply(marginalDesign, lapply, function(x) x%*%D)
        X <- lapply(marginalDesign, lapply, `%*%`, D)
        S <- crossprod(D, k$smooth[[1]]$margin[[1]]$S[[1]]) %*% D
        S <- tensor.prod.penalties(list(S, S))
        sp <- numeric(2)
        sp[] <- tail(coef(object, unconstrained = FALSE), length(covnames))
        S <- sp[1]*S[[1]] + sp[2]*S[[2]]
        # compute relevant quantities separately
        XxX_XxX <- lapply(X,
                          lapply, function(X) {
                            crossprod(row_tensor_square(X))
                          })
        XxX_Y2xY1 <- Map(function(X, y, d) {
          y <- split(y, d[[ARtime]])
          Map( function(X1, X2, y1, y2) {
            kronecker(crossprod(X2, y2), crossprod(X1, y1))}
          , X[-length(X)], X[-1], y[-length(y)], y[-1])
          }, X, Residuals, covariate)
        # combine
        XxX_XxX <- Reduce(`+`, unlist(XxX_XxX, FALSE))
        XxX_Y2xY1 <- Reduce(`+`, unlist(XxX_Y2xY1, FALSE))

        c1 <- solve(XxX_XxX + S, XxX_Y2xY1)

        stop("Continue with new manual fitting of lag 1 variance here.")


        # old version non-manually:

        make_lag1_data <- function(x1, x2, r1, r2) {
          dims <- c(nrow(x1), nrow(x2))
          idx <- list(seq_len(dims[1]), seq_len(dims[2]))
          idx <- expand.grid(V1 = idx[[1]],
                             V2 = idx[[2]])
          d <- cbind(x1[idx$V1, , drop = FALSE],
                     structure(x2[idx$V2, , drop = FALSE], names = paste0(names(x2), "_")))

          d$residuals2 <- kronecker(r2, r1)
          attr(d, "dims") <- dims
          d
        }

        grps <- split(grps, grps) # for reordering

        covariate_comb <- Map(function(x, r, gr) {
          art <- x[[ARtime]]
          x$grps <- gr
          x <- split(x, art)
          r <- split(r, art)
          covco <- Map(make_lag1_data,
                       x1 = x[-length(x)], x2 = x[-1],
                       r1 = r[-length(x)], r2 = r[-1])
          dims <- lapply(covco, attr, "dims")
          covco <- do.call(rbind, covco)
          attr(covco, "dims") <- dims
          covco
        }, covariate, Residuals, grps)
        cov_dims <- lapply(covariate_comb, attr, "dims")
        covariate_comb <- do.call(rbind, covariate_comb)

        args$k <- k_pre$smooth[[1]]$margin[[1]]$df

browser()

        # fit lag 1 covariance model
          kform <- as.formula(paste("residuals2 ~ 0 + ti(",
                                    paste(c(covnames, covnames_), collapse = ","),
                                    ", bs = args$xt$bsmargin,
                              k = args$k, m = args$m, mc = c(FALSE, FALSE))"))

          if(is.null(attr(object, "model_prefit_lag1"))) {
            k_pre1 <- gam(kform, data = covariate_comb, fit = FALSE)
            # store prefit object
            f <- attr(object, "lme_env")
            if(!is.null(f$lmeSt))
              attr(f$lmeSt$corStruct, "model_prefit_lag1") <- k_pre1
          } else {
            k_pre1 <- attr(object, "model_prefit_lag1")
            # update response
            k_pre1$y <- covariate_comb$residuals2
          }

          # fit model
          k1 <- gam(G = k_pre1,
                   sp = rep(tail(coef(object, unconstrained = FALSE), length(covnames)), 2) )

          if(attr(object, "verbose")) {
            plot(k1, asp = 1)
            par(opar)
          }

          # store model object
          f <- attr(object, "lme_env")
          if(!is.null(f$lmeSt)) {
            attr(f$lmeSt$corStruct, "model_lag1") <- k1
          }
  } else {

    # covariate_comb (without residuals) is also required when not fitting

    make_lag1_data <- function(x1, x2) {
      dims <- c(nrow(x1), nrow(x2))
      idx <- list(seq_len(dims[1]), seq_len(dims[2]))
      idx <- expand.grid(V1 = idx[[1]],
                         V2 = idx[[2]])
      d <- cbind(x1[idx$V1, , drop = FALSE],
                 structure(x2[idx$V2, , drop = FALSE], names = paste0(names(x2), "_")))

      attr(d, "dims") <- dims
      d
    }

    grps <- split(grps, grps) # for reordering

    covariate_comb <- Map(function(x, gr) {
      art <- x[[ARtime]]
      x$grps <- gr
      x <- split(x, art)
      covco <- Map(make_lag1_data,
                   x1 = x[-length(x)], x2 = x[-1])
      dims <- lapply(covco, attr, "dims")
      covco <- do.call(rbind, covco)
      attr(covco, "dims") <- dims
      covco
    }, covariate, grps)
    cov_dims <- lapply(covariate_comb, attr, "dims")
    covariate_comb <- do.call(rbind, covariate_comb)
  }

  # obtain fitted values for lag 0 --------------

  # extract design and coefficient matrices
  c0 <- if(refit) k$coefficients else attr(object, "model_lag0")$coefficients
  Z <- if(refit) k$smooth[[1]]$Z else attr(object, "model_lag0")$smooth[[1]]$Z
  c0 <- list(smooth = Z%*%c0[-1], nugget = c0[1])
  c0$smooth <- matrix(c0$smooth, ncol = sqrt(length(c0$smooth)))

  get_lag0_fit <- function(X) {
    # return prediction
    pred <- X %*% tcrossprod(c0$smooth, X)
    diag(pred) <- diag(pred) + c0$nugget
    pred
  }

  val0 <- lapply(marginalDesign, function(x) {
    # return inner list of covariance matrices
    lapply(x, get_lag0_fit)
  })


  # obtain fitted values for lag 1 ---------

  covariate_comb$fitted.values <- if(refit) k1$fitted.values else
    attr(object, "model_lag1")$fitted.values

  val1 <- split(covariate_comb,
                paste(covariate_comb$grps, covariate_comb$grps_))
  val1 <- Map(function(x, dims) {
    x <- split(x, paste(x[[ARtime]], x[[paste(ARtime, "_")]]))
    Map(function(x, dims) matrix(x$fitted.values, nrow = dims[1], ncol = dims[2]),
        x, dims)
  }, val1, cov_dims)

  # extend to overlapping blocks
  val01 <- list()
  for(i in seq_along(val1)) {
    val01[[i]] <- list()
    for(j in seq_along(val1[[i]])) {
      val01[[i]][[j]] <- rbind(
        cbind(val0[[i]][[j]], val1[[i]][[j]]),
        cbind(t(val1[[i]][[j]]), val0[[i]][[j+1]])
        )
    }
  }

  ## compute factor and determinant
  # => to do so: compute precision matrix

  # compute inverse variance matrices
  my_solve <- function(x) {
    x_ <- try(solve(x, silent = TRUE))
    if(inherits(x_, "try-error"))
      x_ <- ginv(x)
    x_
  }

  grp_ids <- structure(seq_along(val1), names = names(val1))
  precision <- Map(function(id) {
    dims0 <- sapply(val0[[id]], nrow)
    M <- bandSparse(n = sum(dims0), k = 0:max(dims0), diagonals = lapply(sum(dims0) - 0:max(dims0), rep, x = 0), symmetric = TRUE)
    this <- cumsum(c(1,dims0))
    for(i in seq_along(val01[[id]])) {
      M[this[i]:(this[i+2]-1), this[i]:(this[i+2]-1)] <- M[this[i]:(this[i+2]-1), this[i]:(this[i+2]-1)] + my_solve(val01[[id]][[i]])

      if(i > 0) {
        M[this[i]:(this[i+1]-1), this[i]:(this[i+1]-1)] <- M[this[i]:(this[i+1]-1), this[i]:(this[i+1]-1)] +  my_solve(val0[[id]][[i]])
      }
    }
    forceSymmetric(M)
  }, grp_ids)

  # so far covariance/correlation matrix not returned
  fac <- lapply(precision, function(x) {
    ret <- try(as.matrix(chol(x)), silent = TRUE)
    if(inherits(ret, "try-error")) {
        xe <- eigen(x)
        xe$values[xe$values < 0] <- 0
        ret <- xe$values * t(xe$vectors)
    }
    ret
    })
  attr(fac, "logDet") <- sum(sapply(precision, function(x) {
    dt <- det(x)
    if(dt < .Machine$double.eps) 0 else log(dt) # compute log determinant of covariance matrix
  }))

  if(!corr)
    return(fac)

  # otherwise complete covariance matrix

  c1 <- if(refit) matrix(k1$coefficients, ncol = sqrt(length(k1$coefficients))) else
    matrix(attr(object, "model_lag1")$coefficients,
           ncol = sqrt(length(attr(object, "model_lag1")$coefficients)))
  # get coefficients of orthogonal basis
  A0_ <- my_solve(attr(marginalDesign, "orthogonalizeDesign_inv") %*% tcrossprod(
    c0$smooth, attr(marginalDesign, "orthogonalizeDesign_inv") ))

  A <- list()
  A1_right <- tcrossprod( c1, attr(marginalDesign, "orthogonalizeDesign_inv") )
  # A1_left:
  A[[1]] <- attr(marginalDesign, "orthogonalizeDesign_inv") %*% c1
  Alpha <- A0_ %*% A1_right

  maxARtime <- max(sapply(val0, length))
  i <- 2
  while(i < maxARtime) {
    A[[i]] <- A[[i-1]] %*% Alpha
    i <- i+1
  }

  A[[1]] <- c1

  complete_cov <- function(v0, v1, Xs) {
    dims0 <- sapply(Xs, nrow)
    # make matrix of matrices
    dimslong <- expand.grid(nr = dims0, nc = dims0)
    M <- matrix(
      apply(dimslong, 1, function(x) matrix(nrow = x[1], ncol = x[2]), simplify = FALSE),
      nrow = length(dims0), ncol = length(dims0))
    # fill matrix
    diag(M) <- v0
    for(i in seq_along(v1)) {
      M[[i,i+1]] <- v1[[i]]; M[[i+1,i]] <- t(v1[[i]]) }
    k <- 2 #start at second off-diagonal
    while(k < ncol(M)) {
      for(i in seq_len(ncol(M)-k)) {
        M[[i, i+k]][] <- Xs[[i]] %*% tcrossprod( A[[k]], Xs[[i+k]] )
        M[[i+k, i]][] <- t(M[[i, i+k]])
      }
      k <- k+1
    }
    do.call(rbind, lapply(1:nrow(M), function(i) do.call(cbind, M[i, ])))
  }
  val <- Map(complete_cov, val0, val1, marginalDesign)

  attr(val, "factor") <- fac
  val
}

#' @import nlme
#' @export
corFactor.corFunAR1 <- function(object, ...) {
  if(!is.null(aux <- attr(object, "factor"))) {
    return(aux)
  }
  corMatrix(object, ..., corr = FALSE)
}

