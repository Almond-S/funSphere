

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
corFunAR1 <- function(value = 0, form = ~ 1, fixed = FALSE, # first: version with same smoothing parameter for auto- and cross-covariance
                      working_correlation = corCAR1,
                      working_control = list(), verbose = FALSE) {
    # store original formula
  form0 <- form

  # and remove the AR time
  l1 <- length(form)
  l2 <- length(form[[l1]])
  l3 <- length(form[[l1]][[l2]])
  stopifnot(form[[l1]][[1]] == as.name("|"))
  if(l3 == 1) {
    ARtime <- form[[l1]][[l2]]
    form[[l1]] <- form[[l1]][[2]]
  } else {
    ARtime <- form[[l1]][[l2]][[l3]]
    form[[l1]][[l2]] <- form[[l1]][[l2]][[2]]
  }
  ARform <- ~ t
  environment(ARform) <- environment(form)
  ARform[[2]] <- ARtime
  attr(value, "ARformula") <- ARform
  attr(value, "ARtime") <- as.character(ARtime)

  attr(value, "lme_formula") <- form

  value <- corSmooth(value, working_correlation = working_correlation,
             working_control = working_control,
             form = form0, fixed = fixed, verbose = verbose)
  class(value) <- c("corFunAR1", class(value))
  value
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


#' @export
#' @import nlme
#'
Initialize.corFunAR1 <- function(object, data, ...) {

  # first initialize covariance smoothing
  object <- NextMethod()

  # change groups and Dim to truly independent groups
  attr(object, "inner_groups") <- attr(object, "groups")
  attr(object, "inner_Dim") <- attr(object, "Dim")

  ngroups <- length(getGroupsFormula(object, TRUE))
  attr(object, "groups") <- if(ngroups>1)
    getGroups(data, formula(object), level = ngroups - 1) else
      factor(rep(1, attr(object, "lag0Dim")$N))
  attr(object, "Dim") <- Dim(object, attr(object, "groups"))
  attr(object, "groups") <- ordered(attr(object, "groups"),
                                    levels = unique(attr(object, "groups")))

  grouptable <- split(attr(object, "groups"), attr(object, "inner_groups"))
  grouptable <- sapply(grouptable, function(x) as.character(x[1]))
  grouptable <- ordered(grouptable, levels = unique(grouptable))
  attr(object, "grouptable") <- grouptable

  e <- environment(attr(object, "solvePLS"))
  # re-organize list of design matrices into groups
  X <- split(e$X, grouptable)
  S <- e$S
  nrowX <- lapply(X, sapply, nrow)
  ncolX <- sqrt(ncol(S))

  ### TODO: GET TEMPORAL ORDER STRAIGHT HERE!!!!

  X2xX1tX2xX1 <- array(0, dim = dim(S))
  for(i in seq_len(length(X) - 1)) {
    X2xX1tX2xX1 <- X2xX1tX2xX1 + get_X2xX1tX2xX1(X[[i]], X[[i+1]])
  }

  DRsolve2 <- get_demmlerreinsch_solver(X2xX1tX2xX1, S)

  # prepare fitting function
  attr(object, "solveLag1PLS") <- function(sp, y) { # U is a basis trafo matrix
    stopifnot(is.list(y) & length(y) == length(X))
    # compute "Xy"
    X2xX1tY2xY1 <- matrix(0, nrow = nrow(X2xX1tX2xX1))
    # again omit X with only one observation
    for(i in seq_len(length(X) - 1)) {
      X2xX1tY2xY1 <- X2xX1tY2xY1 + get_X2xX1tY2xY1(X[[i]], X[[i+1]], y[[i]], y[[i+1]])
    }
    # solve PLS to get coefficient matrix
    matrix(DRsolve2(sp, X2xX1tY2xY1), ncol = ncolX, nrow = ncolX)
  }

  object
}

#' @export
#' @import mgcv nlme Matrix
# #' @rdname nlme::corMatrix.corStruct
#'
corMatrix.corFunAR1 <- function(object, covariate = getCovariate(object),
                                corr = TRUE, # named to be consistent with other corMatrix methods
                                # -> if corr = FALSE, Cholesky factor of precision is computed
                                covariance = TRUE,
                                orthodat = NULL, ...) {

  if(!covariance)
    stop("Currently only covariance matrices and no correlation
                       matrices are computed.")

  coefs <- attr(object, "coefficients")
  # check whether models needs to be refit
  refit <- is.null(coefs)
  if(!refit) {
    oldpars <- attr(coefs, "sp")
    refit <- !(all.equal(oldpars, coef(object, unconstrained = FALSE)) == TRUE)
  }

  if(refit) {
    # get residuals
    Residuals <- attr(object, "residuals") # assigned by update.corDynamic_init / update.corSmooth
    if(is.null(Residuals))
      Residuals <- c(attr(object, "get_residuals")())

    grps <- getGroups(object)
    dm <- Dim(object)

    ## first estimate lag 0 covariance via corSmooth -------------------------

    # set groups & Dim to inner lag 0 structure
    attr(object, "groups") <- attr(object, "inner_groups")
    attr(object, "Dim") <- attr(object, "inner_Dim")

    val0 <- corMatrix.corSmooth(object)
    ecoefs <- attr(val0, "coefficients")
    sigma2 <- attr(val0, "sigma2noise")

    # get design matrices transformed to the level of ecoefs$values
    XU <- lapply( environment(attr(object, "solveLag1PLS"))$X, `%*%`, ecoefs$vectors)

    ## reorganize into groups ------------------------------------------------

    gt <- attr(object, "grouptable")
    Residuals <- split(Residuals, gt)
    val0 <- split(val0, gt)
    XU <- split(XU, gt)

    ## then estimate lag 1 covariance analogously ----------------------------

    if(TRUE) { # Option 1: do complete fit and then project
      coefs <- attr(object, "solveLag1PLS")(
        sp = coef(object, unconstrained = FALSE),
        y = Residuals
      )
      # project coefficients to the level of ecoefs$values for the transformed TP basis
      coefs <- crossprod(ecoefs$vectors, coefs) %*% ecoefs$vectors
    } else { # Option 2: first transform basis and then solve PLS without preparation
      # get transformed penalty matrix
      S <- lapply(attr(object, "smooth")$S, function(s)
        crossprod(ecoefs$vectors, s) %*% ecoefs$vectors)
      S <- tensor.prod.penalties(rep(S, 2))[[1]]

      XU2xXU1tXU2xXU1 <- array(0, dim = dim(S))
      XU2xXU1tY2xY1 <- matrix(0, nrow = nrow(XU2xXU1tXU2xXU1))
      for(i in seq_len(length(X) - 1)) {
        XU2xXU1tXU2xXU1 <- XU2xXU1tXU2xXU1 + get_X2xX1tX2xX1(XU[[i]], XU[[i+1]])
        XU2xXU1tY2xY1 <- XU2xXU1tY2xY1 + get_X2xX1tY2xY1(XU[[i]], XU[[i+1]], y[[i]], y[[i+1]])
      }

      coefs <- matrix(solve(XU2xXU1tXU2xXU1 + coef(object, unconstrained = FALSE) * S,
                            XU2xXU1tY2xY1),
                      ncol = ncol(ecoefs$vectors))
    } # end coefficient computation

    browser()
    # get estimated lag 1 covariance surfaces
    val1 <- Map(function(X1, X2) X1 %*% tcrossprod(coefs, X2),
                XU[-length(XU)], XU[-1])

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


  # compute precision matrix ------------------------------------------------

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

