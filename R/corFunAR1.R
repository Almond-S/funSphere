

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

  # ## handle formula: store original and create formula without time
  # form0 <- form
  # # and remove the AR time
  # l1 <- length(form)
  # l2 <- length(form[[l1]])
  # l3 <- length(form[[l1]][[l2]])
  # stopifnot(form[[l1]][[1]] == as.name("|"))
  # if(l3 == 1) {
  #   ARtime <- form[[l1]][[l2]]
  #   form[[l1]] <- form[[l1]][[2]]
  # } else {
  #   ARtime <- form[[l1]][[l2]][[l3]]
  #   form[[l1]][[l2]] <- form[[l1]][[l2]][[2]]
  # }
  # ARform <- ~ t
  # environment(ARform) <- environment(form)
  # ARform[[2]] <- ARtime
  #
  # attr(value, "ARformula") <- ARform
  # attr(value, "ARtime") <- as.character(ARtime)

  # Initialize corSmooth for lag 0 covariances
  value <- corSmooth(value, working_correlation = working_correlation,
             working_control = working_control,
             form = form, fixed = fixed, verbose = verbose)

  # attr(value, "formula") <- form

  if(is.null(attr(value, "dynamic"))) {
    # attr(value, "original_formula") <- form0
    class(value) <- c("corFunAR1", class(value))
    # attr(attr(value, "dynamic"), "lme_formula") <- form
  } else {
    class(attr(value, "dynamic")) <- c("corFunAR1", class(attr(value, "dynamic")))
    # attr(attr(value, "dynamic"), "original_formula") <- form0
    # attr(attr(value, "dynamic"), "lme_formula") <- form
  }

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
  ## first initialize covariance smoothing for lag 0 covariances
  # form <- formula(object)
  # attr(object, "formula") <- attr(object, "original_formula")
  object <- NextMethod() # Initialize.corSmooth(object, data, ...) #
  # attr(object, "formula") <- attr(object, "formula")

  ## change groups and Dim to truly independent groups
  attr(object, "inner_groups") <- attr(object, "groups")
  attr(object, "inner_Dim") <- attr(object, "Dim")
  ngroups <- length(getGroupsFormula(object, TRUE))
  grps <-  if(ngroups>1)
    getGroups(data, formula(object), level = ngroups - 1) else
      factor(rep(1, attr(object, "lag0Dim")$N))
  attr(object, "groups") <- grps <- ordered(grps,
                                            levels = unique(grps))
  attr(object, "Dim") <- Dim(object, grps)

  ## handle formula: store original and create formula without time
  attr(object, "inner_formula") <- form <- attr(object, "formula")
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

  attr(object, "ARformula") <- ARform
  attr(object, "ARtime") <- as.character(ARtime)
  attr(object, "formula") <- form

  # store information which inner_group belongs to which group
  grouptable <- split(grps, attr(object, "inner_groups"))
  grouptable <- sapply(grouptable, function(x) as.character(x[1]))
  grouptable <- ordered(grouptable, levels = unique(grouptable))
  attr(object, "grouptable") <- grouptable

  # get AR time grid splitted into groups
  tgrid <- model.frame(attr(object, "ARformula"), data)[[1]]
  if(any(as.integer(tgrid) != tgrid))
    stop(paste("The AR process time variable ", attr(object, "ARtime"), " must be integer valued.", sep = "'"))
  tgrid <- as.integer(tgrid)
  tgrid <- split(tgrid, attr(object, "groups"))
  tgrid <- lapply(names(tgrid), function(x) paste(x, min(tgrid[[x]]):max(tgrid[[x]]), sep = "/"))
  attr(object, "ARtimegrid") <- tgrid

  e <- environment(attr(object, "solvePLS"))
  # re-organize list of design matrices into groups
  # X <- split(e$X, grouptable)
  X <- e$X
  S <- e$S
  ncolX <- sqrt(ncol(S))

  attr(object, "findTime") <- . <- function(x, group, time) x[[tgrid[[group]][[time]]]]

  X2xX1tX2xX1 <- array(0, dim = dim(S))
  for(g in seq_along(tgrid)) {
    for(i in seq_len(length(tgrid[[g]]) - 1)) {
      if(!is.null(.(X, g, i)) & !is.null(.(X, g, i+1)))
        X2xX1tX2xX1 <- X2xX1tX2xX1 + get_X2xX1tX2xX1(
          .(X, g, i), .(X, g, i+1))
    }
  }

  DRsolve2 <- get_demmlerreinsch_solver(X2xX1tX2xX1, S)

  # prepare fitting function
  attr(object, "solveLag1PLS") <- function(sp, y) { # U is a basis trafo matrix
    stopifnot(is.list(y) & length(y) == length(X))
    # compute "Xy"
    X2xX1tY2xY1 <- matrix(0, nrow = nrow(X2xX1tX2xX1))
    # again omit X with only one observation
    for(g in seq_along(tgrid)) {
      for(i in seq_len(length(tgrid[[g]]) - 1)) {
        if(!is.null(.(X, g, i)) & !is.null(.(X, g, i+1)))
          X2xX1tY2xY1 <- X2xX1tY2xY1 + get_X2xX1tY2xY1(
            .(X, g, i), .(X, g, i+1),
            .(y, g, i), .(y, g, i+1))
      }
    }
    # solve PLS to get coefficient matrix
    matrix(DRsolve2(sp, X2xX1tY2xY1), ncol = ncolX)
  }

  object
}

#' Inversion of 2x2 block matrix
#'
#' @param M list of matrices of \code{dim(M) == 2} with matching columns / rows
#' @param Dinv list of inverses of \code{M[[1,1]]} and \code{M[[2,2]]} in that order.
#'
#' @export
blockinv <- function(M, .solve = solve, Dinv = lapply(M[c(1,4)], .solve)) {
  R <- list()
  R[[1]] <- .solve(M[[1,1]] - M[[1,2]] %*% Dinv[[2]] %*% M[[2,1]])
  R[[4]] <- .solve(M[[2,2]] - M[[2,1]] %*% Dinv[[1]] %*% M[[1,2]])
  R[[2]] <- - R[[4]] %*% M[[2,1]] %*% Dinv[[1]]
  R[[3]] <- - R[[1]] %*% M[[1,2]] %*% Dinv[[2]]
  dim(R) <- c(2,2)
  R
  }

#' Inversion of symmetric 2x2 block matrix
#'
#' @param bdiag list of two square matrices \code{dim(M) == 2} forming the block diagonal
#' @param odiag matrix with \code{nrow(odiag) == nrow(bdiag[[1]])} and
#' \code{ncol(odiag) == ncol(bdiag[[2]])} forming the off diagonal block
#' @param bdiaginv list of inverses of \code{bdiag} (to allow their pre-computation).
#'
#' @export
blockinv_symm <- function(bdiag, odiag, .solve = solve, bdiaginv = lapply(bdiag, .solve)) {
  ret <- list(bdiag = list())
  ret$bdiag[[1]] <- .solve(bdiag[[1]] - odiag %*% tcrossprod( bdiaginv[[2]], odiag ))
  ret$bdiag[[2]] <- .solve(bdiag[[2]] - crossprod(odiag, bdiaginv[[1]]) %*% odiag)
  ret$odiag <- - ret$bdiag[[1]] %*% odiag %*% bdiaginv[[2]]
  ret
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
    if(!is.list(Residuals))
      Residuals <- split(Residuals, attr(object, "inner_groups"))

    grps <- getGroups(object)
    dm <- Dim(object)

    ## first estimate lag 0 covariance via corSmooth -------------------------

    # set groups & Dim to inner lag 0 structure
    attr(object, "groups") <- attr(object, "inner_groups")
    attr(object, "Dim") <- attr(object, "inner_Dim")
    attr(object, "formula") <- attr(object, "inner_formula")

    val0 <- corMatrix.corSmooth(object)
    ecoefs <- attr(val0, "coefficients")
    sigma2 <- attr(val0, "sigma2noise")

    # get design matrices transformed to the level of ecoefs$values
    XU <- lapply( environment(attr(object, "solveLag1PLS"))$X, `%*%`, ecoefs$vectors)

    ## reorganize into groups ------------------------------------------------

    # gt <- attr(object, "grouptable")
    # Residuals <- split(Residuals, gt)
    # val0 <- split(val0, gt)
    # XU <- split(XU, gt)

    # find inner_group in group corresponding to time
    tgrid <- attr(object, "ARtimegrid")
    . <- attr(object, "findTime")
    gt <- attr(object, "grouptable")

    dim0 <- lapply(lapply(split(sapply(val0, nrow), gt), append, 0, 0), cumsum)
    idx <- function(g, i) {
      this <- which(names(dim0[[g]]) == tgrid[[g]][i])
      (dim0[[g]][this-1]+1):dim0[[g]][this]
    }
    namesval0 <- split(names(val0), gt)

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
      for(g in seq_along(tgrid))
      for(i in seq_len(length(tgrid[[g]]) - 1)) {
        if(!is.null(.(X, g, i)) & !is.null(.(X, g, i+1))) {
          XU2xXU1tXU2xXU1 <- XU2xXU1tXU2xXU1 + get_X2xX1tX2xX1(.(XU, g, i), .(XU, g, i+1))
          XU2xXU1tY2xY1 <- XU2xXU1tY2xY1 + get_X2xX1tY2xY1(.(XU, g, i), .(XU, g, i+1),
                                                           .(y, g, i), .(y, g, i+1))
        }
      }

      coefs <- matrix(solve(XU2xXU1tXU2xXU1 + coef(object, unconstrained = FALSE) * S,
                            XU2xXU1tY2xY1),
                      ncol = ncol(ecoefs$vectors))
    } # end coefficient computation

    # store eigen decomposition of coefficients
    f <- attr(object, "lme_env")
    if(!is.null(f$lmeSt))
      attr(f$lmeSt$corStruct, "coefficients_lag1") <- coefs

    # get estimated lag 1 covariance surfaces
    val1 <- lapply(seq_along(tgrid), function(g) {
      lapply(seq_len(length(tgrid[[g]])-1), function(i) {
        if(!is.null(.(XU, g, i)) & !is.null(.(XU, g, i+1)))
          .(XU, g, i) %*% tcrossprod(coefs, .(XU, g, i+1))
      })
    })


    # # extend to overlapping blocks
    # val01 <- list()
    # for(i in seq_along(val1)) {
    #   val01[[i]] <- list()
    #   for(j in seq_along(val1[[i]])) {
    #     val01[[i]][[j]] <- rbind(
    #       cbind(.(val0, i, j), val1[[i]][[j]]),
    #       cbind(t(val1[[i]][[j]]), .(val0, i, j+1))
    #     )
    #   }
    # }

  if(!corr) {

    # compute precision matrix ------------------------------------------------

    # compute inverse variance matrices
    my_solve <- function(x) {
      x_ <- try(solve(x, silent = TRUE))
      if(inherits(x_, "try-error"))
        x_ <- ginv(x)
      x_
    }

    tlens <- lapply(tgrid, function(g) {
      l <- Dim(object)$len[g]
      names(l) <- g
      l[is.na(l)] <- 0
      l
      })

    prec0 <- structure(
      lapply(seq_along(tgrid), function(tg) {
      structure(
        lapply(seq_along(tgrid[[tg]]), function(i) my_solve(.(val0, tg, i))),
        names = tgrid[[tg]])
    }), names = names(tgrid))

    .. <- function(x, group, times) x[tgrid[[group]][times]]

    prec_diag <- prec0; prec_odiag <- val1
    for(g in 1:length(tgrid)) {
      for(i in 1:length(val1[[g]])) {
        M <- blockinv_symm(..(val0, g, i:(i+1)), val1[[g]][[i]],
                           .solve = my_solve, bdiaginv = prec0[[g]][i:(i+1)])
        prec_odiag[[g]][[i]][] <- M$odiag
        prec_diag[[g]][[i+1]][] <- M$bdiag[[2]]
        if(i == 1)
          prec_diag[[g]][[i]][] <- M$bdiag[[1]] else
            prec_diag[[g]][[i]][] <- prec_diag[[g]][[i]] + M$bdiag[[1]] - prec0[[g]][[i]]
      }
    }

    # precision still arranged like tgrid => re-arrange to data format and combine to matrix
    precision <- list()
    for(g in 1:length(tgrid)) {
      precision[[g]] <- bdiag(prec_diag[[g]][namesval0[[g]]])
      for(i in 1:(length(tgrid[[g]])-1)) {
        id <- list(idx(g, i), idx(g, i+1))
        id <- id[order(sapply(id, `[`, 1))]
        precision[[g]][id[[1]], id[[2]]] <- prec_odiag[[g]][[i]]
      }
      precision[[g]] <- forceSymmetric(precision[[g]])
    }

  # obtain 'Cholesky' factor ------------------------------------------------
    e <- lapply(precision, eigen, symmetric = TRUE)
    for(g in seq_along(e)) {
      e[[g]]$values <- pmin(e[[g]]$values, 1/sigma2)
      e[[g]]$values <- pmax(e[[g]]$values,  # TODO: get better lower bound
                            min(e[[g]]$values[e[[g]]$values>0]))
    }
    fac <- lapply(e, function(x) sqrt(x$values) * t(x$vectors))
    # log determinant of factor
    lD <- -1/2*sum(log(unlist(lapply(e, `[[`, "values"))))
    attr(fac, "logDet") <- lD

    # store eigen decomposition of coefficients
    f <- attr(object, "lme_env")
    if(!is.null(f$lmeSt))
      attr(f$lmeSt$corStruct, "fac") <- fac

    return(fac)
  } else { # now if(corr)

    # otherwise complete covariance matrix ------------------------------------

    # list of coefficient matrices in off-diagonal 1,2, ...
    A <- list()
    Alpha <- 1/ecoefs$values * coefs
    A[[1]] <- coefs

    maxARtime <- max(sapply(tgrid, length))
    i <- 2
    while(i < maxARtime) {
      A[[i]] <- A[[i-1]] %*% Alpha
      i <- i+1
    }

    complete_cov <- function(v0, v1, Xs, tg) {
      dims0 <- sapply(Xs, nrow)
      # make matrix of matrices
      dimslong <- expand.grid(nr = dims0, nc = dims0)
      M <- matrix(
        apply(dimslong, 1, function(x) matrix(nrow = x[1], ncol = x[2]), simplify = FALSE),
        nrow = length(dims0), ncol = length(dims0))
      # fill matrix
      diag(M) <- v0[names(Xs)]
      for(i in seq_along(v1)) {
        this <- c(which(names(Xs) == tg[i]), which(names(Xs) == tg[i+1]))
        M[[this[1],this[2]]] <- v1[[i]]; M[[this[2],this[1]]] <- t(v1[[i]]) }
      k <- 2 #start at second off-diagonal
      while(k < ncol(M)) {
        for(i in seq_len(ncol(M)-k)) {
          this <- c(which(names(Xs) == tg[i]), which(names(Xs) == tg[i+k]))
          M[[this[1], this[2]]][] <- Xs[[this[1]]] %*% tcrossprod( A[[k]], Xs[[this[2]]] )
          M[[this[2], this[1]]][] <- t(M[[this[1], this[2]]])
        }
        k <- k+1
      }
      do.call(rbind, lapply(1:nrow(M), function(i) do.call(cbind, M[i, ])))
    }
    val <- Map(complete_cov, split(val0, gt), val1, split(XU, gt), tgrid)

    # TODO: this fix should be removed after coming up with a more sophisticated solution
    ev <- lapply(val, eigen)
    for(g in seq_along(val)) {
      ev[[g]]$values <- pmax(ev[[g]]$values, min(ev[[g]]$values[ev[[g]]$values>0]))
      val[[g]] <- tcrossprod(sweep(ev[[g]]$vectors, 2, ev[[g]]$values, `*`), ev[[g]]$vectors)
    }

    val
    }
  } else { # i.e. if(!refit)

    if(corr) attr(object, "covariance") else
      attr(object, "fac")
  }
}

#' @import nlme
#' @export
corFactor.corFunAR1 <- function(object, ...) {
  if(!is.null(aux <- attr(object, "factor"))) {
    return(aux)
  }
  corMatrix(object, ..., corr = FALSE)
}

