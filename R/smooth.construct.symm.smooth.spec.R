###################################################################################
# author: Almond Stoecker heavily building on sparseFLMM::smooth.constructor.symm.smooth.spec
# NOTE: this constructor builds on the wrapper function smoothCon provided
# by Simon Wood in package mgcv.
###################################################################################
# description: smooth construct class for smoothing with our symmetry constraint
###################################################################################

##############
# what is what
##############
# k: number of basis functions
# bsmargin: type of penalty for both directions
# m: splines and difference order
# kroneckersum: which penalty matrix should be used
## TRUE to specify a Kronecker sum penalty of the form: S1 \otimes I + I \otimes S2
## FALSE to specify a Kronecker product penalty of the form: S1 \otimes S2

###################
# constraint matrix
###################

#' Construct symmetry constraint matrix for bivariate symmetric smoothing.
#'
#' This function can be used to construct a symmetry constraint matrix that imposes
#' a (skew-)symmetry constraint on (cyclic) spline coefficients in symmetric bivariate smoothing problems and is especially
#' designed for constructing objects of the class "symm.smooth", see \code{\link[sparseFLMM]{smooth.construct.symm.smooth.spec}}.
#'
#' @details Imposing a symmetry constraint to the spline coefficients in order to obtain a reduced coefficient vector is
#' equivalent to right multiplication of the bivariate design matrix
#' with the symmetry constraint matrix obtained with function \code{make_summation_matrix}.
#' The penalty matrix of the bivariate smooth needs to be adjusted to the reduced coefficient vector
#' by left and right multiplication with the symmetry constraint matrix.
#' This function is used in the constructor function \code{\link[sparseFLMM]{smooth.construct.symm.smooth.spec}}.
#'
#'
#' @param k number of marginal basis functions.
#' @param skew logical, should the basis be constraint to skew-symmetry instead
#' of symmetry.
#' @seealso \code{\link[mgcv]{smooth.construct}} and \code{\link[mgcv]{smoothCon}} for details on constructors
#' @export
#' @author Jona Cederbaum, Almond Stoecker
#' @return A basis transformation matrix of dimension \eqn{k^2 \times G} with
#' \eqn{G<k^2} depending on the specified constraint.
#' @references Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
make_summation_matrix <- function(k, skew = FALSE){
  ind_mat <- matrix(1:k^2, nrow = k, ncol = k) # index square
  pairs <- cbind(c(ind_mat), c(t(ind_mat))) # all pairs using transposed = mirror
  cons <- pairs[pairs[, 1]<pairs[, 2], , drop = FALSE] # pairs to use
  C <- diag(k^2) # initialize matrix

  # generate summation matrix for skew-symmetric / symmetric case
  if(skew) {
    C[, cons[, 1]] <- C[, cons[, 1]] - C[, cons[, 2]]
  } else {
    C[, cons[, 1]] <- C[, cons[, 1]] + C[, cons[, 2]]
  }
  if(skew) ind_vec <- cons[,1] else
    ind_vec <- pairs[pairs[, 1] <= pairs[, 2], 1, drop = FALSE]
  C <- C[, ind_vec]

  C
}

######################
# constructor function
######################
#' Symmetric bivariate smooths constructor
#'
#' The \code{symm} class is a smooth class that is appropriate for symmetric smooths, e.g. of covariance functions,
#' using tensor-product smooths in a \code{gam} formula. A constraint matrix is constructed
#' (see \code{\link[sparseFLMM]{make_summation_matrix}}) to impose
#' a (skew-)symmetry constraint on the smooth's coefficients,
#' which considerably reduces the number of coefficients that have to be estimated.
#'
#' @details By default a symmetric bivariate B-spline smooth \eqn{g} is specified,
#' in the sense that \eqn{g(s, t) = g(t, s)}.
#' In contrast to the original implementation of the function in the package
#' \code{sparseFLMM}, this implementation also works for more general smooths
#' and any even number of arguments, i.e. \eqn{g(s1, s2, ..., t1, t2, ...) = g(t1, t2, ..., s1, s2, ...)}.
#' By setting
#' \code{s(..., bs = "symm", xt = list(skew = TRUE))}, a skew-symmetric (or anti-smmetric)
#' smooth with \eqn{g(s, t) = -g(t, s)} can be specified instead.
#'
#' The underlying procedure is the following: First, the marginal spline design matrices and the corresponding
#' marginal difference penalties are built. Second, the tensor product of the marginal design matrices is built
#' and the bivariate penalty matrix is set up. Third, the constraint matrix is applied
#' to the tensor product design matrix and to the penalty matrix.
#'
#' @param object is a smooth specification object or a smooth object.
#' @param data a data frame, model frame or list containing the values
#'  of the (named) covariates at which the smooth term is to be evaluated.
#' @param knots an optional data frame supplying any knot locations
#'  to be supplied for basis construction.
#' @seealso \code{\link[mgcv]{smooth.construct}} and \code{\link[mgcv]{smoothCon}} for details on constructors
#' @export
#' @author Jona Cederbaum, Almond Stoecker
#' @return An object of class "symm.smooth". See \code{\link[mgcv]{smooth.construct}} for the elements it will contain.
#' @references Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
#' Stoecker, Pfeuffer, Steyer, Greven (2022): Elastic Full Procrustes Analysis via Hermitian Covariance Smoothing.
#' @example test/smooth.construct.symm.smooth.spec_test.R
smooth.construct.symm.smooth.spec <- function(object, data, knots){
  ##############
  # check inputs
  ##############
  if(object$dim %% 2)
    stop("Sorry, even number of terms required.")

  if(!(object$xt$bsmargin %in% c("ps", "gp2")))
    warning("Only tested with bsmargin = 'ps' or = 'gp2', yet.
            `gp` for instance does not work." )

  #############################
  # set defaults if no optional
  # arguments are given
  #############################
  if (is.null(object$xt))
    object$xt <- list(skew = FALSE)
  if(is.null(object$xt$skew))
    object$xt$skew <- FALSE
  if(is.null(object$xt$bsmargin))
    object$xt$bsmargin <- "tp"
  if(is.null(object$xt$kroneckersum))
    object$xt$kroneckersum <- TRUE

    xids <- matrix(seq_len(object$dim), ncol = 2)

    # if(length(unique(x)) < object$bs.dim)
    #   warning("basis dimension is larger than number of unique covariates")

    #############
    # check knots
    #############
    for(i in 1:nrow(xids)) {
      term <- object$term[xids[i,]]
      isn <- sapply(term, function(i) is.null(knots[i]))
      if(sum(isn) == 1)
        knots[[term[which(!isn)]]] <- knots[[term[which(isn)]]]
      if(sum(isn) == 2 & !identical(knots[[term[1]]],
                                    knots[[term[2]]]))
          stop("number of specified knots is not equal for both margins")
    }
    if(length(knots) < object$dim) {
      message("Some knots are not provided.
              Knots from first marginal smoother are used also for the second.")
    }

    ##############################
    # build marginal design matrix
    # and marginal penalties
    ##############################
    smooths <- list()
    for(i in 1:2) {
      smooths[[i]] <- smooth.construct(eval(as.call(list(as.symbol("s"),
                                                         as.symbol(object$term[xids[,i]]),
                                                         bs = object$xt$bsmargin,
                                                         k = object$bs.dim,
                                                         xt = object$xt,
                                                         m = object$p.order))),
                                       data = data,
                                       knots = knots)
      if(i==1) {
        for(j in 1:nrow(xids)) {
          term <- object$term[xids[i,]]
          if(is.null(knots[[term[2]]])) {
            k1 <- smooths[[1]]$knots
            if(is.null(k1))
              k1 <- smooths[[1]]$knt
            if(is.null(k1))
              stop("Knots not available in first smoother. Please manually specify all knots.")
            if(is.list(k1))
              knots[[term[2]]] <- k1[[term[1]]] else
                knots[[term[2]]] <- k1
          }
        }
        object$bs.dim <- smooths[[1]]$bs.dim
      }
    }

    ############################
    # build tensor product model
    # matrix and penalty matrix
    ############################
    object$X <- tensor.prod.model.matrix(X = lapply(smooths, `[[`, "X"))

    Sm <- lapply(smooths, function(sm) sm$S[[1]])

    if(object$xt$kroneckersum){
      S <- tensor.prod.penalties(Sm)
      S <- S[[1]] + S[[2]]
    } else{
      S <- Sm[[1]]%x%Sm[[2]]
    }

    ################################################
    # constraint equal coefficients by summation
    # of columns of X and adaption of penalty matrix
    ################################################

    Z <- make_summation_matrix(k = smooths[[1]]$bs.dim,
                               skew = object$xt$skew)
    object$margin <- smooths
    bs.dim <- ncol(Z)

  #########################
  # make symm.smooth object
  #########################

  object$X <- object$X %*% Z
  object$S <- list(crossprod(Z, S) %*% Z)
  object$Z <- Z
  # object$bs.dim <- bs.dim
  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <- bs.dim - object$rank
  # no sum-to-zero constraint for skew-symm bases:
  if(object$xt$skew) object$C <- matrix(0, 0, bs.dim)

  class(object) <- "symm.smooth"
  object
}


##########################
# predict method function
##########################
# needed for functions plot.gam(), model.matrix()
# also needed when bam() is used instead of gam()
# NOTE: the object here is: gam$smooth[[i]] of class symm.smooth which
# can also be generated using smooth.construct()

#' Predict matrix method for (skew-)symmetric bivariate smooths.
#'
#' @param object is a \code{symm.smooth} object created by \code{\link{smooth.construct.symm.smooth.spec}},
#' see \code{\link[mgcv]{smooth.construct}}.
#' @param data see \code{\link[mgcv]{smooth.construct}}.
#' @seealso \code{\link[mgcv]{Predict.matrix}} and \code{\link[mgcv]{smoothCon}} for details on constructors.
#' @export
#' @author Jona Cederbaum, Almond Stoecker
Predict.matrix.symm.smooth <- function (object, data) {

  # almost identical to earlier version of Predict.matrix.symm.smooth
  # only also allowing for the skew-symmetric option
  # in make_summation_matrix

    m <- length(object$margin)
    X <- list()
    for (i in 1:m) {
      term <- object$margin[[i]]$term
      dat <- list()
      for (j in 1:length(term)) {
        dat[[term[j]]] <- data[[term[j]]]
      }
      X[[i]] <- PredictMat(object$margin[[i]], dat, n = length(dat[[1]]))
    }
    X <- tensor.prod.model.matrix(X)
    if(is.null(object$Z)) {
      Z <- make_summation_matrix(k = object$margin[[1]]$bs.dim, skew = object$xt$skew)
    } else {
      Z <- object$Z
    }

    X %*% Z

}
