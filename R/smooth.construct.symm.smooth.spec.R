###################################################################################
# author: Almond Stoecker heavily building on sparseFLMM::smooth.constructor.symm.smooth.spec
# NOTE: this constructor builds on the wrapper function smoothCon provided
# by Simon Wood in package mgcv.
###################################################################################
# description: smooth construct class for smoothing with our symmetry constraint
# NOTE: this class is so far applicable to auto-covariances only.
# It is implemented for tensor product P-splines and
# allows for two different penalty types.
# So far, it assumes the same number and type of basis functions in each direction.
###################################################################################

######################
# underlying procedure
######################
# 1.) For each auto-covariance, we first build the marginal spline design matrices and the corresponding
# marginal difference penalties.
# 2.) The tensor product of the marginal design matrices is built and the bivariate penalty matrix is set up.
# 3.) The constraint matrix is applied to the tensor product design matrix and to the penalty matrix.

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
#' The \code{symm} class is a smooth class that is appropriate for symmetric bivariate smooths, e.g. of covariance functions,
#' using tensor-product smooths in a \code{gam} formula. A constraint matrix is constructed
#' (see \code{\link[sparseFLMM]{make_summation_matrix}}) to impose
#' a (skew-)symmetry constraint on the (cyclic) spline coefficients,
#' which considerably reduces the number of coefficients that have to be estimated.
#'
#' @details By default a symmetric bivariate B-spline smooth \eqn{g} is specified,
#' in the sense that \eqn{g(s, t) = g(t, s)}. By setting
#' \code{s(..., bs = "symm", xt = list(skew = TRUE))}, a skew-symmetric (or anti-smmetric)
#' smooth with \eqn{g(s, t) = -g(t, s)} can be specified instead.
#' In both cases, the smooth can also be constraint to be cyclic
#' with the property \eqn{g(s, t) = g(s + c, t) = g(s, t + c)}
#' for some fixed constant \eqn{c} via specifying \code{xt = list(cyclic = TRUE)}.
#' Note that this does not correspond to specifying a tensor-product smooth from
#' cyclic marginal B-splines as given by the \code{cp}-smooth.
#' In the cyclic case, it is recommended to explicitly specify the range of the domain
#' of the smooth via the \code{knots} argument, as this determines the period and
#' often deviates from the observed range.
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
#' Submitted on arXiv.
#' @example tests/smooth.construct.symm.smooth.spec_example.R
smooth.construct.symm.smooth.spec <- function(object, data, knots){
  ##############
  # check inputs
  ##############
  if(object$dim %% 2)
    stop("Sorry, even number of terms required.")

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

    xids <- matrix(object$dim, ncol = 2)
    # x1 <- data[object$term[xids[,1]]]
    # x2 <- data[[object$term[xids[,2]]]]

    if(length(unique(x)) < object$bs.dim)
      warning("basis dimension is larger than number of unique covariates")

    #############
    # check knots
    #############
    k1 <- list()
    for(i in 1:nrow(xids)) {
      k1[[i]] <- if(is.null(knots[[object$term[1]]]))
        knots[[object$term[2]]] else knots[[object$term[1]]]
      k2 <- knots[[object$term[2]]]
      if(!is.null(k2)) {
        if(!identical(k1[[i]], k2))
          stop("number of specified knots is not equal for both margins")
      }
      # if(is.null(k1[[i]])) k1[[i]] <- range(data[object$term])
    }

    object$knots <- k1
    names(object$knots) <- object$term[xids[1,]]

    ##############################
    # build marginal design matrix
    # and marginal penalties
    ##############################
    smooths <- list()
    for(i in 1:2) smooth.construct(eval(as.call(list(as.symbol("s"),
                                                  as.symbol(object$term[xids[,i]]),
                                                  bs = object$xt$bsmargin,
                                                  k = object$bs.dim,
                                                  m = object$p.order))),
                                   data = data,
                                    knots = object$knots)
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

    Z <- make_summation_matrix(k = object$bs.dim,
                               skew = object$xt$skew)
    object$margin <- smooths
    object$m <- m
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

  if(length(object$term) == 2) {
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
      Z <- make_summation_matrix(k = object$bs.dim, skew = object$xt$skew)
    } else {
      Z <- object$Z
    }

    X %*% Z
  }

}
