gpT <- function(x,defn) {
  ## T matrix for Kamman and Wand Matern Spline...
  ## defn[1] < 0 signals no linear terms
  if (defn[1]<0) x[,1]*0+1 else cbind(x[,1]*0+1,x)
} ## gpT

gpE <- function(x,xk,defn = NA, distance = function(x, y) sqrt(rowSums((x - y)^2))) {
  ## Get the E matrix for a Kammann and Wand Matern spline.
  ## rho is the range parameter... set to K&W default if not supplied
  ind <- expand.grid(x=1:nrow(x),xk=1:nrow(xk))
  ## get d[i,j] the Euclidian distance from x[i] to xk[j]...
  E <- matrix(distance(x[ind$x, , drop = FALSE], xk[ind$xk, , drop = FALSE]), nrow(x), nrow(xk))
  rho <- -1; k <- 1
  sign.type <- 1
  if ((length(defn)==1&&is.na(defn))||length(defn)<1) { type <- 3 } else
    if (length(defn)>0) {
      type <- abs(round(defn[1]))
      sign.type <- sign(defn[1])
    }
  if (length(defn)>1) rho <- defn[2]
  if (length(defn)>2) k <- defn[3]

  if (rho <= 0) rho <- max(E) ## approximately the K & W choise
  E <- E/rho
  if (!type%in%1:5||k>2||k<=0) stop("incorrect arguments to GP smoother")
  if (type>2) eE <- exp(-E)
  E <- switch(type,
              (1 - 1.5*E + 0.5 *E^3)*(E <= 1), ## 1 spherical
              exp(-E^k), ## 2 power exponential
              (1 + E) * eE, ## 3 Matern k = 1.5
              eE + (E*eE)*(1+E/3), ## 4 Matern k = 2.5
              eE + (E*eE)*(1+.4*E+E^2/15) ## 5 Matern k = 3.5
  )
  attr(E,"defn") <- c(sign.type*type,rho,k)
  E
} ## gpE

#' Gaussian process smooth for general distances
#'
#' @param object object of class \code{gp2.smooth.spec} as returned by \code{s(..., bs = "gp2")}.
#' @param data data
#' @param knots knots
#'
#' @details This smooth is mostly a copy of the \code{gp.smooth} with two differences:
#' 1. via \code{s(..., xt = list(distance = d))} a distance function \code{d} can
#' be specified with two arguments each taking a vector with all supplied covariates.
#' For instance, the euclidean distance (default) is specified as
#' \code{d <- function(x,y) sqrt(sum((x-y)^2))}.
#' 2.supplied terms \code{x} and \code{knots} are not centered by substracting \code{colmeans(x)}.
#'
#' @importFrom mgcv smooth.construct
#' @return fitting object of class \code{gp2.smooth}.
#' @method smooth.construct gp2.smooth.spec
#' @export
#'
#' @example test/smooth.construct.gp2.smooth.spec_test.R
smooth.construct.gp2.smooth.spec <- function(object,data,knots)
  ## The constructor for a Kamman and Wand (2003) Matern Spline, and other GP smoothers.
  ## See also Handcock, Meier and Nychka (1994), and Handcock and Stein (1993).
{ ## deal with possible extra arguments of "gp" type smooth
  xtra <- list()

  ## object$p.order[1] < 0 signals stationary version
  if ((length(object$p.order)==1&&is.na(object$p.order))||length(object$p.order)<1) {
    stationary <- FALSE
  } else {
    stationary <- object$p.order[1] < 0
  }

  if (is.null(object$xt$max.knots)) xtra$max.knots <- 2000
  else xtra$max.knots <- object$xt$max.knots
  if (is.null(object$xt$seed)) xtra$seed <- 1
  else xtra$seed <- object$xt$seed

  ## now collect predictors
  x <- array(0,0)

  for (i in 1:object$dim) {
    xx <- data[[object$term[i]]]
    if (i==1) n <- length(xx) else
      if (n!=length(xx)) stop("arguments of smooth not same dimension")
    x<-c(x,xx)
  }

  if (is.null(knots)) { knt <- 0; nk <- 0}
  else {
    knt <- array(0,0)
    for (i in 1:object$dim) {
      dum <- knots[[object$term[i]]]
      if (is.null(dum)) { knt <- 0; nk <- 0; break} # no valid knots for this term
      knt <- c(knt,dum)
      nk0 <- length(dum)
      if (i > 1 && nk != nk0)
        stop("components of knots relating to a single smooth must be of same length")
      nk <- nk0
    }
  }
  if (nk>n) { ## more knots than data - silly.
    nk <- 0
    warning("more knots than data in an ms term: knots ignored.")
  }

  xu <- uniquecombs(matrix(x,n,object$dim),TRUE) ## find the unique `locations'
  if (nrow(xu) < object$bs.dim-1) stop(
    "A term has fewer unique covariate combinations than specified maximum degrees of freedom")
  ## deal with possibility of large data set
  if (nk==0) { ## need to create knots
    nu <- nrow(xu)  ## number of unique locations
    if (n > xtra$max.knots) { ## then there *may* be too many data
      if (nu > xtra$max.knots) { ## then there is really a problem
        rngs <- mgcv:::temp.seed(xtra$seed)
        #seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
        #if (inherits(seed,"try-error")) {
        #  runif(1)
        #  seed <- get(".Random.seed",envir=.GlobalEnv)
        #}
        #kind <- RNGkind(NULL)
        #RNGkind("default","default")
        #set.seed(xtra$seed) ## ensure repeatability
        nk <- xtra$max.knots ## going to create nk knots
        ind <- sample(1:nu,nk,replace=FALSE)  ## by sampling these rows from xu
        knt <- as.numeric(xu[ind,])  ## ... like this
        mgcv:::temp.seed(rngs)
        #RNGkind(kind[1],kind[2])
        #assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
      } else {
        knt <- xu; nk <- nu
      } ## end of large data set handling
    } else { knt <- xu;nk <- nu } ## just set knots to data
  }

  x <- matrix(x,n,object$dim)
  knt <- matrix(knt,nk,object$dim)

  ## centre the covariates...

  object$shift <- 0 #colMeans(x)
  x <- sweep(x,2,object$shift)
  knt <- sweep(knt,2,object$shift)

  ## Get distance function
  distance <- object$xt$distance
  if(is.null(distance))
    distance <- function(x, y) sqrt(rowSums((x - y)^2))

  ## Get the E matrix...
  E <- gpE(knt,knt,object$p.order, distance)
  object$gp.defn <- attr(E,"defn")

  def.k <- c(10,30,100)
  dd <- ncol(knt)
  if (object$bs.dim[1] < 0) { ## default basis dimension
    if(dd > 3) stop("No default basis dimension for GP domain of >3 dimensions.
                    Please, specify `k`.")
    object$bs.dim <- min(ncol(knt) + 1 + def.k[dd], nrow(E))
  }
  if (object$bs.dim < ncol(knt)+2) {
    object$bs.dim <- ncol(knt)+2
    warning("basis dimension reset to minimum possible")
  }
  object$null.space.dim <- if (stationary) 1 else ncol(knt) + 1

  k <- object$bs.dim - object$null.space.dim

  if (k < nk) {
    er <- slanczos(E,k,-1) ## truncated eigen decomposition of E
    D <- diag(c(er$values,rep(0,object$null.space.dim))) ## penalty matrix
  } else { ## no point using eigen-decomp
    D <- matrix(0,object$bs.dim,object$bs.dim)
    D[1:k,1:k] <- E  ## penalty
    er <- list(vectors=diag(k)) ## U is identity here
  }
  rm(E)

  object$S <- list(S=D)

  object$UZ <- er$vectors ## UZ - (original params) = UZ %*% (working params)

  object$knt = knt ## save the knots
  object$df <- object$bs.dim
  object$rank <- k
  object$distance <- distance
  class(object)<-"gp2.smooth"

  object$X <- Predict.matrix.gp2.smooth(object,data)

  object
} ## end of smooth.construct.gp2.smooth.spec

#' @importFrom mgcv Predict.matrix
#' @method Predict.matrix gp2.smooth
#' @export
Predict.matrix.gp2.smooth <- function(object,data)
  # prediction method function for the gp (Matern) smooth class
{ nk <- nrow(object$knt) ## number of 'knots'

## get evaluation points....
for (i in 1:object$dim) {
  xx <- data[[object$term[i]]]
  if (i==1) { n <- length(xx)
  x <- matrix(xx,n,object$dim)
  } else {
    if (n!=length(xx)) stop("arguments of smooth not same dimension")
    x[,i] <- xx
  }
}
x <- sweep(x,2,object$shift) ## apply centering

if (n > nk) { ## split into chunks to save memory
  n.chunk <- n %/% nk
  k0 <- 1
  for (i in 1:n.chunk) { ## build predict matrix in chunks
    ind <- 1:nk + (i-1)*nk
    Xc <- gpE(x=x[ind,,drop=FALSE],xk=object$knt,object$gp.defn, object$distance)
    Xc <- cbind(Xc%*%object$UZ,gpT(x=x[ind,,drop=FALSE],object$gp.defn))
    if (i == 1) X <- matrix(0,n,ncol(Xc))
    X[ind,] <- Xc
  } ## finished size nk chunks

  if (n > ind[nk]) { ## still some left over
    ind <- (ind[nk]+1):n ## last chunk
    Xc <- gpE(x=x[ind,,drop=FALSE],xk=object$knt,object$gp.defn, object$distance)
    Xc <- cbind(Xc%*%object$UZ,gpT(x=x[ind,,drop=FALSE],object$gp.defn))
    X[ind,] <- Xc
    #X <- rbind(X,Xc);rm(Xc)
  }
} else {
  X <- gpE(x=x,xk=object$knt,object$gp.defn, object$distance)
  X <- cbind(X%*%object$UZ,gpT(x=x,object$gp.defn))
}
X
} ## end of Predict.matrix.gp2.smooth
