

# test linear array model for computation of inner product matrices

X <- splines::bs(seq(0,1, len = 40), df = 5)
matplot(X, t = "l")

XX <- get_XxXtXxX(X)

# compare with naive Kronecker product basis matrix
XxX <- kronecker(X, X)
# identify diagonal
d <- c(matrix(diag(nrow = 40), ncol = 1))
image(matrix((XxX + d)[, 20], ncol = 40), asp = 1)
# remove diagonal entries
XxX_ <- XxX[!d, ]
XxX_tXxX_ <- crossprod(XxX_)

all.equal(XX, XxX_tXxX_)
# => nice!

# simulate y
set.seed(3490)
theta <- rnorm(ncol(X))
y <- X %*% theta + rnorm(ncol(X))

XY <- get_XxXtYxY_noDiag(X, y)

YtY <- kronecker(y, y)
# remove diagonal
YtY_ <- YtY[!d, , drop = FALSE]

all.equal(YtY_, crossprod(XxX_, YtY_))



