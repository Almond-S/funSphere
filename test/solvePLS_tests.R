

X <- splines::bs(seq(0,1, len = 10), df = 5)
matplot(X, t = "l")

# test computation of design inner product --------------------------------

XxXtXxX <- get_XxXtXxX(X)

# compare with naive Kronecker product basis matrix
XxX <- kronecker(X, X)
# identify diagonal
d <- c(matrix(diag(nrow = nrow(X)), ncol = 1))
image(matrix((XxX + d)[, 20], ncol = nrow(X)), asp = 1)
# remove diagonal entries
XxX_ <- XxX[!d, ]
XxX_tXxX_ <- crossprod(XxX_)

stopifnot(all.equal(XxXtXxX, XxX_tXxX_))
# => nice!

# test computation of inner product with response -------------------------

# simulate y
set.seed(3490)
theta <- rnorm(ncol(X))
y <- X %*% theta + rnorm(ncol(X))

XxXtYxY <- get_XxXtYxY_noDiag(X, y)

YxY <- kronecker(y, y)
# remove diagonal
YxY_ <- YxY[!d, , drop = FALSE]

XxX_tYxY_ <- crossprod(XxX_, YxY_)

stopifnot(all.equal(c(XxXtYxY), c(XxX_tYxY_)))
# nice!


# test Demmler Reinsch analogue ------------------------------------------

# construct second order difference penalty
D <- diff(diag(ncol(XxXtXxX)), differences = 2)
S <- crossprod(D)

# fit via Demmler Reinsch
drsolve <- get_demmlerreinsch_solver(XxXtXxX, S)
coefs <- drsolve(.1, XxXtYxY)

# naive fit
coefs_ <- solve(XxXtXxX + .1*S, c(XxXtYxY))

stopifnot(all.equal(c(coefs), c(coefs_)))
# Perfect!!!

# use constrained basis -------------------------------------------------
Z <- sparseFLMM::make_summation_matrix(ncol(X))
drsolve. <- get_demmlerreinsch_solver(XxXtXxX, S, Q)
coefs. <- drsolve.(.1, XxXtYxY)

coefs_trafo <- solve(crossprod(Z, XxXtXxX + .1*S) %*% Z, c(crossprod(Z, XxXtYxY)))

stopifnot(all.equal(c(coefs.), c(Z %*% coefs_trafo)))
# a little bit more shrinkage seems to be applied with transformed basis?!?
plot(coefs)
points(coefs., pch = 4, col = "red")
