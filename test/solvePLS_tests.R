

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
drsolve. <- get_demmlerreinsch_solver(XxXtXxX, S, Z)
coefs. <- drsolve.(.1, XxXtYxY)

coefs_trafo <- solve(crossprod(Z, XxXtXxX + .1*S) %*% Z, c(crossprod(Z, XxXtYxY)))

stopifnot(all.equal(c(coefs.), c(Z %*% coefs_trafo)))
# a little bit more shrinkage seems to be applied with transformed basis?!?
plot(coefs)
points(coefs., pch = 4, col = "red")


# test inner product of two matrices --------------------------------------

X2 <- splines::bs(seq(0,1, len = 10), df = 6, degree = 1)
matplot(X2, t = "l")

X2xX1tX2xX1 <- get_X2xX1tX2xX1(X, X2)

X2xX1 <- kronecker(X2, X)
X2xX1tX2xX1_ <- crossprod(X2xX1)

stopifnot(all.equal(X2xX1tX2xX1,X2xX1tX2xX1_))

X2xX1tX2xX1. <- get_X2xX1tX2xX1_weighted(X, X2)

stopifnot(all.equal(X2xX1tX2xX1.,X2xX1tX2xX1_))


# test inner product of with Kronecker response ---------------------------

y2 <- sample(y, length(y))

X2xX1tY2xY1 <- get_X2xX1tY2xY1(X, X2, y, y2)

X2xX1tY2xY1_ <- crossprod(X2xX1, kronecker(y2,y))

stopifnot(all.equal(X2xX1tY2xY1, X2xX1tY2xY1_))

# test block inversion ----------------------------------------------------

set.seed(3490)
M <- sample(1:25, 25)

M <- list(matrix(M[1:9], ncol = 3), matrix(M[10:15], ncol = 3),
          matrix(M[16:21], nrow = 3), matrix(M[22:25], ncol = 2))
dim(M) <- c(2,2)
M_ <- blockinv(M, Dinv = lapply(M[c(1,4)], solve))

N <- do.call(rbind, apply(M, 1, do.call, what = cbind))
N_ <- solve(N)
all.equal(N_[1:3,1:3], M_[[1,1]])
all.equal(N_[1:3,4:5], M_[[1,2]])
all.equal(N_[4:5,1:3], M_[[2,1]])
all.equal(N_[4:5,4:5], M_[[2,2]])

# ... and symmetric version
M[[2,1]] <- t(M[[1,2]])
M. <- blockinv_symm(bdiag = list(M[[1,1]], M[[2,2]]), odiag = M[[1,2]])
M_ <- blockinv(M)
all.equal(M_[[1,1]], M.$bdiag[[1]])
all.equal(M_[[2,2]], M.$bdiag[[2]])
all.equal(M_[[1,2]], M.$odiag)


