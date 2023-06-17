

library(mgcv)
library(funSphere)

dat <- nlme::Soybean

# set up vanilla gam
gam_prefit <- gam(weight ~ s(Time), data = dat, fit = FALSE)

# set up covariance smooth
cov_smooth <- smooth.construct(s(Time, k = 10, bs = "ts"), data = dat, knots = NULL)

# get grouping structure
id <- model.frame(~Plot, dat)[[1]]


# manually test computation of covariance components ----------------------

# fit initial gam
meanfit <- gam(G = gam_prefit)

# prepare covariance fit
res <- split(meanfit$residuals, id)
idx <- split(seq_along(id), id)
X <- lapply(idx, function(id) cov_smooth$X[id, , drop = FALSE])
cov_symm_fit <- cov_symm_setup(X, cov_smooth$S)
cov_symm_fit <- cov_symm_fit(res, return.fun = TRUE)


# Hyperparameter tuning on a test set ---------------------------------------

set.seed(390849)
train <- rep(FALSE, length(X))
train[sample(seq_along(X), 40)] <- TRUE
cov_symm_train <- cov_symm_setup(X[train], cov_smooth$S)
cov_symm_train <- cov_symm_train(res[train], return.fun = TRUE)

track <- new.env()
track$sp <- list()
track$err <- list()

SSE <- function(logsp) {
  theta <- cov_symm_train(exp(logsp))
  res2_test <- lapply(res[-train], tcrossprod)
  pred_test <- lapply(X[-train], wtcrossprod, w = theta)
  error <- unlist(Map(function(res2, pred) sum((res2 - pred)^2) -
                        sum((diag(res2) - diag(pred))^2), res2_test, pred_test))
  error <- mean(error)
  track$sp[[length(track$sp)+1]] <- exp(logsp)
  track$err[[length(track$err)+1]] <- error
  error
}

opt <- optimize(SSE, c(-5,5))
plot(log(unlist(track$sp)), unlist(track$err), t = "b")


# fit on all data and predict ---------------------------------------------

theta <- cov_symm_fit(exp(opt$minimum))

# get estimated coefficients
griddat <- data.frame(Time = seq(min(dat$Time), max(dat$Time), len = 100))
covsurf <- predict_square_smooths(list(theta), cov_smooth, griddat)[[1]]
image(covsurf, asp = 1)

# first compute correct Grammian
Xgrid <- Predict.matrix(cov_smooth, griddat)
coveigen <- predict_square_smooths(list(theta), cov_smooth, griddat,
                                   decompose = TRUE,
                                   Gramian = crossprod(Xgrid) / 100)[[1]]
matplot(coveigen$v[, 1:4], lwd = log(coveigen$d[1:4]), t = "l")


# estimate noise variance -------------------------------------------------

comp_size <- function(theta) {
  pd <- lapply(X, wtcrossprod, w = theta)
  mean(unlist(lapply(pd, diag)))
}
size_ <- comp_size(theta)
nugget <- mean(meanfit$residuals^2) - size_
# completely of

# instead ensure variance smaller or equal empirical ----------------------

# total variance decomposition
size <- mean(meanfit$residuals^2)
cov_gap <- function(logsp) {
  th <- cov_symm_fit(exp(logsp))
  sz_ <- comp_size(th)
  er <- (sz_ - size)^2
  track$sp[[length(track$sp)+1]] <- exp(logsp)
  track$err[[length(track$err)+1]] <- er
  er
}

track <- new.env()
track$sp <- list()
track$err <- list()

opt <- optimize(cov_gap, c(-5,5))
plot(unlist(track$sp), unlist(track$err), t = "b")

# get new nugget
theta <- cov_symm_fit(exp(opt$minimum))
nugget <- mean(meanfit$residuals^2)-comp_size(theta)

# get estimated coefficients
griddat <- data.frame(Time = seq(min(dat$Time), max(dat$Time), len = 100))
covsurf <- predict_square_smooths(list(theta), cov_smooth, griddat)[[1]]
image(covsurf, asp = 1)

# first compute correct Grammian
Xgrid <- Predict.matrix(cov_smooth, griddat)
coveigen <- predict_square_smooths(list(theta), cov_smooth, griddat,
                                   decompose = TRUE,
                                   Gramian = crossprod(Xgrid) / 100)[[1]]
matplot(coveigen$v[, 1:4], lwd = log(coveigen$d[1:4]), t = "l")

