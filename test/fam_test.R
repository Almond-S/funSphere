
library(mgcv)

dat <- nlme::Soybean

# set up vanilla gam
gam_prefit <- gam(weight ~ s(Time), data = dat, fit = FALSE)

# set up covariance smooth
cov_smooth <- smooth.construct(s(Time), data = dat, knots = NULL)

# get grouping structure
id <- model.frame(~Plot, dat)[[1]]


init <- list()
# fit initial gam
init$mean <- gam(G = gam_prefit)

# prepare covariance fit
res <- split(init$mean$residuals, id)
idx <- split(seq_along(id), id)
X <- lapply(idx, function(id) cov_smooth$X[id, , drop = FALSE])
cov_symm_fit <- cov_symm_setup(X, cov_smooth$S)
# get estimated coefficients
theta0 <- cov_symm_fit(res, .1)

# prepare lag 1 covariance fit
idgrid <- levels(id)
idgrid <- data.frame(a = head(idgrid, -1), b = tail(idgrid, -1))
# filter id grid keeping only pairs contained int the data
idcombs <- idgrid[idgrid$a %in% names(res) & idgrid$b %in% names(res), ]
cov_cross_fit <- cov_cross_setup(X[idcombs$a], X[idcombs$b], cov_smooth$S)
theta1 <- cov_cross_fit(res[idcombs$a], res[idcombs$b], sp = .1)

# plot estimates on a grid
pdat <- list(Time = seq(min(dat$Time), max(dat$Time), len = 40))
m <- fam_fit(gam_prefit, cov_smooth, id, cov_sp = .1,
             quadrature_dat = pdat, truncate = FALSE, verbose = TRUE)
pdat <- c(pdat, m$predict(newdata = pdat, decompose = FALSE))
pdat2 <- c(pdat["Time"], m$predict(newdata = pdat, decompose = TRUE))

plot(pdat$Time, pdat$mean, t = "l")
for(i in tail(names(pdat), -2)) {
  image(pdat[[i]], asp = 1, main = i)
  contour(pdat[[i]], add = T)
}
# plot decomposed
par(mfrow = c(1,3))
for(i in tail(names(pdat2), -2)) {
  matplot(pdat2[[i]]$u, main = paste(i, "u", sep = ": "), t = "l", lwd = 10*1/seq_along(pdat2[[i]]$d))
  matplot(pdat2[[i]]$v, main = paste(i, "v", sep = ": "), t = "l", lwd = 10*1/seq_along(pdat2[[i]]$d))
  barplot(pdat2[[i]]$d, main = "d", col = seq_along(pdat2[[i]]$d))
}
