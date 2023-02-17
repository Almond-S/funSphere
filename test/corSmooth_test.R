
# library(funSphere)
library(mgcv)


# example not working -----------------------------------------------------

dat <- nlme::Earthquake

m0 <- gam(accel ~ s(distance), data = dat, fit = FALSE)
plot(gam(G = m0))

m <- gamm(accel ~ s(distance), data = dat,
          correlation = corSmooth(10, form = ~ distance | Quake,
                                  s_xt = list(bsmargin = "tp"),
                                  fix = T, G = m0, s_m = 0))

# another example ---------------------------------------------------------

dat2 <- nlme::Spruce

m0 <- gam(logSize ~ s(days), data = dat2, fit = FALSE)
plot(gam(G = m0))

m <- gamm(logSize ~ s(days), data = dat2,
          correlation = corSmooth(form = ~ days | Tree,
                                  s_xt = list(bsmargin = "tp"),
                                  G = m0, s_m = 0))
plot(m$gam)
k <- attr(m$lme$modelStruct$corStruct, "covariance_model")
