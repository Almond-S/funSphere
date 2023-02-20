
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

dat <- nlme::Spruce

m0 <- gam(logSize ~ s(days), data = dat, fit = FALSE)
plot(gam(G = m0))

m <- gamm(logSize ~ s(days), data = dat,
          correlation = corSmooth(form = ~ days | Tree,
                                  s_xt = list(bsmargin = "tp"),
                                  G = m0, s_m = 0))
plot(m$gam)
k <- attr(m$lme$modelStruct$corStruct, "covariance_model")


# and another one ---------------------------------------------------------

dat <- nlme::Soybean

m0 <- gam(weight ~ s(Time), data = dat, fit = FALSE)
plot(gam(G = m0))

m <- gamm(weight ~ s(Time), data = dat,
          correlation = corSmooth(working_correlation = corAR1(form = ~ Time | Plot),
                                  s_xt = list(bsmargin = "tp"),
                                  s_m = 0, verbose = T))
coef(m$lme$modelStruct) <- c(3,4)

plot(m$gam)
cs <- m$lme$modelStruct$corStruct
attr(cs, "fixed") <- TRUE
m0$sp <- m$gam$sp
attr(cs, "G") <- m0
attr(cs, "verbose") <- TRUE
k <- corMatrix(cs, return.model = T)
plot(k)

setdiff(names(m$gam), names(m0))
