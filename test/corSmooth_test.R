
# library(funSphere)
library(mgcv)


# small example -----------------------------------------------------------

dat <- nlme::Earthquake

# make sure there are no dublicates
set.seed(304)
dat$distance <- dat$distance + rnorm(nrow(dat), sd = .01*sd(dat$distance))

m0 <- gam(accel ~ s(distance), data = dat, fit = FALSE)
plot(gam(G = m0))

m <- gamm(accel ~ s(distance), data = dat,
          correlation = corSmooth(c(.1), form = ~ distance | Quake,
                                  working_correlation = corExp,
                                  s_xt = list(bsmargin = "tp"),
                                  s_m = c(0,2)))

{opar <- par(mfrow = c(1,2))
  plot(gam(G = m0), main = "gam")
  plot(m$gam, main = "gamm: corSmooth")
  par(opar)}


# another example ---------------------------------------------------------

dat <- nlme::Spruce

m0 <- gam(logSize ~ s(days), data = dat, fit = FALSE)
plot(gam(G = m0))

m <- gamm(logSize ~ s(days), data = dat, method = "REML",
          correlation = corSmooth(.1, form = ~ days | Tree,
                                  working_correlation = corAR1,
                                  s_xt = list(bsmargin = "tp")))
{opar <- par(mfrow = c(1,2))
  plot(gam(G = m0), main = "gam")
  plot(m$gam, main = "gamm: corSmooth")
par(opar)}


# and another one ---------------------------------------------------------

dat <- nlme::Soybean

m0 <- gam(weight ~ s(Time), data = dat, fit = FALSE)


m <- gamm(weight ~ s(Time), data = dat, method = "REML",
          correlation = corSmooth(value = 0.01, form =  ~ Time | Plot,
                                    working_correlation = corGaus,
                                  s_xt = list(bsmargin = "tp"),
                                  s_m = c(0,2), verbose = T))
m1 <- gamm(weight ~ s(Time), data = dat,
           correlation = corGaus(form = ~ Time | Plot))

ylim <- range(dat$weight) - coef(m$gam)["(Intercept)"]
{opar <- par(mfrow = c(2,2))
  plot(gam(G = m0), main = "gam", ylim = ylim)
  plot(m$gam, main = "gamm: CorSmooth", ylim = ylim)
  plot(m1$gam, main = "gamm: CorAR1", ylim = ylim)
par(opar)}

mat <- corMatrix(m$lme$modelStruct$corStruct)
e <- eigen(mat[[1]])
e$values

mat_ <- e$vectors %*% (e$values * t(e$vectors))
image(mat_)
image(mat[[1]])

fac <- attr(mat, "fac")[[1]]
image(tcrossprod(solve(fac)))

