
# library(funSphere)
library(mgcv)
library(tictoc)


# small example -----------------------------------------------------------

dat <- nlme::Earthquake

# make sure there are no dublicates
set.seed(304)
dat$distance <- dat$distance + rnorm(nrow(dat), sd = .01*sd(dat$distance))

m0 <- gam(accel ~ s(distance), data = dat)
plot(m0)

tic()
m <- gamm(accel ~ s(distance), data = dat,
          correlation = corSmooth(.1,
                                  form = ~ s(distance, bs = "tp", k = 4) | Quake,
                                  working_correlation = corGaus))
toc()
# old version of corSmooth, fitted on laptop: 14.09 sec - new one down to .7 sec!

{opar <- par(mfrow = c(1,2))
  plot(m0, main = "gam", ylim = range(dat$accel) - m0$coefficients[1])
  plot(m$gam, main = "gamm: corSmooth", ylim = range(dat$accel)- m0$coefficients[1])
  par(opar)}


# another example ---------------------------------------------------------

dat <- nlme::Spruce

m0 <- gam(logSize ~ s(days), data = dat, fit = FALSE)
plot(gam(G = m0))

m <- gamm(logSize ~ s(days), data = dat, method = "REML",
          correlation = corSmooth(.1, form = ~ s(days, bs = "tp") | Tree,
                                  working_correlation = corAR1))
{opar <- par(mfrow = c(1,2))
  plot(gam(G = m0), main = "gam")
  plot(m$gam, main = "gamm: corSmooth")
par(opar)}


# and another one ---------------------------------------------------------

dat <- nlme::Soybean

m0 <- gam(weight ~ s(Time), data = dat, fit = FALSE)


m <- gamm(weight ~ s(Time), data = dat, method = "REML",
          correlation = corSmooth(value = 0.01,
                                  form =  ~ s(Time, bs = "tp", k = 3) | Plot,
                                    working_correlation = corCAR1, verbose = T))
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
image(mat_, asp = 1)
image(mat[[1]], asp = 1)


