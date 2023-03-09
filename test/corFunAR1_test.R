

# library(funSphere)
library(mgcv)

dat <- nlme::Soybean
set.seed(490)
dat$Time <- dat$Time + rnorm(nrow(dat), sd = .1*sd(dat$Time))

# the time series time (only three different times but hopefully enough for a starter)
dat$Year <- as.integer(dat$Year)

dat <- dat[order(dat$Time), ]

cs <- corFunAR1(form = weight ~ s(Time) | Variety / Year)
attr(cs, "lme_formula")

covar <- getCovariate(cs, data = dat)
resp <- getResponse(cs, data = dat)
cs <- Initialize(cs, dat)

# compare to standard
car <- corCAR1(form = ~ Time)
car <- Initialize(car, dat)

FunMat <- corMatrix(cs, corr = F)
Matrix::image(crossprod(FunMat[[1]]))
persp(crossprod(FunMat[[1]]))

# try model fit

dat$Time <- rnorm(nrow(dat), dat$Time, sd = .0001)

m0 <- gamm(weight ~ s(Time), data = dat, method = "REML",
           correlation = corSmooth(value = .01,
                                   form = ~ s(Time, k = 3) | Variety / Year,
                                   working_correlation = corCAR1, verbose = T),
           control = lmeControl(maxIter = 1, returnObject = T))

m <- gamm(weight ~ s(Time), data = dat, method = "REML",
          correlation = corFunAR1(value = .01,
                                  form = ~ Time + Year | Variety,
                                  working_correlation = corCAR1,
                                  s_xt = list(bsmargin = "tp"),
                                  s_m = c(0,2), verbose = T),
          control = lmeControl(maxIter = 1, returnObject = T))
plot(m$gam)
