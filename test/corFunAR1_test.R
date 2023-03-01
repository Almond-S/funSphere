

# library(funSphere)
library(mgcv)

dat <- nlme::Soybean
set.seed(490)
dat$Time <- dat$Time + rnorm(nrow(dat), sd = .1*sd(dat$Time))

# the time series time (only three different times but hopefully enough for a starter)
dat$Year <- as.integer(dat$Year)

cs <- corFunAR1(form = weight~ Time + Year | Variety, s_xt = list(bsmargin = "tp"))
covar <- getCovariate(cs, data = dat)
resp <- getResponse(cs, data = dat)
cs <- Initialize(cs, dat)

FunMat <- corMatrix(cs, corr = F)
Matrix::image(FunMat[[1]])



