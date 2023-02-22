

# library(funSphere)
library(mgcv)

dat <- nlme::Soybean

# the time series time (only three different times but hopefully enough for a starter)
dat$Year <- as.numeric()
