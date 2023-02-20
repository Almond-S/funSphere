
library(funSphere)
library(mgcv)

dat <- nlme::Soybean

m <- gamm(weight ~ s(Time), data = dat,
          correlation = corDynamic(corAR1(form = ~ Time | Plot)))

m0 <- gamm(weight ~ s(Time), data = dat,
          correlation = corAR1(form = ~ Time | Plot))
plot(m0$gam)
