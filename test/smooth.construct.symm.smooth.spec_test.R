
library(mgcv)
library(funSphere)

# check out example dataset ------------------------------------------------

data(world.cities, package = "maps")

library(ggplot2)
library(patchwork)

caps <- world.cities[as.logical(world.cities$capital), ]

world <- map_data("world")

(p <- ggplot(data = NULL) +
    geom_map(
      data = world, map = world,
      aes(map_id = region),
      color = "black", fill = "lightgray", size = 0.1
    ) +
    coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 90), expand = F) +
    coord_fixed() +
    geom_point(data = caps, aes(long, lat, col = sqrt(pop)))
)


# first fit strange symmetric 2D GP -------------------------------------
m0 <- gam(sqrt(pop) ~ s(long, lat, bs = "symm",
                        xt = list(bsmargin = "gp2"), m = -3), data = caps)
predgrid <- expand.grid(long = seq(-180, 180, 10), lat = seq(-180, 180, 10))
predgrid$pop0 <- predict(m0, predgrid)^2

p + geom_raster(data = predgrid, aes(long, lat, fill = sqrt(pop0)), alpha = .96, show.legend = F)

pop0 <- matrix(predgrid$pop0, nrow = sqrt(nrow(predgrid)))
all.equal(pop0, t(pop0))


# now try covariance data -------------------------------------------------
cits <- world.cities[, -2]
set.seed(8934)
n_i <- 30
n_sample <- 50
cits$sample <- NA
ids <- seq_len(nrow(cits))
for(i in n_sample)
  cits$sample[sample(ids, size = n_i, replace = FALSE)] <- i
cits <- cits[!is.na(cits$sample), ]
cits$pop <- as.double(cits$pop)

covcaps <- make_cov_data(cits, response_var = "pop")
nrow(covcaps)

m1 <- gam(sqrt(poppop) ~ s(long1, lat1, long2, lat2, bs = "symm",
                        xt = list(bsmargin = "gp2"), m = -3), data = covcaps)

predgrid1 <- expand.grid(
  long1 = seq(-180, 180, 10), lat1 = seq(-90, 90, len = 5),
  long2 = seq(-180, 180, 10), lat2 = 0 )
predgrid1$poppop <- predict(m1, predgrid1)^2

p + geom_raster(data = predgrid1,
                aes(long1, long2,
                    fill = sqrt(poppop)),
                alpha = .96, show.legend = F) +
  facet_wrap(~lat1)

predgrid_diag <- expand.grid(long1 = seq(-180, 180, 10),
                             lat1 = seq(-90, 90, 10))
predgrid_diag[c("long2", "lat2")] <- predgrid_diag[, 1:2]

predgrid_diag$poppop <- predict(m1, predgrid_diag)^2

p + geom_raster(data = predgrid_diag,
                aes(long1, lat1,
                    fill = sqrt(poppop)),
                alpha = .9, show.legend = F)


