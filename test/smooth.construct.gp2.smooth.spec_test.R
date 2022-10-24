
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
  geom_point(data = caps, aes(long, lat, col = sqrt(pop)))
)


# first example fitting standard GP ---------------------------------------

predgrid <- with(caps, expand.grid(long = seq(-180, 180, by = 5),
                                   lat = seq(-90, 90, by = 5)))

m0 <- gam( sqrt(pop) ~ s(long, lat, bs = "gp", m = -3), data = caps)

predgrid$pop0 <- predict(m0, newdata = predgrid)^2
p + geom_raster(data = predgrid, aes(long, lat, fill = sqrt(pop0)), alpha = .8, show.legend = F)

# now the same with Euclidean GP2 ---------------------------------------

m1 <- gam( sqrt(pop) ~ s(long, lat, bs = "gp2", m = -3), data = caps)
predgrid$pop1 <- predict(m1, newdata = predgrid)^2
p + geom_raster(data = predgrid, aes(long, lat, fill = sqrt(pop1)), alpha = .8, show.legend = F)

all.equal(predgrid$pop0, predgrid$pop1)


# and finally on the sphere -----------------------------------------------

m2 <- gam( sqrt(pop) ~ s(long, lat, bs = "gp2", m = -3,
                         xt = list(distance = geosphere::distVincentySphere)), data = caps)
predgrid$pop2 <- predict(m2, newdata = predgrid)^2
(p + geom_raster(data = predgrid, aes(long, lat, fill = sqrt(pop1)), alpha = .8, show.legend = F)) /
  (p + geom_raster(data = predgrid, aes(long, lat, fill = sqrt(pop2)), alpha = .8, show.legend = F))

