
library(mgcv)

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
image(pop0, asp = 1)
all.equal(pop0, t(pop0))

