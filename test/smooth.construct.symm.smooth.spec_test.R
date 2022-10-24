
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
    coord_fixed()
) +
  geom_point(data = caps, aes(long, lat, col = sqrt(pop)))


# first fit strange symmetric 2D GP -------------------------------------
m0 <- gam(sqrt(pop) ~ s(long, lat, bs = "symm",
                        xt = list(bsmargin = "gp2"), m = -3), data = caps)
predgrid <- expand.grid(long = seq(-180, 180, 10), lat = seq(-180, 180, 10))
predgrid$pop0 <- predict(m0, predgrid)^2

p + geom_raster(data = predgrid, aes(long, lat, fill = sqrt(pop0)), alpha = .96) +
  scale_fill_viridis_c()

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
  long1 = seq(-180, 180, 10), lat1 = seq(-90, 90, 5),
  long2 = c(-90, 0), lat2 = c(-20,0) )
predgrid2 <- predgrid1
names(predgrid2) <- names(predgrid1)[c(3:4,1:2)]

predgrid1$poppop <- predict(m1, predgrid1)^2
predgrid2$poppop <- predict(m1, predgrid2)^2

(p + geom_raster(data = predgrid1,
                aes(long1, lat1,
                    fill = sqrt(poppop)),
                alpha = .96) +
  scale_fill_viridis_c() +
  facet_grid(lat2~long2, labeller = label_both)) +
(p + geom_raster(data = predgrid2,
                aes(long2, lat2,
                    fill = sqrt(poppop)),
                alpha = .96) +
  scale_fill_viridis_c() +
  facet_grid(lat1~long1, labeller = label_both))

predgrid_sym <- expand.grid(
  long1 = seq(-180, 180, 10), lat1 = c(-20,0),
  long2 = seq(-180, 180, 10), lat2 = c(-20,0) )
predgrid_sym$poppop <- predict(m1, predgrid_sym)^2

ggplot(predgrid_sym, aes(long1, long2, fill = sqrt(poppop))) +
  geom_raster() + geom_contour(aes(z = sqrt(poppop)), col = "white") +
  scale_fill_viridis_c() +
  facet_grid(lat1 ~ lat2, labeller = label_both) +
  coord_fixed()

predgrid_diag <- expand.grid(long1 = seq(-180, 180, 10),
                             lat1 = seq(-90, 90, 10))
predgrid_diag[c("long2", "lat2")] <- predgrid_diag[, 1:2]

predgrid_diag$poppop <- predict(m1, predgrid_diag)^2

p + geom_raster(data = predgrid_diag,
                aes(long1, lat1,
                    fill = sqrt(poppop)),
                alpha = .9) +
  scale_fill_viridis_c()


