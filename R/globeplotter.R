
#' prepare utilities for plotting Goode holomosine projections (with tmap)
#'
#' @param type the ocean (default) or land version of the projection
#'
#' @return list containing different utlities
#' @import stars sf
#' @export
goode <- function(type = c("ocean", "land"), longrat = 30, latgrat = longrat/2) {
  type <- match.arg(type)

  # relevant coordinate reference systems (CRS)
  crs_default <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  crs_goode <- switch(type,
                      land = "+proj=igh +towgs84=0,0,0",
                      ocean = "+proj=igh_o +lon_0=-160")

  # Prepare polygon of Goode projection cuts
  pm <- .01*c(1,-1)
  goode_cuts <- switch(type,
                       # code from https://wilkelab.org/practicalgg/articles/goode.html
                       land = cbind(
                         longs = c(
                           rep(180, 181), # right side down
                           rep(80+pm, each = 91), # third cut bottom
                           rep(-20+pm, each = 91), # second cut bottom
                           rep(100+pm, each = 91), # first cut bottom
                           rep(-180, 181), # left side up
                           rep(40-pm, each = 91), # cut top
                           180 # close
                         ),
                         lats = c(
                           90:-90, # right side down
                           -90:0, 0:-90, # third cut bottom
                           -90:0, 0:-90, # second cut bottom
                           -90:0, 0:-90, # first cut bottom
                           -90:90, # left side up
                           90:0, 0:90, # cut top
                           90 # close
                         )
                       ),
                       # for Goode ocean version based on Goode 1925
                       ocean = cbind(
                         longs <- c(
                           rep(20, 181), # right side down
                           rep(-70+pm, each = 91), # second cut bottom
                           rep(140+pm, each = 91), # first cut bottom
                           rep(20+pm[1], 181), # left side up
                           rep(110-pm, each = 91), # first cut top
                           rep(-100-pm, each = 91), # second cut top
                           20
                         ),
                         lats = c(
                           90:-90, # right side down
                           -90:0, 0:-90, # second cut bottom
                           -90:0, 0:-90, # first cut bottom
                           -90:90, # left side up
                           90:0, 0:90, # first cut top
                           90:0, 0:90, # second cut top
                           90 # close
                         )
                       )
  )
  # transform cuts
  bg <- st_transform(
    st_sfc(
      st_polygon(
        list(goode_cuts)),
      crs = crs_default
      ), crs = crs_goode)

  ## create graticule
  if(length(longrat) == 1)
    longrat <- seq(-180, 180, by = longrat)
  if(length(latgrat) == 1)
    latgrat <- seq(-90, 90, by = latgrat)

  grat <- st_intersection(
    st_transform(
      st_graticule(lon = longrat, lat = latgrat),
      crs = crs_goode),
    bg)

  # goodize polygon
  goodize_polygon <- function(x) {
    if(type == "ocean")
      x <- st_break_antimeridian(x, lon_0 = -160)
    x <- st_transform(x, crs = crs_goode)
    x <- st_make_valid(x)
    # restrict polygons to background polygon
    x <- st_intersection(x, bg)
  }

  # goodize raster
  goodize_raster <- function(x) st_crop(
    st_warp(
      st_transform(
        x,
        crs = crs_goode),
      crs = crs_goode),
    bg)

  # goodize coordinates
  goodize_coordinates <- function(x, longitude = "LONGITUDE", latitude = "LATITUDE") {
    coos <- as.matrix(x[, c(longitude, latitude)])
    coos <- st_sf(
      st_sfc(
        st_multipoint(coos)
        )
      )
    st_crs(coos) <- crs_default
    # warp raster to new CRS
    coos <- st_transform(coos, crs = crs_goode)

    coos <- as.data.frame(unclass(coos[[1]][[1]]))
    names(coos) <- c("LONGITUDE", "LATITUDE")
    x[c(longitude, latitude)] <- coos
    x
  }

  list(background = bg, graticule = grat,
       goodize_polygon = goodize_polygon,
       goodize_raster = goodize_raster,
       goodize_coordinates = goodize_coordinates,
       crs_goode = crs_goode, crs_default = crs_default)
}


#' Plot dataset with prediction in Goode holomosine projections
#'
#' @param x data.frame containing the `variable` to plot and
#' respective `longitude`s and `latitude`s as columns.
#' @param variable character, name of variable to illustrate.
#' @param dims character vector, names of longitude and latitude columns.
#' @param method plot method. Defaults to 'Goode_ocean', plotting the data
#' with R package `tmap` in the ocean-focussed version of the
#' the intersected holomosine projection of Goode (1925).
#'
#' @return
#' @export
#'
#' @examples
globeplotter_goode <- function(pred = NULL, pred_variables = "TEMP",
                        dat = NULL, dat_variables = NULL,
                        dims = c("LONGITUDE", "LATITUDE"),
                        return_fun = FALSE,
                        method = c("Goode_ocean", "Goode_land"),
                        pred_raster_args = list(palette = "-YlGnBu", n = 9, legend.reverse = TRUE,
                                             title = "Water temperature"),
                        dat_dots_args = list(col = "darkred"),
                        bgpolygon_borders_args = list(col = "darkgrey"),
                        land_raster_args = list("elevation", palette = "-Greys", n = 9,
                                                legend.show = FALSE),
                        world_lines_args = list(),
                        graticule_lines_args = list(col = "darkgrey")) {
  goode <- c("Goode_ocean", "Goode_land")
  method <- match.arg(method)

  # collection of simple features which define the land boundaries
  world <- st_as_sf(maps::map("world", plot = FALSE), fill = FALSE)
  data("land", package = "tmap")

  if(method %in% goode) {
    # Based on https://github.com/r-tmap/tmap-book/blob/master/code/crs_examples.R
    require(tmap)

    # relevant coordinate reference systems (CRS)
    crs_default <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    crs_goode <- switch(method,
                        Goode_land = "+proj=igh +towgs84=0,0,0",
                        Goode_ocean = "+proj=igh_o +lon_0=-160")

    # Prepare polygon of Goode projection cuts
    pm <- .01*c(1,-1)
    goode_cuts <- switch(method,
                         # code from https://wilkelab.org/practicalgg/articles/goode.html
                         Goode_land = cbind(
                           longs = c(
                             rep(180, 181), # right side down
                             rep(80+pm, each = 91), # third cut bottom
                             rep(-20+pm, each = 91), # second cut bottom
                             rep(100+pm, each = 91), # first cut bottom
                             rep(-180, 181), # left side up
                             rep(40-pm, each = 91), # cut top
                             180 # close
                           ),
                           lats = c(
                             90:-90, # right side down
                             -90:0, 0:-90, # third cut bottom
                             -90:0, 0:-90, # second cut bottom
                             -90:0, 0:-90, # first cut bottom
                             -90:90, # left side up
                             90:0, 0:90, # cut top
                             90 # close
                           )
                         ),
                         # for Goode ocean version based on Goode 1925
                         Goode_ocean = cbind(
                           longs <- c(
                             rep(20, 181), # right side down
                             rep(-70+pm, each = 91), # second cut bottom
                             rep(140+pm, each = 91), # first cut bottom
                             rep(20+pm[1], 181), # left side up
                             rep(110-pm, each = 91), # first cut top
                             rep(-100-pm, each = 91), # second cut top
                             20
                           ),
                           lats = c(
                             90:-90, # right side down
                             -90:0, 0:-90, # second cut bottom
                             -90:0, 0:-90, # first cut bottom
                             -90:90, # left side up
                             90:0, 0:90, # first cut top
                             90:0, 0:90, # second cut top
                             90 # close
                           )
                         )
                         )
    bg <- list(goode_cuts) %>%
      st_polygon() %>%
      st_sfc(
        crs = crs_default
      ) %>%
      st_transform(crs = crs_goode)

    # apply Goode projection
    if(method == "Goode_ocean")
      world_goode <- world %>% st_break_antimeridian(lon_0 = -160)
    world_goode <- st_transform(world_goode, crs = crs_goode)
    world_goode <- world_goode %>% st_make_valid()
    # restrict polygons to background polygon
    world_goode <- st_intersection(world_goode, bg)

    ## create graticule
    grat <- sf::st_graticule(lon = seq(-180, 150, by = 30), lat = seq(-90, 90, by = 30)) %>%
      st_transform(crs = crs_goode) %>%
      st_intersection(bg)

    ## create and transform spatial data of dat
    if(!is.null(dat)) {
      dat_st <- st_multipoint(as.matrix(dat[, dims])) %>% st_sfc() %>% st_sf()
      st_crs(dat_st) <- crs_default
      # warp raster to new CRS
      dat_st <- dat_st %>% st_transform(crs = crs_goode)

      if(!is.null(dat_variables))
        dat_st[, dat_variables] <- dat[, dat_variables]
    }

    ## warp and crop also land accordingly
    land_st <- land["elevation"] %>% st_transform(crs = crs_goode) %>%
      st_warp(crs = crs_goode) %>% st_crop(bg)

    ## create plot
    .pl_defaults <- list(
      pred = pred,
      pred_variables = pred_variables,
      pred_raster_args = pred_raster_args,
      dat_dots_args = dat_dots_args,
      bgpolygon_borders_args = bgpolygon_borders_args,
      land_raster_args = land_raster_args,
      world_lines_args = world_lines_args,
      graticule_lines_args = graticule_lines_args
    )

    pl <- function(pred = .pl_defaults$pred,
                   pred_variables = .pl_defaults$pred_variables,
                   pred_raster_args = .pl_defaults$pred_raster_args,
                   dat_dots_args = .pl_defaults$dat_dots_args,
                   bgpolygon_borders_args = .pl_defaults$bgpolygon_borders_args,
                   land_raster_args = .pl_defaults$land_raster_args,
                   world_lines_args = .pl_defaults$world_lines_args,
                   graticule_lines_args = .pl_defaults$graticule_lines_args) {

      ## create and transform spatial data of pred
      if(!is.null(pred)) {
        pred_st <- st_as_stars(pred, dims = c("LONGITUDE", "LATITUDE"),
                               coords = crs_default)
        st_crs(pred_st) <- crs_default
        # warp raster to new CRS
        pred_st <- pred_st %>% st_transform(crs = crs_goode) %>%
          st_warp(crs = crs_goode)
        # crop to background
        pred_st <- st_crop(pred_st, bg)
      }

      p <- tm_shape(bg) + do.call(tm_borders, bgpolygon_borders_args)

      if(!is.null(pred) & !is.null(pred_raster_args))
        p <- p + tm_shape(pred_st) +
        do.call(tm_raster, c(list(pred_variables), pred_raster_args))

      if(!is.null(dat) & !is.null(dat_dots_args))
        p <- p + tm_shape(dat_st) +
        do.call(tm_dots, c(as.list(dat_variables), dat_dots_args))

      if(!is.null(land_raster_args))
        p <- p + tm_shape(land_st) +
        do.call(tm_raster, land_raster_args)

      if(!is.null(world_lines_args))
        p <- p + tm_shape(world_goode) +
        do.call(tm_lines, world_lines_args)

      if(!is.null(graticule_lines_args))
        p <- p + tm_shape(grat) +
        do.call(tm_lines, graticule_lines_args)

      p + tm_layout(frame = FALSE, legend.outside = TRUE)
    }
  }

  if(return_fun)
    pl else
      pl(pred = pred, pred_variables = pred_variables)
} # plot_goode





#' Plot data.frame with predictions in interactive 3D chart
#'
#' @param pred
#' @param pred_variables
#' @param dat
#' @param dat_variables
#' @param dims
#' @param return_fun
#' @param method
#' @param pred_args
#' @param dat_args
#' @param land_args
#' @param world_args
#' @param graticule_args
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
globeplotter_3D <- function(pred = NULL, pred_variables = "TEMP",
                       dat = NULL, dat_variables = NULL,
                       dims = c("LONGITUDE", "LATITUDE"),
                       return_fun = FALSE,
                       pred_args = list(palette = "-YlGnBu", n = 9, legend.reverse = TRUE,
                                               title = "Water temperature"),
                       dat_args = list(col = "darkred"),
                       land_args = list("elevation", palette = "-Greys", n = 9,
                                               legend.show = FALSE),
                       world_args = list(),
                       graticule_args = list(col = "darkgrey"), ...) {

  # collection of simple features which define the land boundaries
  world <- st_as_sf(maps::map("world", plot = FALSE), fill = FALSE)
  data("land")

  # Plotly utilities --------------------------------------------------

  # hide all the axes
  empty_axis <- list(
    showgrid = FALSE,
    zeroline = FALSE,
    showticklabels = FALSE,
    title = ""
  )

  # helper function for converting polar -> cartesian
  degrees2radians <- function(degree) degree * pi / 180

  lon_shift <- 0*pi

  # set custom color scale
  colorscale <- data.frame(
    breaks = seq(0, 1, length.out = 9)
  ) %>% mutate(
    colors = rev(scales::colour_ramp(brewer.pal(9, "YlGnBu"))(breaks))
  )
  colorscale_land <- data.frame(
    breaks = seq(0, 1, length.out = 9)
  ) %>% mutate(
    colors = rev(scales::colour_ramp(brewer.pal(9, "Greys"))(breaks))
  )

  # prepare land data
  land_pl <- expand.grid(
    lon = st_get_dimension_values(land, which = "x"),
    lat = st_get_dimension_values(land, which = "y"))
  land_pl$elevation <- c(land$elevation)
  # in cartesian coos
  land_pl <- land_pl %>% mutate(
    x = cos(degrees2radians(lon) + lon_shift) * cos(degrees2radians(lat)),
    y = sin(degrees2radians(lon) + lon_shift) * cos(degrees2radians(lat)),
    z = sin(degrees2radians(lat))
  )
  # Remove z with NA elevation - otherwise oceans will be filled
  land_pl <- land_pl %>% mutate(
    z = ifelse(is.na(elevation), NA, z)
  )
  land_pl <- as.list(land_pl) %>% lapply(matrix, nrow = nrow(land$elevation))

  earthradius <- 6378137

  # combine surrounding into function
  add_maplayout <- function(p, zoom = .9, rot = pi/3) add_surface(p,
                                                                  data = land_pl,
                                                                  x = ~(1+elevation/earthradius)*x, y = ~(1+elevation/earthradius)*y, z = ~(1+elevation/earthradius)*z,
                                                                  surfacecolor = ~elevation,
                                                                  hoverinfo = "none",
                                                                  colorscale = colorscale_land,
                                                                  showscale = FALSE,
                                                                  showlegend = FALSE,
                                                                  connectgaps = FALSE
  ) %>% add_sf(
    data = world,
    x = ~ 1.001 * cos(degrees2radians(x)+lon_shift) * cos(degrees2radians(y)),
    y = ~ 1.001 * sin(degrees2radians(x)+lon_shift) * cos(degrees2radians(y)),
    z = ~ 1.001 * sin(degrees2radians(y)),
    showlegend = FALSE,
    color = I("black"), hoverinfo = "none", size = I(1)
  ) %>% layout(
    title = TeX("\\theta"),
    #title = TeX("$\\hat{\\mu}(u) = \\sum_{i=1}^n\\sum_{j=1}^{r_i} \\alpha_{i,j} \\psi(\\langle u, u_{i,j} \\rangle)$"),#"Model prediction for the mean",
    scene = list(
      xaxis = empty_axis,
      yaxis = empty_axis,
      zaxis = empty_axis,
      aspectratio = list(x = 1, y = 1, z = 1),
      camera = list(eye = list(x = 2*zoom*cos(rot), y = 1.15*zoom*sin(rot), z = 0.2*zoom))
    )
  )

  # Prepare data ------------------------------------------------------------

  if(!is.null(dat)) {
    dat_ <- list()
    dat_$x <- cos(degrees2radians(dat[[dims[1]]]) + lon_shift) * cos(degrees2radians(dat[[dims[2]]]))
    dat_$y <- sin(degrees2radians(dat[[dims[1]]]) + lon_shift) * cos(degrees2radians(dat[[dims[2]]]))
    dat_$z <- sin(degrees2radians(dat[[dims[2]]]))
    stopifnot(length(dat_variables) <= 1)
    if(length(dat_variables) == 1)
      dat_$Reponse <- dat[[dat_variables]]
  }

  # create plot function ----------------------------------------------------

  .pl_defaults <- list(
    pred = pred,
    pred_variables = pred_variables,
    pred_args = pred_args,
    dat_args = dat_args,
    land_args = land_args,
    world_args = world_args,
    graticule_args = graticule_args
  )

  pl <- function(pred = .pl_defaults$pred,
                 pred_variables = .pl_defaults$pred_variables,
                 pred_args = .pl_defaults$pred_args,
                 dat_args = .pl_defaults$dat_args,
                 land_args = .pl_defaults$land_args,
                 world_args = .pl_defaults$world_args,
                 graticule_args = .pl_defaults$graticule_args) {


    # create plot
    p <- plot_ly()

    # prepare prediction data
    if(!is.null(pred) & !is.null(pred_args)) {
      pred_ <- list()
      pred_$x <- cos(degrees2radians(pred[[dims[1]]]) + lon_shift) * cos(degrees2radians(pred[[dims[2]]]))
      pred_$y <- sin(degrees2radians(pred[[dims[1]]]) + lon_shift) * cos(degrees2radians(pred[[dims[2]]]))
      pred_$z <- sin(degrees2radians(pred[[dims[2]]]))
      stopifnot(length(pred_variables) == 1)
      pred_$Prediction <- pred[[pred_variables]]
      pred_ <- pred_ %>% lapply(matrix, nrow = length(unique(pred[[dims[1]]])))
      rm(pred)

      p <- p %>% add_surface(
        data = pred_,
        x = ~x, y = ~y, z = ~z,
        surfacecolor = ~Prediction,
        text = ~round(Prediction, 2),
        hoverinfo = "text",
        colorbar = list(title = "Temperature [°C]"),
        colorscale = colorscale,
        showlegend = FALSE,
        connectgaps = FALSE
      )
    }

    if(!is.null(dat) & !is.null(dat_args)) {
      p <- p %>% add_trace(
        data = dat_,
        x = ~x, y = ~y, z = ~z,
        color = I("darkred"),
        text = if(!is.null(dat_[["Response"]])) ~round(Response, 2),
        name = "Argo floats",
        hoverinfo = "text",
        showlegend = TRUE,
        marker = list(size = 2)
      )
    }
      p %>%
      add_maplayout(zoom = 1)
  }

  if(return_fun)
    return(pl) else
      pl()
}
