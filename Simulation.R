library(here)
library(sf)
library(terra)
library(spatstat)
library(scampr)
library(dplyr)

cape_town <- st_read(here("sa_census_2011_cape_town_sp_gis_data",
                          "SA Census 2011 Cape Town SP",
                          "SP_Cape_Town_2011_Apr2013.shp")) |>
  st_transform(crs = 22234)

#Remove robben island
cape_town <- cape_town[cape_town$SP_NAME != "Robben Island SP", ]

#Load the covariates
pop_density <- terra::rast("Cape Town population density.tif")
facility_dist <- rast("Distance to healthcare facility.tif")
park_dist <- rast("Distance to parks.tif")

#Reduce spatial resolution
pop_density <- project(pop_density, "EPSG:22234")
agg_factor <- 1000 / res(pop_density)
pop_density_1km <- aggregate(pop_density, fact= (agg_factor), fun = "sum") # Use aggregate for downscaling

facility_dist <- project(facility_dist, "EPSG:22234")
fac_dist_1km <- aggregate(facility_dist, fact= (agg_factor), fun = "mean")

park_dist <- project(park_dist, "EPSG:22234")
park_dist_1km <- aggregate(park_dist, fact= (agg_factor), fun = "mean")

#Observation window
matrx <- as.matrix(pop_density_1km, wide = TRUE)
matrx_flipped <- matrx[nrow(matrx):1, ]

#image <- as.im(matrx_flipped)

extent_image <- ext(pop_density_1km)

mask_matrix_binary <- ifelse(is.nan(matrx_flipped), FALSE, TRUE)
mask_matrix_binary <- !is.nan(matrx_flipped)



t_win <- owin(mask = mask_matrix_binary,  xrange = c(extent_image$xmin, extent_image$xmax),
              yrange = c(extent_image$ymin, extent_image$ymax))

#log intensity

#Beta coefficients
b0 <- 0.00000002
b1 <- 0.000005
b2 <- 0.00000001
b3 <- 0.00000002

#Rescale the variables
# pop <- pop_density_1km/10000
# fac <- fac_dist_1km/25000
# park <- park_dist_1km/15000

eta <-  log(pop_density_1km) + (b2*park_dist_1km) + (b3*fac_dist_1km)
#eta_scale <- b0 + (b1*scale(pop_density_1km)) + (b2*scale(park_dist_1km)) + (b2*scale(fac_dist_1km))

eta_df <- as.data.frame(eta, xy = TRUE)

log_intensity = vec2im(exp(eta_df$sum), eta_df$x, eta_df$y)
#log_intensity <- as.im(log_intensity, W = t_win)

#rescale
Y <- spatstat.geom::rescale(log_intensity, 20, "metres")

test <- rpoispp(lambda = (Y))
test
plot(test)
plot(density(test, sigma = 15))

persp(density(test, sigma = 3000), theta = 345, phi = 55)
plot(density(test, sigma = 3000))
plot(log_intensity_pos)

plot(test, cex = 0.01, alpha = 0.01)

matern_sim = spatstat.random::rLGCP(model = "matern",
                                    mu = (Y),
                                    var = 1, scale = 300, nu = 0.5,
                                    saveLambda = TRUE)

matern_sim

exp_sim = rLGCP(model = "exp",
                mu = (log_intensity_m),
                var = 1, scale = 300,
                saveLambda = TRUE)
