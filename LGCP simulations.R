library(parallel)
library(sf)
library(terra)
library(spatstat)
library(scampr)
library(dplyr)

set.seed(123)

#Load boundary
cape_town <- st_read(here("sa_census_2011_cape_town_sp_gis_data",
                          "SA Census 2011 Cape Town SP",
                          "SP_Cape_Town_2011_Apr2013.shp")) |>
  st_transform(crs = 22234)

cape_town <- cape_town[cape_town$SP_NAME != "Robben Island SP", ]
cape_town <- st_make_valid(cape_town)

#Load the covariates
pop_density <- terra::rast("Cape Town population density.tif") |>
  project("EPSG:22234")

facility_dist <- rast("Distance to healthcare facility.tif") |>
  project("EPSG:22234")

facility_dist <- facility_dist/10000

#Population density as the baseline (offset). We replace zeros with an arbitrary small number
baseline <- classify(scale(pop_density), cbind(-Inf, 0, 0.0001), include.lowest = TRUE)

#Load the parameters
parameter_space <- list(beta_0 = c(-10, -12, -14),
                      beta_1 = c(-0.1, -0.5, -1),
                      scale = c(100, 400, 800),
                      var = c(1, 1.5, 2)) |>
  expand.grid()

source("parameter sampling.R")

parameters <- parameter_sample(df = parameter_space, size = 300)

#Simulate LGCPs with matern correlation from the given parameters.
sim_list <- mclapply(seq_len(nrow(parameters)), function(x) {

  params <- unlist(parameters[x, ])  # Extract row as vector

  eta <- log(baseline) + params[1] + (params[2] * scale(facility_dist))
  intensity <- exp(eta)

  # Convert to an image format
  intensity_df <- as.data.frame(intensity, xy = TRUE)
  intensity_im <- vec2im(intensity_df$sum, intensity_df$x, intensity_df$y)

  # Simulate the log-Gaussian Cox process
  matern_sim <- spatstat.random::rLGCP(
    model = "matern",
    mu = log(intensity_im),
    var = params[4],
    scale = params[3],
    nu = 0.5,
    saveLambda = FALSE
  )

  return(matern_sim)
}, mc.cores = 5)

qs2::qs_save(sim_list, "Simulations.rds")

parameters$n_events <- sapply(sim_list, \(X){
  X$n
})

parameters$index <- 1:nrow(parameters)

mod <- kppm(sim_list[[65]] ~ cov_fac_im +  offset(log(baseline_im)), "LGCP",
model = "matern", nu = 0.5)

library(qs2)

mod_list <- qs2::qs_read("500 factorial design models.qs2")

install.packages("~/Downloads/maptools_1.1-8(1).tar.gz", repos = NULL, type = "source")
