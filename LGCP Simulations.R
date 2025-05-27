library(parallel)
library(sf)
library(terra)
library(spatstat)
library(scampr)
library(dplyr)
library(qs2)
library(here)

#Load the covariates
pop_density <- terra::rast("Cape Town population density.tif")
pop_density <- project(pop_density, "EPSG:22234")

facility_dist <- rast("Distance to healthcare facility.tif")/10000
facility_dist <- project(facility_dist, "EPSG:22234")

#Population density as the baseline (offset)
baseline <- classify(pop_density, cbind(-Inf, 0, 0.0001), include.lowest = TRUE)

#Load the parameters
parameter_space <- list(beta_0 = c(-13, -14, -15),
                        beta_1 = c(-0.1, -0.5, -1),
                        scale = c(100, 400, 800),
                        var = c(1, 1.5, 2)) |>
  expand.grid()

set.seed(123)

#Sample the parameters for the simulation from the parameter space
source(here("Helper functions", "parameter sampling.R"))

parameters <- parameter_sample(df = parameter_space, size = 500)

#Simulate LGCPs with matern correlation from the given parameters.
sim_list <- mclapply(seq_len(nrow(parameters)), function(i) {

  params <- unlist(parameters[i, ])  # Extract row as vector

  eta <- log(baseline) + params[1] + (params[2] * facility_dist)
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
}, mc.cores = 8)

qs_save(sim_list, "500 factorial design simulations.qs2")

parameters$n_events <- sapply(sim_list, function(X){
  X$n
})
parameters$index <- seq(1:nrow(parameters))

qs_save(parameters, "500 Parameters.qs2")

