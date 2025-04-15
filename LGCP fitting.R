library(spatstat)
library(here)
library(terra)
library(scampr)
library(qs2)

#Load the covariates
pop_density <- terra::rast("Cape Town population density.tif")
pop_density <- project(pop_density, "EPSG:22234")

facility_dist <- rast("Distance to healthcare facility.tif")
facility_dist <- project(facility_dist, "EPSG:22234")

#Population density as the baseline (offset)
baseline <- classify(scale(pop_density), cbind(-Inf, 0, 0.0001), include.lowest = TRUE)

#Load the point simulated point pattern data
mod_list <- qs2::qs_read(here("500 factorial design simulations.qs2"))


