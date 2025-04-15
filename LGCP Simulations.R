library(parallel)
library(sf)
library(terra)
library(spatstat)
library(scampr)
library(dplyr)

# #Load boundary
# cape_town <- st_read(here("sa_census_2011_cape_town_sp_gis_data",
#                           "SA Census 2011 Cape Town SP",
#                           "SP_Cape_Town_2011_Apr2013.shp")) |>
#   st_transform(crs = 22234)

# cape_town <- cape_town[cape_town$SP_NAME != "Robben Island SP", ]
# cape_town <- st_make_valid(cape_town)

#Load the covariates
pop_density <- terra::rast("Cape Town population density.tif")
pop_density <- project(pop_density, "EPSG:22234")

facility_dist <- rast("Distance to healthcare facility.tif")
facility_dist <- project(facility_dist, "EPSG:22234")

#Population density as the baseline (offset)
baseline <- classify(scale(pop_density), cbind(-Inf, 0, 0.0001), include.lowest = TRUE)

#Load the parameters
parameter_space <- list(beta_0 = c(-12, -13, -14),
                        beta_1 = c(-0.1, -0.5, -1),
                        scale = c(100, 400, 800),
                        var = c(1, 1.5, 2)) |>
  expand.grid()

set.seed(123)

source("parameter sampling.R")

parameters <- parameter_sample(df = parameter_space, size = 500)

saveRDS(parameters, "500 Parameters.rds")

#Simulate LGCPs with matern correlation from the given parameters.
sim_list <- mclapply(seq_len(nrow(parameters)), function(i) {

  params <- unlist(parameters[i, ])  # Extract row as vector

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
}, mc.cores = 7)

saveRDS(sim_list, "500 factorial design simulations.rds")

qs_save(sim_list, "500 factorial design simulations.qs2")


n_events <- sapply(sim_list, function(X){
  X$n
})

parameters$n_events <- n_events
parameters$index <- seq(1:nrow(parameters))

#Fit the models
cov_fac <- as.data.frame(scale(facility_dist), xy = TRUE)
cov_fac_im <- vec2im(cov_fac$sum, cov_fac$x, cov_fac$y)

baseline_im <- (as.data.frame((baseline), xy = TRUE))
baseline_im <- vec2im(baseline_im$sum, baseline_im$x, baseline_im$y)


#sim_list <- readRDS("300 factorial design simulations.rds")

mod_list <- lapply(sim_list, function(X){

  mod <- kppm(X ~ cov_fac_im + offset(log(baseline_im)), "LGCP",
              model = "matern", nu = 0.5)


  return(mod)
})

saveRDS(mod_list, "/home/shared/shared-data/Spatial Sim results Eke/500 factorial design models.rds")


qs_save(mod_list[[1]], "/home/shared/shared-data/Spatial Sim results Eke/500 factorial design models.qs2")

mod <- kppm(sim_list[[35]] ~ cov_fac_im + offset(log(baseline_im)), "LGCP",
            model = "matern", nu = 0.5)

#
# library(spatstat)
# set.seed(123)
#
# # Define number of bootstrap samples
# n_boot <- 100
#
# # Store bootstrap estimates
# boot_results <- matrix(NA, n_boot, length(mod$modelpar))
#
# for (i in 1:n_boot) {
#   # Simulate a new point pattern from the fitted model
#   sim_data <- simulate(mod, nsim = 1)[[1]]
#   print(1)
#   # Refit LGCP model
#   fit_boot <- kppm(sim_data ~ 1, clusters = "LGCP")
#   print(2)
#   # Store estimated hyperparameters (variance and scale)
#   boot_results[i, ] <- mod$modelpar
#   print(3)
# }
#
# # Compute confidence intervals (e.g., 95% percentile-based)
# ci_lower <- apply(boot_results, 2, quantile, probs = 0.025)
# ci_upper <- apply(boot_results, 2, quantile, probs = 0.975)
#
# # Display results
# data.frame(Parameter = names(fit$modelpar),
#            Estimate = fit$modelpar,
#            CI_Lower = ci_lower,
#            CI_Upper = ci_upper)
#

b1 <- sapply(mod_list, function(X){
  coefs <- coef(X[[1]])

})


b1 <- vector(mode = "list", length = length(mod_list))

for (i in seq_along(mod_list)){

  coefs <- coef(mod_list[[i]])

  b1[[i]] <- (coefs[2])

}



b1 <- unlist(b1)
beta_1 <- unlist(parameters[2])

bias_b1 <- mean(b1 - beta_1)

b0 <- vector(mode = "list", length = length(mod_list))

for (i in seq_along(mod_list)){

  coefs <- coef(mod_list[[i]])

  b0[[i]] <- (coefs[1])

}

b0 <- unlist(b0)

beta_0 <- unlist(parameters[1])

bias_b0 <- mean(b0 - beta_0)

sizes <- unlist(sizes)
sizes_mb <- sizes/1024^2
