library(spdep)
library(terra)
library(here)
library(qs2)
library(parallel)
library(INLA)
library(dplyr)

#Load boundaries for subplace
cape_town_sp <- st_read(here("sa_census_2011_cape_town_sp_gis_data",
                             "SA Census 2011 Cape Town SP",
                             "SP_Cape_Town_2011_Apr2013.shp")) |>
  st_transform(crs = 22234)

cape_town_sp <- cape_town_sp[cape_town_sp$SP_NAME != "Robben Island SP", ]
cape_town_sp <- st_make_valid(cape_town_sp)

#Create boundaries for main-place
cape_town_mp <- cape_town_sp %>%
  group_by(MP_NAME) %>%
  summarise(st_union(geometry))

#Load the spatially continuous covariates
pop_density <- terra::rast("Cape Town population density.tif") |>
  project("EPSG:22234")

facility_dist <- rast("Distance to healthcare facility.tif") |>
  project("EPSG:22234")

#convert it to 10's of kilometers
facility_dist <- facility_dist/10000

#Population density as the baseline (offset)
baseline <- classify(scale(pop_density), cbind(-Inf, 0, 0.0001), include.lowest = TRUE)

#Aggregate the covariates according to the boundary for sub-place
boundary_sp <- vect(cape_town_sp)

cape_town_sp$pop_density <- extract(baseline, boundary_sp, fun = sum, na.rm = TRUE)[,2]
cape_town_sp$fac_dist <- extract(facility_dist, boundary_sp, fun = mean, na.rm = TRUE)[,2]

#Aggregate the covariates according to the boundary for main-place
boundary_mp <- vect(cape_town_mp)

cape_town_mp$pop_density <- extract(baseline, boundary_mp, fun = sum, na.rm = TRUE)[,2]
cape_town_mp$fac_dist <- extract(facility_dist, boundary_mp, fun = mean, na.rm = TRUE)[,2]

#Load the parameters the data was simulated from
parameters <- qs_read(here("500 Parameters.qs2"))

#Load the simulated data
sim_list <- qs2::qs_read(here("500 factorial design simulations.qs2"))

#### aggregate the simulated point pattern data to the areal level (sub place ####
#and main place)
source(here("Helper functions", "Aggregate counts.R"))

#For sub-place
sim_aggregated_list_sp <- lapply(sim_list, function(X, boundary = cape_town_sp){

  points <- as.data.frame(X)
  points_sf <- st_as_sf(points, crs = 22234, coords = c("x", "y"))

  ct_boundary <- boundary
  boundary$counts <- agg_counts(boundary, points_sf)
  boundary$ID <- seq(1, nrow(boundary))
  return(boundary)

})

#For main-place
sim_aggregated_list_mp <- lapply(sim_list, function(X, boundary = cape_town_mp){

  points <- as.data.frame(X)
  points_sf <- st_as_sf(points, crs = 22234, coords = c("x", "y"))

  ct_boundary <- boundary
  boundary$counts <- agg_counts(boundary, points_sf)
  boundary$ID <- seq(1, nrow(boundary))
  return(boundary)

})

#### Fit BYM models for subplace boundaries ####
W_nb_sp <- poly2nb(cape_town_sp)
W_mat_sp <- nb2mat(W_nb_sp, style="B")

g_sp <- INLA::inla.read.graph(W_mat_sp) #incase INLA is needed

#Alternate priors
prior_list <-  list(
  prec = list(prior = "pc.prec", param = c(0.25 , 0.01)),
  phi = list(prior = "pc",    param = c(0.5, 0.5)))

bym_models_sp <- lapply(sim_aggregated_list_sp, function(X) {
  tryCatch({
    inla(
      counts ~ fac_dist + f(ID, model = "bym2", graph = g_sp),
      family = "poisson",
      data = X, E = pop_density,
      control.predictor = list(compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE)
    )
  }, error = function(e1) {
    message("INLA failed: ", conditionMessage(e1), " -- trying again with different priors")

    tryCatch({
      inla(
        counts ~ fac_dist + f(ID, model = "bym2", graph = g_sp, hyper = prior_list),
        family = "poisson",
        data = X, E = pop_density,
        control.predictor = list(compute = TRUE),
        control.compute = list(dic = TRUE, waic = TRUE)
      )
    }, error = function(e2) {
      message("INLA failed with different priors: ", conditionMessage(e2),
              " -- trying with different BYM specification")

      tryCatch({
        inla(
          counts ~ fac_dist + f(ID, model = "bym", graph = g_sp),
          family = "poisson",
          data = X, E = pop_density,
          control.predictor = list(compute = TRUE),
          control.compute = list(dic = TRUE, waic = TRUE)
        )
      }, error = function(e3) {
        message("INLA failed again: ", conditionMessage(e3))
        return(NA)
      })

    })
  })
})

form <- counts ~ fac_dist + offset(log(pop_density))

bym_models_sp_alt <- lapply(sim_aggregated_list_sp, function(X) {

  tryCatch({
    bym_mod <- S.CARbym(
      form, family = "poisson",
      data = X, n.sample = 15000, burnin = 5000, thin = 50,
      W = W_mat_sp)

    return(bym_mod)
  },
  error = function(e) {
    message("Error: ", conditionMessage(e))
    return(NA)
  })
})

bym_mod <- S.CARbym(
  form, family = "poisson",
  data = sim_aggregated_list_sp[[1]], n.sample = 15000, burnin = 5000, thin = 50,
  W = W_mat_sp)

#### Fit BYM models for main-place boundaries ####
W_nb_mp <- poly2nb(cape_town_mp)
W_mat_mp <- nb2mat(W_nb_mp, style = "B")

g_mp <- INLA::inla.read.graph(W_mat_mp)

bym_models_mp <- lapply(sim_aggregated_list_mp, function(X) {

  tryCatch({
    bym_mod <- inla(
      counts ~ fac_dist + f(ID, model = "bym2", graph = g_mp,
                            hyper = prior_list), family = "poisson",
      data = X, E = pop_density)

    return(bym_mod)
  },
  error = function(e) {
    message("Error: ", conditionMessage(e))
    return(NA)
  })

})

#Extract the estimated fixed effects parameters (betas)

#### Bias ####
betas_list_sp <- lapply(bym_models_sp, function(X){

  b0_sp <- X$summary.fixed[1, 1,]
  b1_sp <- X$summary.fixed[2, 1]

  return(data.frame(cbind(b0_sp, b1_sp)))
})

#Double checking with CARBayes
# betas_list_sp_alt <- lapply(bym_models_sp_alt, function(X){
#
#   b0_sp <- X$summary.results[1, 1]
#   b1_sp <- X$summary.results[2, 1]
#
#   return(data.frame(cbind(b0_sp, b1_sp)))
# })


betas_list_mp <- lapply(bym_models_mp, function(X){

  b0_mp <- X$summary.fixed[1 , 1]
  b1_mp <- X$summary.fixed[2, 1]
  return(data.frame(cbind(b0_mp, b1_mp)))
})

betas_mp <- Reduce(f= rbind, betas_list_mp)
betas_sp <- Reduce(f= rbind, betas_list_sp)
betas_sp_alt <- Reduce(f = rbind, betas_list_sp_alt)

#Bias
bias_b0_sp <- mean(parameters$beta_0 - betas_sp$b0_sp)
bias_b1_sp <- mean(parameters$beta_1 - betas_sp$b1_sp)

# bias_b0_sp_alt <- mean(parameters$beta_0 - betas_sp_alt$b0_sp)
# bias_b1_sp_alt <- mean(parameters$beta_1 - betas_sp_alt$b1_sp)

bias_b0_mp <- mean(parameters$beta_0 - betas_mp$b0_mp)
bias_b1_mp <- mean(parameters$beta_1 - betas_mp$b1_mp)


#RMSE
rmse_b0_mp <- sqrt(mean((parameters$beta_0 - betas_mp$b0_mp)^2))
rmse_b1_mp <- sqrt(mean((parameters$beta_1 - betas_mp$b1_mp)^2))

rmse_b0_sp <- sqrt(mean((parameters$beta_0 - betas_sp$b0_sp)^2))
rmse_b1_sp <- sqrt(mean((parameters$beta_1 - betas_sp$b1_sp)^2))

#Coverage
ci_b0_sp <- lapply(bym_models_sp, function(X){

  as.data.frame(cbind(low_ci_sp = X$summary.fixed$`0.025quant`[1],
               high_ci_sp = X$summary.fixed$`0.975quant`[1]))
}) |>
  Reduce(f = rbind)

cov_b0 <- sum(between(parameters$beta_0, ci_b0_sp$low_ci_sp, ci_b0_sp$high_ci_sp))


ci_b1_sp <- lapply(bym_models_sp, function(X){

  as.data.frame(cbind(low_ci_sp = X$summary.fixed$`0.025quant`[2],
                      high_ci_sp = X$summary.fixed$`0.975quant`[2]))
}) |>
  Reduce(f = rbind)

cov_b1 <- mean(between(parameters$beta_1, ci_b1_sp$low_ci_sp, ci_b1_sp$high_ci_sp))

