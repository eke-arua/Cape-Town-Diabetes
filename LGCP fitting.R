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
sim_list <- qs2::qs_read(here("500 factorial design simulations.qs2"))

#Fit the models
mod_list <- lapply(sim_list, function(X){

  mod <- kppm(X ~ cov_fac_im + offset(log(baseline_im)), "LGCP",
              model = "matern", nu = 0.5)


  return(mod)
})


#Extract
b1 <- sapply(mod_list, function(X){
  coefs <- coef(X[[1]])

})


b1 <- vector(mode = "list", length = length(mod_list))

for (i in seq_along(mod_list)){

  coefs <- coef(mod_list[[i]])

  b1[[i]] <- (coefs[2])

}

