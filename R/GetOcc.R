#' GetOccProb
#'
#' Predict annual probability of occupancy using parameter estimates from top model
#' @param spp Vector containing alpha codes for all species of interest
#' @param alpha Four letter alpha code for species of interest
#' @export

GetOccProb <- function(alpha = NULL, path){
  if(is.null(alpha)){
    spp_list <- read.csv(path)
    spp <- as.character(spp_list$spp)

    ### Register cores
    cores <- parallel::detectCores()
    if(length(spp) < cores) cores <- length(spp)
    doParallel::registerDoParallel(cores = cores)

    ### Run posterior predictive checks in parallel
    occ_run <- foreach::foreach(i = 1:length(spp), .combine = c,
                                .packages = c("dplyr", "BayesCorrOcc")) %dopar%{
                                  ## Read count data
                                  dat <-  readRDS(paste0("inst/output/", spp[i], "/bbs_data.rds"))

                                  ## Read jags output
                                  sim_list <- readRDS(paste0("inst/output/", spp[i], "/jags_fit.rds"))

                                  years <- seq(from = dat$start_year, to = dat$end_year)

                                  ## Get climate values for all cells/years
                                  covs <- raster.to.array(alpha = spp[i], years)

                                  ## Get model matrix for s(lat, lon) of raster cells
                                  xy <- data.frame(lon = covs$xy$x, lat = covs$xy$y)
                                  org.data <- data.frame(z = rep(1, length(dat$lat)), lon = dat$lon, lat = dat$lat)
                                  org.mod <- mgcv::gam(z ~ s(lon, lat, k = 60, bs = 'ds', m = c(1, 0.5)), data = org.data, family = "binomial")
                                  X <- predict(object = org.mod, type = "lpmatrix", newdata = xy)

                                  ## For each posterior sample, estimate occupancy for each cell
                                  r.psi <- array(0, dim = c(dim(sim_list$sims.list$xpsi)[1], dim(covs$climate)[1], length(years)))

                                  for(ii in 1:dim(sim_list$sims.list$xpsi)[1]){
                                    for (tt in 1:length(years)) {
                                      lpsi <- X %*% sim_list$sims.list$b[ii,,tt] +
                                        covs$climate[,,tt] %*% (sim_list$sims.list$g[ii,] * sim_list$sims.list$betaT[ii,])
                                      r.psi[ii, , tt] <- plogis(lpsi)
                                    }
                                  }

                                  occ <- list(occ = r.psi, xy = xy)
                                  saveRDS(occ, file = paste0("inst/output/", spp[i], "/occ.rds"))
                                  return(spp[i])
                                }
    return(occ_run)
  }else{
    ## Read count data
    dat <-  readRDS(paste(path, alpha, "bbs_data.rds", sep = "/"))

    ## Read jags output
    sim_list <- readRDS(paste(path, alpha, "jags_fit.rds", sep = "/"))

    years <- seq(from = dat$start_year, to = dat$end_year)

    ## Get climate values for all cells/years
    covs <- raster.to.array(alpha, years)

    ## Get model matrix for s(lat, lon) of raster cells
    xy <- data.frame(lon = covs$xy$x, lat = covs$xy$y)
    org.data <- data.frame(z = rep(1, length(dat$lat)), lon = dat$lon, lat = dat$lat)
    org.mod <- mgcv::gam(z ~ s(lon, lat, k = 60, bs = 'ds', m = c(1, 0.5)), data = org.data, family = "binomial")
    X <- predict(object = org.mod, type = "lpmatrix", newdata = xy)

    ## For each posterior sample, estimate occupancy for each cell
    r.psi <- array(0, dim = c(dim(sim_list$sims.list$xpsi)[1], dim(covs$climate)[1], length(years)))

    for(ii in 1:dim(sim_list$sims.list$xpsi)[1]){
      for (tt in 1:length(years)) {
        lpsi <- X %*% sim_list$sims.list$b[ii,,tt] +
          covs$climate[,,tt] %*% (sim_list$sims.list$g[ii,] * sim_list$sims.list$betaT[ii,])
        r.psi[ii, , tt] <- plogis(lpsi)
      }
    }

    occ <- list(occ = r.psi, xy = xy)
    saveRDS(occ, file =paste(path, alpha, "occ.rds", sep = "/"))
  }

}


#' raster.to.array
#'
#' Convert bioclim rasters in array w/ dim n.cell x n.vars x n.years

raster.to.array <- function(alpha, years, path) {
  ## Read count data
  dat <-  readRDS(paste(path, alpha, "bbs_data.rds", sep = "/"))

  index <- c(1,2,8,12,18)
  scale.values <- read.csv(paste(path, alpha, "clim_scale.csv", sep = "/"))

  xy.values <- data.frame(Latitude = dat$lat, Longitude = dat$lon)
  mu.lat <- mean(xy.values$Latitude)
  sd.lat <- sd(xy.values$Latitude)

  mu.lon <- mean(xy.values$Longitude)
  sd.lon <- sd(xy.values$Longitude)



  for (ii in seq_along(years)){
    # get climate data within masked area
    for (jj in seq_along(index)) {
      range <- raster::mask(BayesCorrOcc::NA_biovars[[paste0("biovars",as.character(years[ii]))]][[index[jj]]], dat$hull)
      assign(paste0("bio.", index[jj]), range)    # mask the climate data
      assign(paste0("all.values.",index[jj]), raster::getValues(get(paste0("bio.",index[jj]))))                                                 # extract climate values
      assign(paste0("values.",index[jj]), get(paste0("all.values.",index[jj]))[!is.na(get(paste0("all.values.",index[jj])))])                   # remove NAs
      assign(paste0("center.values.",index[jj]), (get(paste0("values.",index[jj])) - scale.values[jj,"mean"])/scale.values[jj,"sd"])     # center and scale
      if (ii+jj==2)  climate <- array(NA, dim=c(length(get(paste0("values.",index[jj]))), 10, length(years)))
      climate[, jj, ii] <- get(paste0("center.values.",index[jj]))
    }  # end jj loop
      climate[, (length(index) + 1):(length(index)*2), ii] <- climate[,1:length(index),ii]^2
  } # end ii loop


  dimnames(climate) <- list(1:(dim(climate)[1]),
                            c("tmp","dtr","Twet","Prec","Pwarm","sq_tmp","sq_dtr","sq_Twet","sq_Prec","sq_Pwarm"),
                            c(years))


  ## Scale lat/long
  xy <- as.data.frame(raster::rasterToPoints(bio.1))
  xy <- dplyr::select(xy, -bio1)

  raster_covs <- list(climate = climate, xy = xy)
  return(raster_covs)
}   # end function
