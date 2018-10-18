#' OccSummary
#'
#' Estimate mean & CI psi for all raster cells
#' @param spp Vector containing alpha codes for all species of interest
#' @param alpha Four letter alpha code for species of interest
#' @return Data frame containing annual estimates of the following indices:
#' @export

OccSummary <- function(spp = NULL, alpha = NULL, path){
  if(!is.null(spp)){
    cores <- parallel::detectCores()
    if(length(spp) < cores) cores <- length(spp)
    doParallel::registerDoParallel(cores = cores)

    ### Run posterior predictive checks in parallel
    occ_summary <- foreach::foreach(i = 1:length(spp), .combine = c,
                                .packages = c("dplyr", "BayesCorrOcc")) %dopar%{
                                  occ <- readRDS(paste0(path, "/", spp[i], '/occ.rds'))
                                  dat <- readRDS(paste0(path, "/", spp[i], "/bbs_data.rds"))

                                  years <- seq(from = dat$start_year, to = dat$end_year)


                                  Psi <- LCI <- UCI <- NULL
                                  for(tt in 1:dat$nYears){
                                    tPsi <- apply(occ$occ[,,tt], 2, function(x) mean(x))
                                    tLCI <- apply(occ$occ[,,tt], 2, function(x) quantile(x, probs = 0.025))
                                    tUCI <- apply(occ$occ[,,tt], 2, function(x) quantile(x, probs = 0.975))

                                    Psi <- c(Psi, tPsi)
                                    LCI <- c(LCI, tLCI)
                                    UCI <- c(UCI, tUCI)
                                  }

                                  psi <- data.frame(Year = rep(years, each = dim(occ$occ)[2]),
                                                    Latitude = rep(occ$xy$lat, dat$nYears),
                                                    Longitude = rep(occ$xy$lon, dat$nYears),
                                                    Psi = Psi, LCI = LCI, UCI = UCI)

                                  saveRDS(psi, file = paste0(path, "/", spp[i], "/psi.rds"))
                                  return(spp[i])
                                }
    return(occ_summary)
  }

  if(!is.null(alpha)){
    occ <- readRDS(paste0(path, "/", alpha, '/occ.rds'))
    dat <- readRDS(paste0(path, "/", alpha, "/bbs_data.rds"))

    years <- seq(from = dat$start_year, to = dat$end_year)


    Psi <- LCI <- UCI <- NULL
    for(tt in 1:dat$nYears){
      tPsi <- apply(occ$occ[,,tt], 2, function(x) mean(x))
      tLCI <- apply(occ$occ[,,tt], 2, function(x) quantile(x, probs = 0.025))
      tUCI <- apply(occ$occ[,,tt], 2, function(x) quantile(x, probs = 0.975))

      Psi <- c(Psi, tPsi)
      LCI <- c(LCI, tLCI)
      UCI <- c(UCI, tUCI)
    }

    psi <- data.frame(Year = rep(years, each = dim(occ$occ)[2]),
                      Latitude = rep(occ$xy$lat, dat$nYears),
                      Longitude = rep(occ$xy$lon, dat$nYears),
                      Psi = Psi, LCI = LCI, UCI = UCI)

    saveRDS(psi, file = paste0(path, "/", alpha, "/psi.rds"))
  }
}
