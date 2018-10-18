#' GetIndices
#'
#' Estimate indices of range dynamics from annual occupancy estimates
#' @param alpha Vector of alpha codes for species of interest
#' @return Data frame containing annual estimates of the following indices:
#' @return    avg.psi = proportion of area occupied (i.e., range size)
#' @return    s.lat = southern range limit
#' @return    s.core = southern core range limit
#' @return    n.lat = northern range limit
#' @return    n.core = northern core range limit
#' @return    w.lat = western range limit
#' @return    w.core = western core range limit
#' @return    e.lat = eastern range limit
#' @return    e.core = eastern core range limit
#' @return    avg.lat = occupancy-weighted mean breeding latitude
#' @return    avg.lon = occupancy-weighted mean breeding longitude
#' @export

GetIndices <- function(spp = NULL, alpha = NULL, path){
  if(!is.null(spp)){
    cores <- parallel::detectCores()
    if(length(spp) < cores) cores <- length(spp)
    doParallel::registerDoParallel(cores = cores)

    ### Run posterior predictive checks in parallel
    indices <- foreach::foreach(i = 1:length(spp), .combine = c,
                                    .packages = c("dplyr", "BayesCorrOcc")) %dopar%{
                                      occ <- readRDS(paste0(path, "/", spp[i], '/occ.rds'))
                                      dat <- readRDS(paste0(path, "/", spp[i], "/bbs_data.rds"))

                                      years <- seq(from = dat$start_year, to = dat$end_year)
                                      indices <- array(NA, dim = c(11, dim(occ$occ)[c(1,3)]))

                                      for(tt in 1:dat$nYears){
                                        indices[1,,tt] <- apply(occ$occ[,,tt], 1, mean)

                                        indices[2,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lat, limit = "south"))

                                        indices[3,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south"))

                                        indices[4,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lat, limit = "north"))

                                        indices[5,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north"))

                                        indices[6,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lon, limit = "west"))

                                        indices[7,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west"))

                                        indices[8,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lon, limit = "east"))

                                        indices[9,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east"))

                                        indices[10,,tt] <- apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE))

                                        indices[11,,tt] <- apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE))
                                      }


                                      saveRDS(indices, file = paste0(path, "/", spp[i], "/indices_post.rds"))

                                      return(spp[i])
                                    }
    return(indices)
  }

  if(!is.null(alpha)){
    occ <- readRDS(paste0(path, "/", alpha, '/occ.rds'))
    dat <- readRDS(paste0(path, "/", alpha, "/bbs_data.rds"))

    years <- seq(from = dat$start_year, to = dat$end_year)

    indices <- array(NA, dim = c(11, dim(occ$occ)[c(1,3)]))

    for(tt in 1:dat$nYears){
      indices[1,,tt] <- apply(occ$occ[,,tt], 1, mean)

      indices[2,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lat, limit = "south"))

      indices[3,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south"))

      indices[4,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lat, limit = "north"))

      indices[5,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north"))

      indices[6,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lon, limit = "west"))

      indices[7,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west"))

      indices[8,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.999, coord = occ$xy$lon, limit = "east"))

      indices[9,,tt] <- apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east"))

      indices[10,,tt] <- apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE))

      indices[11,,tt] <- apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE))
    }


    saveRDS(indices, file = paste0(path, "/", alpha, "/indices_post.rds"))
  }
}


#' range.limit
#'
#' Estimate range limit using cumulative probability method


range.limit <- function(cell.probs, prob, coord, limit){
  xy <- data.frame(x = cell.probs, y = coord)
  x <- dplyr::arrange(xy, y)$x/sum(xy$x, na.rm = TRUE)
  y <- dplyr::arrange(xy, y)$y
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]


    if(limit == "south"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), (1 - prob))$y
    }

    if(limit == "north"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), prob)$y
    }

  if(limit == "west"){
    lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), (1 - prob))$y
  }

  if(limit == "east"){
    lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), prob)$y
  }

  lim
}
