#' GetBioVars
#'
#' Extract scaled annual bioclim variable values for each BBS route
#' @param alpha alpha code for species of interest
#' @param index Integer vector containing the bioclim variables of interest
#' @param ind_name Character vector containing abbreviated name of the bioclim variables of interest
#' @return A .csv files containing the following fields:
#' @return   routeID The unique 8 digit route ID for each route
#' @return   Latitude The latitude for the route
#' @return   Longitude The longitude for the route
#' @return   Ind_Year The value for biovar[Ind] in year Year
#' @export

GetBioVars <- function(alpha, index = c(1, 2, 8, 12, 18), path,
                       ind_name = c("tmp", "dtr", "Twet", "Prec", "Pwarm")){

  dat <- readRDS(paste0(path, "/", alpha, "/bbs_data.rds"))

  rxy <- data.frame(Latitude = dat$lat, Longitude = dat$lon)

  xy <- dplyr::select(rxy, Longitude, Latitude)

  bbs_years <- seq(from = dat$start_year, to = dat$end_year)

  if(sum(bbs_years %in% seq(from = 1971, to = 2014)) < length(bbs_years)) stop("Count data contains years with no climate data")

  for (jj in 1:length(index)) {
    for (ii in 1:dat$nYears) {
      out	<- raster::extract(BBS.SDM::NA_biovars[[ii]][[index[jj]]], xy)
      varname <- paste(ind_name[[jj]], bbs_years[[ii]], sep = "_")
      rxy[[varname]] <- out
    }
  }

  ## For locations w/ NA for climate variables, fill in mean of 8 neighboring cells
  problem_routes <- which(is.na(rxy[, ncol(rxy)]))
  first_climate <- grep(ind_name[1], colnames(rxy))[1]

    if(length(problem_routes) > 0){
      ind <- 0
      min.adj <- -0.5
      max.adj <- 0.5
      while(ind < 1){
        problem_xy <- dplyr::select(dplyr::slice(rxy, problem_routes), Longitude, Latitude)

        fix_routes <- problem_out <- problem_nearby <- NULL
        adj <- c(min.adj,0,max.adj)

        for (jj in 1:length(index)) {
          for (ii in 1:dat$nYears) {
            for (i_lat in 1:3) {
              for (i_lon in 1:3) {
                mjc	<- raster::extract(BBS.SDM::NA_biovars[[ii]][[index[jj]]],
                                       problem_xy + matrix(c(adj[i_lon], adj[i_lat]),
                                                           nrow = length(problem_routes), ncol=2), byrow=T)
                problem_nearby	<- cbind(problem_nearby, mjc)
              }
            }
            problem_out <- rowMeans(problem_nearby, na.rm = TRUE)
            fix_routes <- cbind(fix_routes, problem_out)
          }
        }
        rxy[problem_routes, first_climate:ncol(rxy)] <- fix_routes
        problem_routes2 <- which(is.na(rxy[, ncol(rxy)]))
        if(length(problem_routes2) == 0){
          ind <- 1
          problem_routes <- problem_routes2
        }else{
          problem_routes <- problem_routes2
          min.adj <- min.adj - 0.5
          max.adj <- max.adj + 0.5
        }
      }
    }


  ### Center and scale climate variables
  clim_scale <- matrix(99, nrow = length(ind_name), ncol = 2, dimnames = list(ind_name, c("mean", "sd")))

  for(ii in seq_along(ind_name)){
    clim_scale[ii, "mean"] <- mean(as.matrix(dplyr::select(rxy, grep(ind_name[ii], names(rxy)))), na.rm = TRUE)
    clim_scale[ii, "sd"] <- sd(as.matrix(dplyr::select(rxy, grep(ind_name[ii], names(rxy)))), na.rm = TRUE)

    rxy[,grep(ind_name[ii], names(rxy))] <- (rxy[,grep(ind_name[ii], names(rxy))] - clim_scale[ii, "mean"]) / clim_scale[ii, "sd"]
  }

  write.csv(clim_scale,
            paste(path, alpha, "clim_scale.csv", sep = "/"),
            row.names = FALSE)

  ### Add squared climate variables
  sq_rxy <- dplyr::as_data_frame(rxy[, first_climate:ncol(rxy)] ^ 2)
  colnames(sq_rxy) <- paste("sq", colnames(sq_rxy), sep = "_")
  rxy <- dplyr::bind_cols(rxy, sq_rxy)
  rxy <- dplyr::select(rxy, -Latitude, -Longitude)

  biovars <- NULL
  for(i in 1:dat$nYears){
    tmp <- rxy[,grepl(bbs_years[i], colnames(rxy))]
    biovars <- abind::abind(biovars, tmp, along = 3)
  }

  saveRDS(biovars, file = paste(path, alpha, "biovars.rds", sep = "/"))
}

