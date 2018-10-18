#' GetCorrData
#'
#' Filter raw BBS counts by species, add 0 counts, and latitude/longitude
#' @param bbs List containing the raw BBS counts, weather, and route info obtained from the function `get_BBSn()`
#' @param alpha 4-letter alpha code for the species of interest
#' @param common Common name for the species of interest
#' @param start.year Optional start year (default is 1966)
#' @param end.year Optional end year (default is most recent year of available BBS data)
#' @param countrynum Optional numeric vector containing the countries of interest (840 = US, 124 = Canada, 484 = Mexico)
#' @param statenum Optional numeric vector containing the states of interest (see BBS website for code values)
#' @param strata Optional numeric vector containing the BBS strata of interest (see BBS website for code values)
#' @param bcr Optional numeric vector containing the BCRs of interest (see BBS website for code values)
#' @param occ Counts or presence/absence? (Default == presence/absence)
#' @param lon Optional numeric vector containing maximum and minimum longitude of routes to keep
#' @return A .csv file containing the following fields:
#' @return   routeID The unique 8 digit route ID for each route
#' @return   Year The year that the count was conducted
#' @return   aou The numeric code for the focal species
#' @return   countN Total individuals of the focal species recorded at stop N (50-stop data) or N-9:N (10-stop data)
#' @return   Latitude The latitude for the route
#' @return   Longitude The longitude for the route
#' @return   Stratum The BBS stratum for the route
#' @return   BCR The Bird Conservation Region for the route
#' @export

GetCorrData <- function(bbs_raw = bbs, path, alpha = NULL, common = NULL,
                         start.year = 1966, end.year = 2016,
                         statenum = NULL, countrynum = NULL, strata = NULL, bcr = NULL, buffer = 2,
                         occ = TRUE, lon = NULL){

  if(is.null(alpha) & is.null(common))stop("Please provide either the alpha code or common name of the species of interest")

  if(!is.null(alpha)){
    AOU <- GetAOU(alpha = alpha)
  }else{
    AOU <- GetAOU(common = common)
    alpha <- code_lookup$alpha[code_lookup$common == common]
  }

  years <- seq(from = start.year, to = end.year)
  nYears <- length(years)
  ten.stops <- sum(grepl("count", names(bbs_raw$counts))) > 0

  spp_counts <- dplyr::filter(bbs_raw$counts, aou == AOU)

  if(!is.null(countrynum)){spp_counts <- dplyr::filter(spp_counts, grepl(paste("^", countrynum, sep = ""), routeID))}
  if(!is.null(statenum)){spp_counts <- dplyr::filter(spp_counts, regexpr(as.character(statenum), routeID) == 4)}


  run_atrb <- dplyr::filter(bbs_raw$weather, routeID %in% spp_counts$routeID & RPID %in% c(101, 501))

  if(!is.null(years)){run_atrb <- dplyr::filter(bbs_raw$weather,  routeID %in% spp_counts$routeID)}
  run_atrb$StartTemp <- as.integer(run_atrb$StartTemp)
  run_atrb$EndTemp <- as.integer(run_atrb$EndTemp)
  run_atrb$StartTemp[run_atrb$TempScale %in% c("f", "F", "N")] <- (run_atrb$StartTemp[run_atrb$TempScale %in%  c("f", "F", "N")] - 32)/1.8
  run_atrb$EndTemp[run_atrb$TempScale %in%  c("f", "F", "N")] <- (run_atrb$EndTemp[run_atrb$TempScale %in%  c("f", "F", "N")] - 32)/1.8

  ## Add novice observer dummy variable
  run_atrb <- dplyr::arrange(run_atrb, routeID, Year)
  run_atrb <- dplyr::mutate(run_atrb, novice = as.integer(!duplicated(ObsN)))

  ## Add distance and time-of-detection dummary variable
  run_atrb <- dplyr::mutate(run_atrb, twedt = ifelse(RPID == 101, 0, 1))

  ## Add year day
  run_atrb <- dplyr::mutate(run_atrb, day = lubridate::yday(lubridate::mdy(paste(run_atrb$Month, run_atrb$Day, run_atrb$Year, sep = "/"))))
  run_atrb <- dplyr::select(run_atrb, routeID, Year, RunType, StartTemp, StartWind, StartTime, novice, twedt, ObsN, day)
  run_atrb <- dplyr::filter(run_atrb, RunType == 1)
  run_atrb <- run_atrb[!duplicated(run_atrb),]

  ### Add RunType to count data
  spp_counts <- suppressMessages(dplyr::left_join(spp_counts, run_atrb))

  ### Remove any routes with no route-runs with runtype == 1
  ### Run info for 0 counts (run has weather data but no count data)
  run_atrb2 <- dplyr::filter(run_atrb, routeID %in% spp_counts$routeID)
  count0 <- suppressMessages(dplyr::anti_join(run_atrb2, spp_counts))

  ### Add 0 counts to data & fill in AOU code
  spp_counts_full <- suppressMessages(dplyr::full_join(count0, spp_counts))
  spp_counts_full$aou <- AOU
  spp_counts_full <- dplyr::select(spp_counts_full, -RunType)

  ### Fill in 0 counts
  spp_counts_full[is.na(spp_counts_full)] <- 0

  ### Add longitude and latitude
  route_atrb <- dplyr::select(bbs_raw$routes, routeID, Latitude, Longitude, Stratum, BCR)
  spp_counts_full <- suppressMessages(dplyr::left_join(spp_counts_full, route_atrb))
  spp_counts_full <- spp_counts_full[!duplicated(spp_counts_full),]
  if(!is.null(strata)){spp_counts_full <- dplyr::filter(spp_counts_full, Stratum %in% strata)}
  if(!is.null(bcr)){spp_counts_full <- dplyr::filter(spp_counts_full, BCR %in% bcr)}
  if(!is.null(lon)){spp_counts_full <- dplyr::filter(spp_counts_full, Longitude < max(lon) & Longitude > min(lon))}
  ### Filter only years of interest
  spp_counts_full <- dplyr::filter(spp_counts_full,  Year %in% years)
  spp_counts_full <- RemoveOutliers(raw.counts = spp_counts_full, Write = FALSE)

  ### Buffer routes w/ observed counts
  ## Create convex hull
  bbs_pts <- sp::SpatialPoints(spp_counts_full[,c('Longitude','Latitude')])
  hull <- rgeos::gConvexHull(bbs_pts)
  buff_hull <- rgeos::gBuffer(hull, width = buffer)

  ## Identify routes within convex hull
  all_coords <- sp::SpatialPoints(route_atrb[,c('Longitude','Latitude', "routeID")])
  inside_hull <- !is.na(sp::over(all_coords, buff_hull))
  buff_routes	<- route_atrb[inside_hull == "TRUE",]
  buff_routes <- suppressMessages(dplyr::anti_join(buff_routes, spp_counts, by = "routeID"))
  buff_routes <- dplyr::select(buff_routes, routeID, Latitude, Longitude, Stratum, BCR)

  ## Create data frame containing years that buffered routes were run (0 count)
  buff_run <- dplyr::filter(bbs_raw$weather,  routeID %in% buff_routes$routeID & Year %in% seq(from = start.year, to = end.year) & RunType == 1 & RPID %in% c(101, 501))
  buff_run <- dplyr::distinct(buff_run, routeID, Year, .keep_all = FALSE)

  ## Add lat, long, stratum, & BCR
  buff_run <- suppressMessages(dplyr::left_join(buff_run, buff_routes))

  ### Create data frame with 0 counts for buffered routes
  col_counts <- grep("count|stop", names(spp_counts_full), value = TRUE)
  count_buff <- dplyr::as_data_frame(matrix(0, nrow = nrow(buff_run),
                                            ncol = length(col_counts)))
  names(count_buff) <- col_counts
  count_buff <- dplyr::bind_cols(buff_run, count_buff)
  count_buff$aou <- unique(spp_counts_full$aou)

  spp_counts_full <- dplyr::bind_rows(spp_counts_full, count_buff)


  ### Covert count data to long format
  counts <- dplyr::select(spp_counts_full, routeID, Year, grep("count|stop", names(spp_counts_full)))
  if(ten.stops) {counts <- dplyr::select(counts, -grep("stoptotal", names(counts)))}
  counts <- tidyr::gather(counts, key = "stop", value = "n", -routeID, -Year)


  ### Add column with presence/absence data
  counts <- dplyr::mutate(counts, occ = ifelse(n > 0, 1, 0))



  ### Covert back to wide w/ 1 column for each year/stop (i.e., svy)
  counts <- tidyr::unite(counts, svy, Year, stop, sep = "_")
  counts <- counts[!duplicated(counts),]
  if(occ){
    counts <- dplyr::select(counts, -n)
    counts <- tidyr::spread(counts, key = svy, value = occ)
  }else{
    counts <- dplyr::select(counts, -occ)
    counts <- tidyr::spread(counts, key = svy, value = n)
  }

  nRoutes <- length(unique(counts$routeID))
  nStops <- 50
  if(ten.stops) nStops <- 5

  counts <- dplyr::arrange(counts, routeID)


  if(nYears > 1){
    h <- NULL
    for(t in 1:nYears){
     ht <- counts[,grep(years[t], names(counts))]
     h <- abind::abind(h, ht, along = 3)
    }
  }else{
    h <- dplyr::select(counts, -routeID)
  }


  covs <- dplyr::select(spp_counts_full, routeID, Year, StartTemp, StartWind, StartTime, novice, twedt, ObsN, day)
  temp <- dplyr::select(covs, routeID, Year, StartTemp)
  temp <- temp[!duplicated(temp[,-3]),]
  temp$StartTemp <- scale(as.numeric(temp$StartTemp))[,1]
  temp <- tidyr::spread(temp, key = Year, value = StartTemp)
  temp <- dplyr::select(temp, -routeID)
  temp[is.na(temp)] <- 0

  wind <- dplyr::select(covs, routeID, Year, StartWind)
  wind <- wind[!duplicated(wind[,-3]),]
  wind$StartWind <- (as.numeric(wind$StartWind) - 1)/(max(as.numeric(wind$StartWind), na.rm = TRUE) - 1)
  wind <- tidyr::spread(wind, key = Year, value = StartWind)
  wind <- dplyr::select(wind, -routeID)
  wind[is.na(wind)] <- mean(unlist(wind), na.rm = TRUE)

  time <- dplyr::select(covs, routeID, Year, StartTime)
  time <- time[!duplicated(time[,-3]),]
  time$StartTime <- scale(as.numeric(time$StartTime))[,1]
  time <- tidyr::spread(time, key = Year, value = StartTime)
  time <- dplyr::select(time, -routeID)
  time[is.na(time)] <- 0

  nov <- dplyr::select(covs, routeID, Year, novice)
  nov <- nov[!duplicated(nov[,-3]),]
  nov <- tidyr::spread(nov, key = Year, value = novice)
  nov <- dplyr::select(nov, -routeID)
  nov[is.na(nov)] <- 0

  twedt <- dplyr::select(covs, routeID, Year, twedt)
  twedt <- twedt[!duplicated(twedt[,-3]),]
  twedt <- tidyr::spread(twedt, key = Year, value = twedt)
  twedt <- dplyr::select(twedt, -routeID)
  twedt[is.na(twedt)] <- 0

  obs <- dplyr::select(covs, routeID, Year, ObsN)
  obs <- obs[!duplicated(obs[,-3]),]
  obs <- tidyr::spread(obs, key = Year, value = ObsN)
  obs <- dplyr::select(obs, -routeID)
  obs <- as.matrix(obs)
  obs[obs == 0] <- NA
  obs <- matrix(as.numeric(as.factor(obs)), nrow = nRoutes, ncol = nYears)
  nObs <- max(obs, na.rm = TRUE)
  obs[is.na(obs)] <- nObs + 1

  coord <- dplyr::filter(spp_counts_full, !duplicated(routeID))
  coord <- dplyr::arrange(coord, routeID)
  slat <- scale(coord$Latitude)[,1]
  slon <- scale(coord$Longitude)[,1]
  bcr <- coord$BCR


  alpha <- code_lookup$alpha[code_lookup$AOU == AOU]
  dat <- list(alpha = alpha, h = h, temp = temp, time = time, wind = wind, nov = nov, twedt = twedt, obs = obs,
              slat = slat, slon = slon, lat = coord$Latitude, lon = coord$Longitude, bcr = bcr, nRoutes = nRoutes, nStops = nStops,
              nYears = nYears, start_year = start.year, end_year = end.year, hull = hull)

  dir.create(paste(path, alpha, sep = "/"), showWarnings = FALSE)
  saveRDS(object = dat, file = paste(path, alpha, "bbs_data.rds", sep = "/"))
}


