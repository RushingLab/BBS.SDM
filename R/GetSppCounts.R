#' GetSppCounts
#'
#' Filter raw BBS counts by species, add 0 counts, and latitude/longitude
#' @param bbs List containing the raw BBS counts, weather, and route info obtained from the function `get_BBSn()`
#' @param AOU Numeric AOU code for the species of interest
#' @param years Optional vector containing the years of interest (default in full study period)
#' @param countrynum Optional numeric vector containing the countries of interest (840 = US, 124 = Canada, 484 = Mexico)
#' @param statenum Optional numeric vector containing the states of interest (see BBS website for code values)
#' @param strata Optional numeric vector containing the BBS strata of interest (see BBS website for code values)
#' @param bcr Optional numeric vector containing the BCRs of interest (see BBS website for code values)
#' @param Write Should data be written as .csv file?
#' @param path Path where .csv file should be saved (default is working directory)
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

GetSppCounts <- function(bbs_raw = bbs, AOU,
                         years = seq(from = 1997, to = 2016), 
                         statenum = NULL, countrynum = NULL, strata = NULL, bcr = NULL,
                         Write = TRUE, path = NULL){

  spp_counts <- dplyr::filter(bbs_raw$counts, aou == AOU & Year %in% years)

  if(!is.null(countrynum)){spp_counts <- dplyr::filter(spp_counts, grepl(paste("^", countrynum, sep = ""), routeID))}
  if(!is.null(statenum)){spp_counts <- dplyr::filter(spp_counts, regexpr(as.character(statenum), routeID) == 4)}
  

  run_atrb <- dplyr::filter(bbs_raw$weather, routeID %in% spp_counts$routeID)
  if(!is.null(years)){run_atrb <- dplyr::filter(bbs_raw$weather,  routeID %in% spp_counts$routeID & Year %in% years)}
  run_atrb <- dplyr::select(run_atrb, routeID, Year, RunType)
  run_atrb <- dplyr::filter(run_atrb, RunType == 1)
  run_atrb <- run_atrb[!duplicated(run_atrb),]
  
  ### Add RunType to count data
  spp_counts <- dplyr::left_join(spp_counts, run_atrb)
  
  ### Remove any routes with no route-runs with runtype == 1
  ### Run info for 0 counts (run has weather data but no count data)
  run_atrb2 <- dplyr::filter(run_atrb, routeID %in% spp_counts$routeID)
  count0 <- dplyr::anti_join(run_atrb2, spp_counts)

  ### Add 0 counts to data & fill in AOU code
  spp_counts_full <- dplyr::full_join(count0, spp_counts)
  spp_counts_full$aou <- AOU
  spp_counts_full <- dplyr::select(spp_counts_full, -RunType)

  ### Fill in 0 counts
  spp_counts_full[is.na(spp_counts_full)] <- 0

  ### Add longitude and latitude
  route_atrb <- dplyr::select(bbs_raw$routes, routeID, Latitude, Longitude, Stratum, BCR)
  spp_counts_full <- dplyr::left_join(spp_counts_full, route_atrb)
  spp_counts_full <- spp_counts_full[!duplicated(spp_counts_full),]
  if(!is.null(strata)){spp_counts_full <- dplyr::filter(spp_counts_full, Stratum %in% strata)}
  if(!is.null(bcr)){spp_counts_full <- dplyr::filter(spp_counts_full, BCR %in% bcr)}
  
  if(Write){
    if(is.null(path)){
      write.csv(spp_counts_full,
                "raw_counts.csv",
                row.names = FALSE)
    }else{
      write.csv(spp_counts_full,
                paste(path, "raw_counts.csv", sep = "/"),
                row.names = FALSE)
    }
  }else{
    spp_counts_full
  }
}


