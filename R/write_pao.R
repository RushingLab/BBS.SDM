#' write_pao
#'
#' Write occupancy and climate data for program Presence
#' @param alpha Four letter alpha code for the species of interest
#' @return A .pao file containing the detection histories, covariates, and summary information to input into Presence
#' @export

write_pao <- function(alpha){

  dat <-  readRDS(here::here(paste0("inst/output/", alpha, "/bbs_data.rds")))
  covs <- readRDS(here::here(paste0("inst/output/", alpha, "/biovars.rds")))

  tot_stops <- dim(dat$h)[2]
  n_seasons <- 10


  ### Create nRoutes x nYears*nStops detection history & climate covariates
  det_hist <- NULL
  sitecovs <- NULL
  for(i in 1:n_seasons){
    det_hist <- cbind(det_hist, dat$h[,,i])
    sitecovs <- cbind(sitecovs, covs[,,i])
  }

  ### Add scaled lat/lon
  sitecovs <- cbind(sitecovs, dat$slat, dat$slon)

  ### Name covariates
  vars <- c("tmp", "dtr", "Twet", "Prec", "Pwarm",
            "sq_tmp", "sq_dtr", "sq_Twet", "sq_Prec", "sq_Pwarm")
  years <- seq(from = dat$start_year, to = dat$start_year + 9)
  colnames(sitecovs) <- c(apply(expand.grid(vars, years), 1, paste, collapse="_"), "Lat", "Lon")


  ## Add stop number
  for(ss in 1:tot_stops) {
    sc <- scale(1:tot_stops)[ss];
    sc2 <- (scale(1:tot_stops)[ss]) ^ 2;
    sitecovs <- cbind(sitecovs, sc, sc2)
    colnames(sitecovs)[ncol(sitecovs) - 1] <- paste0('Stop', ss)
    colnames(sitecovs)[ncol(sitecovs)] <- paste0('sq_Stop', ss)
  }



  pname <- paste0("inst/output/", alpha, "/pres_in.pao")

  nss <- rep(tot_stops, n_seasons)

  spp_pao <- RPresence::createPao(data = det_hist, nsurveyseason = nss,
                                    unitcov = sitecovs, survcov = NULL,
                                    title = paste(alpha, "PRESENCE Analysis", sep = " "),
                                    paoname = pname)


  RPresence::writePao(pao = spp_pao)
}

