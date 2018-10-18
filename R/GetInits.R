#' GetInits
#'
#' Estimate model coefficients in Presence; save output for mu/se values & inits
#' @param alpha four letter alpha code for species of interest; if NULL, runs Presence model for all species in parallel
#' @export

GetInits <- function(alpha = NULL, path){
    if(is.null(alpha)){
      ### Read species list
      spp_list <- read.csv(here::here("inst/spp_list.csv"))

      ### Register cores
      cores <- parallel::detectCores()
      if(length(spp_list$spp) < cores) cores <- length(spp_list$spp)
      doParallel::registerDoParallel(cores = cores)

      ### Run models in parallel
      inits <- foreach::foreach(i = 1:length(spp_list$spp), .combine = c,
                       .packages = c("dplyr", "BBS.SDM")) %dopar%{
                         ## Read count data
                         dat <-  readRDS(paste(path, spp_list$spp[i], "bbs_data.rds", sep = "/"))

                         ## Write pao file
                         BBS.SDM::write_pao(alpha = spp_list$spp[i])

                         ## Read pao file
                         pao <- RPresence::readPao(paste(path, spp_list$spp[i], "pres_in.pao", sep = "/"))

                         ## Create design matrices
                         start_yr <- dat$start_year

                         n_surv <- pao$nsurveys  # total num of surveys
                         n_seas <- pao$nseasons	# num seasons
                         n_count <- pao$nsurveyseason[1]

                         years <- seq(from = start_yr, to = start_yr + n_seas - 1)

                         num.betas	<- c(11,     # psi
                                        1,      #th0
                                        1,      #th1
                                        11,     #gamma
                                        11,     #epsilon
                                        5)      #p

                         # 1st design matrix, 1 row for psi, K rows for theta0, K rows for theta1
                         dm1 <- matrix('0', n_surv * 2 + 1, sum(num.betas[1:3]))

                         # psi
                         dm1[1, 1:num.betas[1]] <- c(1, colnames(pao$unitcov)[grepl(start_yr, colnames(pao$unitcov))])

                         # th0
                         dm1[2:(n_surv + 1),(1 + num.betas[1]):sum(num.betas[1:2])] <- rep(1, n_surv)

                         # th1
                         dm1[(n_surv + 2):(2 * n_surv + 1), (1 + sum(num.betas[1:2])):sum(num.betas[1:3])] <- rep(1, n_surv)

                         rownames(dm1) <- c('psi', paste0('th0(', 1:n_surv,')'), paste0('th1(', 1:n_surv,')'))
                         colnames(dm1) <- paste0('a',1:sum(num.betas[1:3]))

                         # gam
                         gam.dm <- c(1, colnames(pao$unitcov)[grepl((start_yr + 1), colnames(pao$unitcov))])
                         for (ii in 2:(length(years) - 1)) {
                           gam.dm <- c(gam.dm, c(1, colnames(pao$unitcov)[grepl((start_yr + ii), colnames(pao$unitcov))]))
                         }

                         dm2 <- matrix(gam.dm, n_seas - 1, num.betas[4], byrow = T,
                                       dimnames = list(paste0('gam', 1:(n_seas-1)), paste0('b', 1:num.betas[4])))


                         # eps
                         dm3 <- matrix(gam.dm, n_seas - 1, num.betas[5], byrow = T,
                                       dimnames = list(paste0('eps', 1:(n_seas-1)), paste0('c', 1:num.betas[4])))

                         # p
                         scale.stop <- c("Stop1", "sq_Stop1")
                         for (ii in 2:n_count) {
                           scale.stop <- c(scale.stop, paste0("Stop", ii), paste0("sq_Stop", ii))
                         }
                         scale.stop <- matrix(scale.stop, nrow = n_surv, ncol = 2, byrow = T)

                         coord.mat <- matrix(rep(c("Lat", "Lon"), n_surv), nrow = n_surv, byrow = TRUE)

                         p1.intercept <- rep(1, n_surv)
                         dm4 <- cbind(p1.intercept, scale.stop, coord.mat)

                         rownames(dm4) <- c(paste0('p1(',1:n_surv,')'))     	# could get rid of parentheses
                         colnames(dm4) <- paste0('d', 1:dim(dm4)[2])

                         # theta.pi
                         dm5 <- matrix(0, n_seas, 1, dimnames = list(paste0('thta0pi', 1:n_seas), NULL)) 		# NOTE THE ZERO when fixing, note no colname


                         dm_list <- list(dm1 = dm1, dm2 = dm2, dm3 = dm3, dm4 = dm4, dm5 = dm5, dm6 = NULL)


                         fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
                         r1 <- dim(dm_list$dm1)[1] + dim(dm_list$dm2)[1] + dim(dm_list$dm3)[1] + dim(dm_list$dm4)[1]
                         rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

                         ## Run model
                         RPresence::runPresence(a = pao, fixedvals = fixedpars,
                                                dm = dm_list, model = "Inits",
                                                noDerived = TRUE, limit.real = TRUE,
                                                modname = paste0("inits_", spp_list$spp[i]))

                         suppressMessages(file.rename(from = here::here(paste0("pres_inits_", tolower(spp_list$spp[i]), ".out")),
                                                      to = paste(path, toupper(spp_list$spp[i]), "inits.out", sep = "/")))

                         inits <- suppressWarnings(GetBetas(alpha = spp_list$spp[i], path = path))
                         saveRDS(inits, paste(path, spp_list$spp[i], "inits.rds", sep = "/"))
                         return(as.character(spp_list$spp[i]))
                       }
      return(inits)
    }else{
      ## Read count data
      dat <-  readRDS(paste(path, alpha, "bbs_data.rds", sep = "/"))

      ## Write pao file
      BBS.SDM::write_pao(alpha = alpha)

      ## Read pao file
      pao <- RPresence::readPao(paste(path, alpha, "pres_in.pao", sep = "/"))

      ## Create design matrices
      start_yr <- dat$start_year

      n_surv <- pao$nsurveys  # total num of surveys
      n_seas <- pao$nseasons	# num seasons
      n_count <- pao$nsurveyseason[1]

      years <- seq(from = start_yr, to = start_yr + n_seas - 1)

      num.betas	<- c(11,     # psi
                     1,      #th0
                     1,      #th1
                     11,     #gamma
                     11,     #epsilon
                     5)      #p

      # 1st design matrix, 1 row for psi, K rows for theta0, K rows for theta1
      dm1 <- matrix('0', n_surv * 2 + 1, sum(num.betas[1:3]))

      # psi
      dm1[1, 1:num.betas[1]] <- c(1, colnames(pao$unitcov)[grepl(start_yr, colnames(pao$unitcov))])

      # th0
      dm1[2:(n_surv + 1),(1 + num.betas[1]):sum(num.betas[1:2])] <- rep(1, n_surv)

      # th1
      dm1[(n_surv + 2):(2 * n_surv + 1), (1 + sum(num.betas[1:2])):sum(num.betas[1:3])] <- rep(1, n_surv)

      rownames(dm1) <- c('psi', paste0('th0(', 1:n_surv,')'), paste0('th1(', 1:n_surv,')'))
      colnames(dm1) <- paste0('a',1:sum(num.betas[1:3]))

      # gam
      gam.dm <- c(1, colnames(pao$unitcov)[grepl((start_yr + 1), colnames(pao$unitcov))])
      for (ii in 2:(length(years) - 1)) {
        gam.dm <- c(gam.dm, c(1, colnames(pao$unitcov)[grepl((start_yr + ii), colnames(pao$unitcov))]))
      }

      dm2 <- matrix(gam.dm, n_seas - 1, num.betas[4], byrow = T,
                    dimnames = list(paste0('gam', 1:(n_seas-1)), paste0('b', 1:num.betas[4])))


      # eps
      dm3 <- matrix(gam.dm, n_seas - 1, num.betas[5], byrow = T,
                    dimnames = list(paste0('eps', 1:(n_seas-1)), paste0('c', 1:num.betas[4])))

      # p
      scale.stop <- c("Stop1", "sq_Stop1")
      for (ii in 2:n_count) {
        scale.stop <- c(scale.stop, paste0("Stop", ii), paste0("sq_Stop", ii))
      }
      scale.stop <- matrix(scale.stop, nrow = n_surv, ncol = 2, byrow = T)

      coord.mat <- matrix(rep(c("Lat", "Lon"), n_surv), nrow = n_surv, byrow = TRUE)

      p1.intercept <- rep(1, n_surv)
      dm4 <- cbind(p1.intercept, scale.stop, coord.mat)
      # dm4 <- rbind(dm4, dm4)
      rownames(dm4) <- c(paste0('p(',1:n_surv,')'))#,paste0('p2(',1:n_surv,')'))     	# could get rid of parentheses
      colnames(dm4) <- paste0('d', 1:dim(dm4)[2])

      # theta.pi
      dm5 <- matrix(0, n_seas, 1, dimnames = list(paste0('thta0pi', 1:n_seas), NULL)) 		# NOTE THE ZERO when fixing, note no colname


      dm_list <- list(dm1 = dm1, dm2 = dm2, dm3 = dm3, dm4 = dm4, dm5 = dm5, dm6 = NULL)


      fixedpars <- data.frame(params = rep("eq", pao$nseasons))
      r1 <- dim(dm_list$dm1)[1] + dim(dm_list$dm2)[1] + dim(dm_list$dm3)[1] + dim(dm_list$dm4)[1]
      rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

      ## Run model
        RPresence::runPresence(a = pao, fixedvals = fixedpars,
                               dm = dm_list,
                               noDerived = TRUE, limit.real = TRUE,
                               modname = "inits", outfile = here::here(paste0("inst/output", alpha, "/inits.out")))

      suppressMessages(file.rename(from = here::here("pres_inits.out"),
                                   to = paste(path, alpha, "inits.out", sep = "/")))

      inits <- suppressWarnings(GetBetas(alpha = alpha, path = path))
      saveRDS(inits, paste(path, alpha, "inits.rds", sep = "/"))
    }
}

