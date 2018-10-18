#' RunMod
#'
#' Run model
#' @param alpha alpha code for species of interest; if NULL, runs model in parallel for all species
#' @export

RunMod <- function(spp = NULL, alpha = NULL, path, nI = 50000, nB = 10000, nC = 2, nT = 50, cores = 50, Parallel = FALSE){
  if(!is.null(spp)){
    ## Check if model has already been run for species
    spp_run <- NULL
    for(i in 1:length(spp)){
      spp_test <- file.exists(paste(path, spp[i], "/jags_fit.rds", sep = "/"))
      if(!spp_test) spp_run <- c(spp_run, as.character(spp[i]))
    }

    if(is.null(spp_run)){
      return(spp)
    }else{
      ### Register core
      if(is.null(cores)){
        cores <- parallel::detectCores()
      }
      if(length(spp_run) < cores) cores <- length(spp_run)
      doParallel::registerDoParallel(cores = cores)

      ### Run models in parallel
      mods_run <- foreach::foreach(i = 1:length(spp_run), .combine = c,
                                   .packages = c("dplyr", "BBS.SDM", "jagsUI")) %dopar%{

                                     ### Read data
                                     dat <-  readRDS(paste(path, spp_run[i], "bbs_data.rds", sep = "/"))
                                     covs <- readRDS(paste(path, spp_run[i], "biovars.rds", sep = "/"))
                                     inits <- readRDS(paste(path, spp_run[i], "inits.rds", sep = "/"))

                                     ### For inits w/ very large (i.e., likely poorly estimated), change to 0 to ensure model doesn't choke
                                     inits$psi.se[abs(inits$psi.betas) > 8] <- 1
                                     inits$psi.betas[abs(inits$psi.betas) > 8] <- 0
                                     inits$p.betas[abs(inits$p.betas) > 8] <- 0
                                     inits$th.betas[abs(inits$th.betas) > 8] <- 0

                                     ### Get data for GAM JAGS model
                                     jagam.data <- data.frame(z = rep(1, length(dat$lat)), x = dat$lon, y = dat$lat)
                                     jagam.mod <- mgcv::jagam(z ~ s(x, y, k = 60, bs = 'ds', m = c(1, 0.5)),
                                                              data = jagam.data, family = "binomial",
                                                              file = system.file("jags", "jagam.jags", package = "BBS.SDM"))


                                     ### Data for JAGS
                                     jags.data <- list(h = dat$h, nStops = dat$nStops, nRoutes = dat$nRoutes, nYears = dat$nYears,
                                                       Xp = dat$wind, nov = dat$nov, obs = dat$obs, nObs = max(dat$obs) - 1,
                                                       Xclim = covs, nPred = dim(covs)[2]/2, twedt = dat$twedt,
                                                       X = jagam.mod$jags.data$X, S1 = jagam.mod$jags.data$S1, zero = jagam.mod$jags.data$zero,
                                                       mu = inits$psi.betas, se = inits$psi.se)


                                     ### Parameters to monitor
                                     jags.params <- c("xpsi", "pi", "lambda", "betaT", "g",
                                                      "alpha0", "alpha1", "alpha2", "alpha3", "sigma.obs",
                                                      "sigma.gam", "rho", "b", "omega", "z")


                                     ### Initial values
                                     y <- dat$h
                                     y[is.na(y)] <- rbinom(n = length(y[is.na(y)]), size = 1, prob = 0.5)
                                     jags.inits <- function(){list(y = y, z = apply(dat$h, c(1, 3), max),
                                                                   alpha0 = inits$p.betas[1], alpha1 = rnorm(1),
                                                                   alpha2 = rnorm(1), alpha3 = rnorm(1),
                                                                   omega = c(rnorm(max(dat$obs) - 1), NA), sigma.obs = runif(1, 0, 5),
                                                                   betaT = inits$psi.betas,
                                                                   g = rbinom(dim(covs)[2], size = 1, prob = 0.5),
                                                                   #xpsi = plogis(inits$th.betas),
                                                                   b = matrix(rep(jagam.mod$jags.ini$b, each = 43), ncol = 43, byrow = TRUE),
                                                                   lambda = jagam.mod$jags.ini$lambda)}


                                     ### Fit model
                                     mod <- system.file("jags", "cor_Occ_dyn.jags", package = "BBS.SDM")
                                     jags.fit <- jagsUI::jags(data = jags.data, parameters.to.save = jags.params,
                                                              inits = jags.inits, model.file = mod,
                                                              n.chains = nC, n.iter = nI, n.adapt = nB/2, n.burnin = nB, n.thin = nT,
                                                              parallel = Parallel, verbose = FALSE)

                                     ### Save output
                                     saveRDS(jags.fit, paste(path, spp_run[i], "jags_fit.rds", sep = "/"))
                                     return(spp_run[i])
                                   }
      return(mods_run)
    }
  }

  if(!is.null(alpha)){
    ### Read data
    dat <-  readRDS(paste(path, alpha, "bbs_data.rds", sep = "/"))
    covs <- readRDS(paste(path, alpha, "biovars.rds", sep = "/"))
    inits <- readRDS(paste(path, alpha, "inits.rds", sep = "/"))

    ### For inits w/ very large (i.e., likely poorly estimated), change to 0 to ensure model doesn't choke
    inits$psi.se[abs(inits$psi.betas) > 8] <- 1
    inits$psi.betas[abs(inits$psi.betas) > 8] <- 0
    inits$p.betas[abs(inits$p.betas) > 8] <- 0
    inits$th.betas[abs(inits$th.betas) > 8] <- 0

    ### Get data for GAM JAGS model
    jagam.data <- data.frame(z = rep(1, length(dat$lat)), x = dat$lon, y = dat$lat)
    jagam.mod <- mgcv::jagam(z ~ s(x, y, k = 60, bs = 'ds', m = c(1, 0.5)),
                             data = jagam.data, family = "binomial",
                             file = system.file("jags", "jagam.jags", package = "BBS.SDM"))

    ### Data for JAGS
    jags.data <- list(h = dat$h, nStops = dat$nStops, nRoutes = dat$nRoutes, nYears = dat$nYears,
                      Xp = dat$wind, nov = dat$nov, obs = dat$obs, nObs = max(dat$obs) - 1,
                      Xclim = covs, nPred = dim(covs)[2]/2, twedt = dat$twedt,
                      X = jagam.mod$jags.data$X, S1 = jagam.mod$jags.data$S1, zero = jagam.mod$jags.data$zero,
                      mu = inits$psi.betas, se = inits$psi.se)


    ### Parameters to monitor
    jags.params <- c("xpsi", "pi", "betaT", "g",
                     "alpha0", "alpha1", "alpha2", "alpha3", "lambda", "sigma.obs",
                     "sigma.gam", "rho", "b", "omega", "z", "psi")


    ### Initial values
    y <- dat$h
    y[is.na(y)] <- rbinom(n = length(y[is.na(y)]), size = 1, prob = 0.5)
    jags.inits <- function(){list(y = y, z = apply(dat$h, c(1, 3), max),
                                  alpha0 = inits$p.betas[1], alpha1 = rnorm(1),
                                  alpha2 = rnorm(1), alpha3 = rnorm(1),
                                  omega = c(rnorm(max(dat$obs) - 1), NA), sigma.obs = runif(1, 0, 5),
                                  betaT = inits$psi.betas,
                                  g = rbinom(dim(covs)[2], size = 1, prob = 0.5),
                                  #xpsi = plogis(inits$th.betas),
                                  b = matrix(rep(jagam.mod$jags.ini$b, each = 43), ncol = 43, byrow = TRUE),
                                  lambda = jagam.mod$jags.ini$lambda)}

    ### Fit model
    mod <- system.file("jags", "cor_Occ_dyn.jags", package = "BBS.SDM")
    jags.fit <- jagsUI::jags(data = jags.data, parameters.to.save = jags.params,
                             inits = jags.inits, model.file = mod,
                             n.chains = nC, n.iter = nI, n.adapt = nB/2, n.burnin = nB, n.thin = nT,
                             parallel = Parallel, verbose = FALSE)

    ### Save output
    saveRDS(jags.fit, paste(path, alpha, "jags_fit.rds", sep = "/"))
  }
}

