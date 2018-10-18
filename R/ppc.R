#' ppc
#'
#' Perform posterior predictive check on fitted model
#' @param alpha alpha code for species of interest
#' @return Bayesian p-value
#' @export

ppc <- function(spp = NULL, alpha = NULL, path){
  if(!is.null(spp)){
    ### Register cores
    cores <- parallel::detectCores()
    if(length(spp) < cores) cores <- length(spp)
    doParallel::registerDoParallel(cores = cores)

    ### Run posterior predictive checks in parallel
    ppc_run <- foreach::foreach(i = 1:length(spp), .combine = c,
                                 .packages = c("dplyr", "BayesCorrOcc")) %dopar%{

          ### Read data & fitted model object
          dat <- readRDS(paste(path, spp[i], "bbs_data.rds", sep = "/"))
          mod <- readRDS(paste(path, spp[i], "jags_fit.rds", sep = "/"))
          bio <- readRDS(paste(path, spp[i], "biovars.rds", sep = "/"))

          ### Get lat/lon splines
          jagam.data <- data.frame(z = rep(1, length(dat$lat)), x = dat$lon, y = dat$lat)
          jagam.mod <- mgcv::jagam(z ~ s(x, y, k = 60, bs = 'ds', m = c(1, 0.5)), data = jagam.data, family = "binomial",
                                   file = system.file("jags", "jagam.jags", package = "BBS.SDM"))

          ### Empty matrices to store observed and simulated X2 estimates
          fit <- matrix(numeric(length = mod$mcmc.info$n.samples*dat$nYears), nrow = dat$nYears)
          fit.new <- matrix(numeric(length = mod$mcmc.info$n.samples*dat$nYears), nrow = dat$nYears)

          ### Store annual psi estimats for each route
          PSI <- array(dim = c(dat$nRoutes, dat$nYears, mod$mcmc.info$n.samples))

          ### Fit in year 1
          ## Observed # of routes w/ each possible detection history
          obs_hist <- BayesCorrOcc::obs_h(dat$h, 1)

          ### For each posterior sample...
          for(ii in 1:mod$mcmc.info$n.samples){
            ### Matrix to store simulated detection histories
            sim_h <- matrix(NA, nrow = dat$nRoutes, ncol = 5)

            ### For each route, estimate psi and p
            PSI[, 1, ii] <- plogis(jagam.mod$jags.data$X %*% mod$sims.list$b[ii,,1] +
                                     bio[,,1] %*% (mod$sims.list$g[ii,] * mod$sims.list$betaT[ii,]))
            p1 <- plogis(mod$sims.list$alpha0[ii] + mod$sims.list$alpha1[ii] * dat$wind[, 1] +
                           mod$sims.list$alpha2[ii] * dat$nov[, 1] + mod$sims.list$alpha3[ii] * dat$twedt[, 1] + mod$sims.list$omega[dat$obs[, 1]])

            ### Estimate expected prob for each possible detection history|psi, p, xpsi
            exp_pr <- BayesCorrOcc::y_probs(psi = PSI[obs_hist$no_hist, 1, ii], xpsi = mod$sims.list$xpsi[ii,], p = p1[obs_hist$no_hist])

            ### Equilibrium stop-level availability
            pi <- mod$sims.list$xpsi[ii, 1] / (mod$sims.list$xpsi[ii, 1] + 1 - mod$sims.list$xpsi[ii, 2])

            ### Simulated availability and detection histories
            sim_y <- matrix(nrow = dat$nRoutes, ncol = 5)
            sim_h <- matrix(nrow = dat$nRoutes, ncol = 5)

            ### Simulate route-level occupancy|psi
            z.new <- rbinom(n = dat$nRoutes, size = 1, prob = PSI[, 1, ii])

            ### Simulate stop-level availability|z.new, xpsi
            sim_y[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * pi)
            sim_h[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, 1] * p1)
            for(k in 2:5){
              sim_y[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * mod$sims.list$xpsi[ii, (sim_y[, k - 1] + 1)])
              sim_h[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, k] * p1)
            }

            ### Counts of each possible detection history of simulated histories
            sim_hist <- BayesCorrOcc::obs_h(sim_h[obs_hist$no_hist,], sim = TRUE)

            ### Expected counts of each possible detection history
            exp_h <- apply(exp_pr, 2, sum)

            ### Observed & simulated X2 statistics
            fit[1, ii] <- sum((obs_hist$obs_h - exp_h)^2/exp_h)
            fit.new[1, ii] <- sum((sim_hist - exp_h)^2/exp_h)
          }


          ### Fit for years 2 - nYears
          for(t in 2:dat$nYears){
            ## Observed # of routes w/ each possible detection history
            obs_hist <- BayesCorrOcc::obs_h(dat$h, t)


            ### For each posterior sample...
            for(ii in 1:mod$mcmc.info$n.samples){
              ### Matrices to store expected and simulated detection histories
              exp_pr <- matrix(NA, nrow = dat$nRoutes, ncol = 32)
              sim_h <- matrix(NA, nrow = dat$nRoutes, ncol = 5)

              ### For each route, estimate psi and p
              PSI[, t, ii] <- plogis(jagam.mod$jags.data$X %*% mod$sims.list$b[i,,t] +
                                       bio[,,t] %*% (mod$sims.list$g[i,] * mod$sims.list$betaT[i,]))
              p1 <- plogis(mod$sims.list$alpha0[ii] + mod$sims.list$alpha1[ii] * dat$wind[, t] +
                             mod$sims.list$alpha2[ii] * dat$nov[, t] + mod$sims.list$alpha3[ii] * dat$twedt[, t] + mod$sims.list$omega[dat$obs[, t]])

              ### Estimate expected prob for each possible detection history|psi, p, xpsi
              exp_pr <- BayesCorrOcc::y_probs(psi = PSI[obs_hist$no_hist, t, ii], xpsi = mod$sims.list$xpsi[ii,], p = p1[obs_hist$no_hist])

              ### Equilibrium stop-level availability
              pi <- mod$sims.list$xpsi[ii, 1] / (mod$sims.list$xpsi[ii, 1] + 1 - mod$sims.list$xpsi[ii, 2])

              ### Simulated availability and detection histories
              sim_y <- matrix(nrow = dat$nRoutes, ncol = 5)
              sim_h <- matrix(nrow = dat$nRoutes, ncol = 5)

              ### Simulate route-level occupancy|psi
              z.new <- rbinom(n = dat$nRoutes, size = 1, prob = PSI[, t, ii])

              ### Simulate stop-level availability|z.new, xpsi
              sim_y[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * pi)
              sim_h[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, 1] * p1)
              for(k in 2:5){
                sim_y[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * mod$sims.list$xpsi[ii, (sim_y[, k - 1] + 1)])
                sim_h[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, k] * p1)
              }


              ### Counts of each possible detection history of simulated histories
              sim_hist <- BayesCorrOcc::obs_h(sim_h[obs_hist$no_hist,], sim = TRUE)

              ### Expected counts of each possible detection history
              exp_h <- apply(exp_pr, 2, sum)

              ### Observed & simulated X2 statistics
              fit[t, ii] <- sum((obs_hist$obs_h - exp_h)^2/exp_h)
              fit.new[t, ii] <- sum((sim_hist - exp_h)^2/exp_h)
        }

    }

    fit_df <- data.frame(fit = c(fit), fit.new = c(fit.new),
                         Year = rep(seq(1:dat$nYears), mod$mcmc.info$n.samples))
    ppc <- list(fit = fit_df, p = sum(fit > fit.new)/length(fit))

    saveRDS(ppc, file = paste(path, spp[i], "ppc.rds", sep = "/"))
    return(spp[i])
    }
   return(ppc_run)
  }

  if(!is.null(alpha)){
    ### Read data & fitted model object
    dat <- readRDS(paste(path, alpha, "bbs_data.rds", sep = "/"))
    mod <- readRDS(paste(path, alpha, "jags_fit.rds", sep = "/"))
    bio <- readRDS(paste(path, alpha, "biovars.rds", sep = "/"))

    ### Get lat/lon splines
    jagam.data <- data.frame(z = rep(1, length(dat$lat)), x = dat$lon, y = dat$lat)
    jagam.mod <- mgcv::jagam(z ~ s(x, y, k = 60, bs = 'ds', m = c(1, 0.5)),
                             data = jagam.data, family = "binomial",
                             file = system.file("jags", "jagam.jags", package = "BBS.SDM"))

    ### Empty matrices to store observed and simulated X2 estimates
    fit <- matrix(numeric(length = mod$mcmc.info$n.samples*dat$nYears), nrow = dat$nYears)
    fit.new <- matrix(numeric(length = mod$mcmc.info$n.samples*dat$nYears), nrow = dat$nYears)

    ### Store annual psi estimats for each route
    PSI <- array(dim = c(dat$nRoutes, dat$nYears, mod$mcmc.info$n.samples))

    ### Fit in year 1
    ## Observed # of routes w/ each possible detection history
    obs_hist <- BayesCorrOcc::obs_h(dat$h, 1)

    ### For each posterior sample...
    for(i in 1:mod$mcmc.info$n.samples){
      ### Matrix to store simulated detection histories
      sim_h <- matrix(NA, nrow = dat$nRoutes, ncol = 5)

      ### For each route, estimate psi and p
      PSI[, 1, i] <- plogis(jagam.mod$jags.data$X %*% mod$sims.list$b[i,,1] +
                              bio[,,1] %*% (mod$sims.list$g[i,] * mod$sims.list$betaT[i,]))
      p1 <- plogis(mod$sims.list$alpha0[i] + mod$sims.list$alpha1[i] * dat$wind[, 1] +
                     mod$sims.list$alpha2[i] * dat$nov[, 1] + mod$sims.list$alpha3[i] * dat$twedt[, 1] + mod$sims.list$omega[dat$obs[, 1]])

      ### Estimate expected prob for each possible detection history|psi, p, xpsi
      exp_pr <- BayesCorrOcc::y_probs(psi = PSI[obs_hist$no_hist, 1, i], xpsi = mod$sims.list$xpsi[i,1,], p = p1[obs_hist$no_hist])

      ### Equilibrium stop-level availability
      pi <- mod$sims.list$xpsi[i, 1, 1] / (mod$sims.list$xpsi[i, 1, 1] + (1 - mod$sims.list$xpsi[i, 1, 2]))

      ### Simulated availability and detection histories
      sim_y <- matrix(nrow = dat$nRoutes, ncol = 5)
      sim_h <- matrix(nrow = dat$nRoutes, ncol = 5)

      ### Simulate route-level occupancy|psi
      z.new <- rbinom(n = dat$nRoutes, size = 1, prob = PSI[, 1, i])

      ### Simulate stop-level availability|z.new, xpsi
      sim_y[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * pi)
      sim_h[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, 1] * p1)
      for(k in 2:5){
        sim_y[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * mod$sims.list$xpsi[i, 1, (sim_y[, k - 1] + 1)])
        sim_h[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, k] * p1)
      }

      ### Counts of each possible detection history of simulated histories
      sim_hist <- BayesCorrOcc::obs_h(sim_h[obs_hist$no_hist,], sim = TRUE)

      ### Expected counts of each possible detection history
      exp_h <- apply(exp_pr, 2, sum)

      ### Observed & simulated X2 statistics
      fit[1, i] <- sum((obs_hist$obs_h^(1/2) - exp_h^(1/2))^2)
      fit.new[1, i] <- sum((sim_hist^(1/2) - exp_h^(1/2))^2)
    }


    ### Fit for years 2 - nYears
    for(t in 2:dat$nYears){
      ## Observed # of routes w/ each possible detection history
      obs_hist <- BayesCorrOcc::obs_h(dat$h, t)


      ### For each posterior sample...
      for(i in 1:mod$mcmc.info$n.samples){
        ### Matrices to store expected and simulated detection histories
        exp_pr <- matrix(NA, nrow = dat$nRoutes, ncol = 32)
        sim_h <- matrix(NA, nrow = dat$nRoutes, ncol = 5)

        ### For each route, estimate psi and p
        PSI[, t, i] <- plogis(jagam.mod$jags.data$X %*% mod$sims.list$b[i,,t] +
                                bio[,,t] %*% (mod$sims.list$g[i,] * mod$sims.list$betaT[i,]))
        p1 <- plogis(mod$sims.list$alpha0[i] + mod$sims.list$alpha1[i] * dat$wind[, t] +
                       mod$sims.list$alpha2[i] * dat$nov[, t] + mod$sims.list$alpha3[i] * dat$twedt[, t] + mod$sims.list$omega[dat$obs[, t]])

        ### Estimate expected prob for each possible detection history|psi, p, xpsi
        exp_pr <- BayesCorrOcc::y_probs(psi = PSI[obs_hist$no_hist, t, i], xpsi = mod$sims.list$xpsi[i,t,], p = p1[obs_hist$no_hist])

        ### Equilibrium stop-level availability
        pi <- mod$sims.list$xpsi[i, t, 1] / (mod$sims.list$xpsi[i, t, 1] + 1 - mod$sims.list$xpsi[i, t, 2])

        ### Simulated availability and detection histories
        sim_y <- matrix(nrow = dat$nRoutes, ncol = 5)
        sim_h <- matrix(nrow = dat$nRoutes, ncol = 5)

        ### Simulate route-level occupancy|psi
        z.new <- rbinom(n = dat$nRoutes, size = 1, prob = PSI[, t, i])

        ### Simulate stop-level availability|z.new, xpsi
        sim_y[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * pi)
        sim_h[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, 1] * p1)
        for(k in 2:5){
          sim_y[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * mod$sims.list$xpsi[i, t, (sim_y[, k - 1] + 1)])
          sim_h[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, k] * p1)
        }


        ### Counts of each possible detection history of simulated histories
        sim_hist <- BayesCorrOcc::obs_h(sim_h[obs_hist$no_hist,], sim = TRUE)

        ### Expected counts of each possible detection history
        exp_h <- apply(exp_pr, 2, sum)

        ### Observed & simulated X2 statistics
        fit[t, i] <- sum((obs_hist$obs_h^(1/2) - exp_h^(1/2))^2)
        fit.new[t, i] <- sum((sim_hist^(1/2) - exp_h^(1/2))^2)
      }

    }

    fit_df <- data.frame(fit = c(fit), fit.new = c(fit.new),
                         Year = rep(seq(1:dat$nYears), mod$mcmc.info$n.samples))
    ppc <- list(fit = fit_df, p = mean(fit > fit.new))

    saveRDS(ppc, file = paste(path, dat$alpha, "ppc.rds", sep = "/"))
  }

}


#' obs_h
#'
#' Count observed number of routes with each detection history
#' @param h detection history matrix
#' @param year Year for with observed counts needed
#' @return List containing observed counts for each detection history & index of NA counts
#' @export


obs_h <- function(h, year, sim = FALSE){
  hists <- c("00000", "00001", "00010", "00011", "00100", "00101", "00110", "00111",
             "01000", "01001", "01010", "01011", "01100", "01101", "01110", "01111",
             "10000", "10001", "10010", "10011", "10100", "10101", "10110", "10111",
             "11000", "11001", "11010", "11011", "11100", "11101", "11110", "11111")

  if(sim){
    x <- as.data.frame(h)
    x2 <- tidyr::unite(x, hist, 1:5, sep = "")
    y <- table(x2)

    non_zero <- hists[which(hists %in% names(y))]

    obs <- integer(length = 32)
    names(obs) <- hists

    for(i in 1:length(non_zero)){
      obs[names(obs) == non_zero[i]] <- y[names(y) == non_zero[i]]
    }
    return(obs)
  }else{
    x <- as.data.frame(h[,,year])
    x2 <- tidyr::unite(x, hist, 1:5, sep = "")
    y <- table(x2)

    non_zero <- hists[which(hists %in% names(y))]

    obs <- integer(length = 32)
    names(obs) <- hists

    for(i in 1:length(non_zero)){
      obs[names(obs) == non_zero[i]] <- y[names(y) == non_zero[i]]
    }

    no_hist <- which(x2 != "NANANANANA")

    obs_H <- list(obs_h = obs, no_hist = no_hist)
    return(obs_H)
  }

}
