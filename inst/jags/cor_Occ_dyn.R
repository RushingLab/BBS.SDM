sink(file="inst/jags/cor_Occ_dyn.jags")
cat("
    model {

    #### Prior distributions

    ## Priors for linear indicator variables -- prob = 1 if quadratic term in model, 0.5 otherwise
    for(ii in 1:nPred){
      g[ii] ~ dbern(0.5)
    }


    ## Priors for quadratic indicator variables - prob = 0.5
    for(ii in (nPred + 1):(nPred * 2)){
      g[ii] ~ dbern(g[ii - nPred] * 0.5)
    }


    ## Priors for betas - Normal(0, tau.beta) if g == 1; Normal(mu, 1/se^2) if g == 0
    for(ii in 1:(nPred * 2)){
      ## Initial occupancy
      betaT[ii] ~ dnorm(bpriorm[ii], tprior[ii])T(-10, 10)
      beta[ii] <- g[ii] * betaT[ii]

      bpriorm[ii] <- (1 - g[ii]) * mu[ii]
      tprior[ii] <- g[ii] * 0.01 + (1 - g[ii]) * pow(se[ii], -2)
    }



    ## Detection priors
    alpha0 ~ dnorm(0, 0.1)T(-10, 10)
    alpha1 ~ dnorm(0, 0.1)T(-10, 10)
    alpha2 ~ dnorm(0, 0.1)T(-10, 10)
    alpha3 ~ dnorm(0, 0.1)T(-10, 10)

    for(ii in 1:nObs){
      omega[ii] ~ dnorm(0, tau.obs)
    }
    omega[(nObs + 1)] <- 0

    tau.obs <- pow(sigma.obs, -2)
    sigma.obs ~ dunif(0, 10)


    ## Spatial correlation priors
    for(tt in 1:nYears){
      lxpsi[tt, 1] ~ dnorm(mu.xpsi[1], tau.xpsi)
      logit(xpsi[tt, 1]) <- lxpsi[tt, 1]
      lxpsi[tt, 2] ~ dnorm(mu.xpsi[2], tau.xpsi)
      logit(xpsi[tt, 2]) <- lxpsi[tt, 2]
      pi[tt] ~dunif(0, 1)#<- xpsi[tt, 1]/(xpsi[tt, 1] + (1 - xpsi[tt, 2]))
    }

    mean.xpsi[1] ~ dunif(0, 1)
    mu.xpsi[1] <- log(mean.xpsi[1]) - log(1 - mean.xpsi[1])

    mean.xpsi[2] ~ dunif(0, 1)
    mu.xpsi[2] <- log(mean.xpsi[2]) - log(1 - mean.xpsi[2])

    tau.xpsi <- pow(sigma.xpsi, -2)
    sigma.xpsi ~ dunif(0, 10)
    # pi ~ dunif(0, 1)
    # xpsi[1] ~ dunif(0,1)
    # xpsi[2] ~ dunif(0, 1)

    #### GAM priors from mgcv::jagam()

    ### YEAR 1
    ## Parametric effect priors
    for (ii in 1:1) {
      b[ii, 1] ~ dnorm(0,0.033)
    }

    ## prior for s(sim_dat$xy$x,sim_dat$xy$y)...
    K1 <- S1[1:59,1:59] * lambda[1]
    b[2:60, 1] ~ dmnorm(zero[2:60], K1)


    ## smoothing parameter priors
    for(ii in 1:1) {
      lambda[ii] ~ dgamma(.05,.005)
      rho[ii] <- log(lambda[ii])
    }

    ### YEARS 2-t
    for(jj in 1:60){
      for(tt in 2:nYears){
        b[jj, tt] ~ dnorm(b[jj, tt - 1], tau.gam)
      }
    }

    tau.gam <- pow(sigma.gam, -2)
    sigma.gam ~ dunif(0, 10)


    ### Likelihood
    for (ii in 1:nRoutes) {
      ## Detection probability
        p[ii, 1, 1] <- 0
        logit(p[ii, 2, 1]) <- alpha0 + alpha1 * Xp[ii, 1] + alpha2 * nov[ii, 1] +  alpha3 * twedt[ii, 1] + omega[obs[ii, 1]]

      ## Initial occupancy
        logit(psi[ii, 1]) <- inprod(X[ii,], b[, 1]) + inprod(Xclim[ii,,1], beta)
        z[ii, 1] ~ dbern(psi[ii, 1])
        z.new[ii, 1] ~ dbern(psi[ii, 1])

      ##  y: local availability at stop 1 -- 0 = locally unavailable,  1 = locally available
        y[ii, 1, 1] ~ dbern(pi[1])
        y.new[ii, 1, 1] ~ dbern(pi[1])
        h[ii, 1, 1] ~ dbern(p[ii, (z[ii, 1] * y[ii, 1, 1] + 1), 1])
        h.new[ii, 1, 1] ~ dbern(p[ii, (z.new[ii, 1] * y.new[ii, 1, 1] + 1), 1])

      ## Availability at stops 2-nStops
        for (jj in 2:nStops) {
          y[ii, jj, 1] ~ dbern(xpsi[1, (z[ii, 1] * y[ii, jj - 1, 1] + 1)])
          y.new[ii, jj, 1] ~ dbern(xpsi[1, (z.new[ii, 1] * y.new[ii, jj - 1, 1] + 1)])
          h[ii, jj, 1] ~ dbern(p[ii, (z[ii, 1] * y[ii, jj, 1] + 1), 1])
          h.new[ii, jj, 1] ~ dbern(p[ii, (z.new[ii, 1] * y.new[ii, jj, 1] + 1), 1])
        } # jj


      for(tt in 2:nYears){
        p[ii, 1, tt] <- 0
        logit(p[ii, 2, tt]) <- alpha0 + alpha1 * Xp[ii, tt] + alpha2 * nov[ii, tt] + alpha3 * twedt[ii, tt] +omega[obs[ii, tt]]

        ## Occupancy
        logit(psi[ii, tt]) <- inprod(X[ii,], b[, tt]) + inprod(Xclim[ii,,tt], beta)
        z[ii, tt] ~ dbern(psi[ii, tt])
        z.new[ii, tt] ~ dbern(psi[ii, tt])

        ##  y: local availability at stop 1 -- 0 = locally unavailable,  1 = locally available
        y[ii, 1, tt] ~ dbern(pi[tt])
        y.new[ii, 1, tt] ~ dbern(pi[tt])
        h[ii, 1, tt] ~ dbern(p[ii, (z[ii, tt] * y[ii, 1, tt] + 1), tt])
        h.new[ii, 1, tt] ~ dbern(p[ii, (z.new[ii, tt] * y.new[ii, 1, tt] + 1), tt])

        ## Availability at stops 2-nStops
        for (jj in 2:nStops) {
          y[ii, jj, tt] ~ dbern(xpsi[tt, (z[ii, tt] * y[ii, jj - 1, tt] + 1)])
          y.new[ii, jj, tt] ~ dbern(xpsi[tt, (z.new[ii, tt] * y.new[ii, jj - 1, tt] + 1)])
          h[ii, jj, tt] ~ dbern(p[ii, (z[ii, tt] * y[ii, jj, tt] + 1), tt])
          h.new[ii, jj, tt] ~ dbern(p[ii, (z.new[ii, tt] * y.new[ii, jj, tt] + 1), tt])
        } # jj
      } # tt
    } # ii

    # for(tt in 1:nYears){
    #   sum.z[tt] <- sum(z[1:nRoutes, tt])
    #   sum.z.new[tt] <- sum(z.new[1:nRoutes, tt])
    #
    #   fit.z[tt] <- pow((sum.z[tt] - sum(psi[1:nRoutes, tt])), 2)/sum(psi[1:nRoutes, tt])
    #   fit.z.new[tt] <- pow((sum.z.new[tt] - sum(psi[1:nRoutes, tt])), 2)/sum(psi[1:nRoutes, tt])
    #
    #   mu.y[tt] <- mean(y[1:nRoutes, 1:nStops, tt])
    #   mu.y.new[tt] <- mean(y.new[1:nRoutes, 1:nStops, tt])
    #
    #   fit.y[tt] <- pow((mu.y[tt] - pi[tt]), 2)/pi[tt]
    #   fit.y.new[tt] <- pow((mu.y.new[tt] - pi[tt]), 2)/pi[tt]
    # }

} # End model
    ", fill=TRUE)
sink()
