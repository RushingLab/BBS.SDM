#' get.betas
#'
#' Retrieve beta values and var-cov matrix from top model

GetBetas <- function(alpha, path) {

  # read top model and pull out betas and vc matrix
  inits.out <- scan(paste(path, alpha, "inits.out", sep = "/"), what = 'character', sep = '\n', quiet = T)

  jj <- grep('std.error', inits.out)
  jj.end <- grep('Variance-Covariance Matrix of Untransformed', inits.out)
  betas <- inits.out[(jj + 1):(jj.end - 1)]

  psi.est <- betas[grep('psi', betas)]
  loc.per <- regexpr("psi", psi.est)
  psi.names <- trimws(substr(psi.est, loc.per + 4, loc.per + 11))
  psi.names <- gsub("[0-9]", "", psi.names)[-1]
  psi.names <- gsub("_$", "", psi.names)
  psi.betas <- as.numeric(substr(psi.est, 41, 50))
  psi.se <- as.numeric(substr(psi.est, 54, 59))
  names(psi.betas) <- trimws(c("int", psi.names))
  names(psi.se) <- trimws(c("int", psi.names))
  psi.na <- is.na(psi.se)
  psi.betas[psi.na] <- 0
  psi.se[psi.na] <- 1

  gam.est <- betas[grep('gam', betas)]
  loc.per <- regexpr("gam", gam.est)
  gam.names <- substr(gam.est, loc.per + 5, loc.per + 12)
  gam.names <- gsub("[0-9]","",gam.names)[-1]
  gam.names <- gsub("_$", "", gam.names)
  gam.betas <- as.numeric(substr(gam.est, 41, 50))
  gam.se <- as.numeric(substr(gam.est, 54, 59))
  names(gam.betas) <- trimws(c("int", gam.names))
  names(gam.se) <- trimws(c("int", gam.names))
  gam.na <- is.na(gam.se)
  gam.betas[gam.na] <- 0
  gam.se[gam.na] <- 1

  eps.est <- betas[grep('eps', betas)]
  loc.per <- regexpr("eps", eps.est)
  eps.names <- substr(eps.est, loc.per + 5, loc.per + 12)
  eps.names <- gsub("[0-9]","",eps.names)[-1]
  eps.names <- gsub("_$", "", eps.names)
  eps.betas <- as.numeric(substr(eps.est, 41, 50))
  eps.se <- as.numeric(substr(eps.est, 54, 59))
  names(eps.betas) <- trimws(c("int", eps.names))
  names(eps.se) <- trimws(c("int", eps.names))
  eps.na <- is.na(eps.se)
  eps.betas[eps.na] <- 0
  eps.se[eps.na] <- 1

  p.betas <- betas[grep('^D', betas)]
  p.names <- substr(p.betas, 11, 18)[-1]
  p.betas <- as.numeric(substr(p.betas, 41, 50))
  names(p.betas) <- trimws(c("int", p.names))

  th.betas <- betas[grep('th', betas)]
  th.names <- substr(th.betas, 6, 8)
  th.betas <- as.numeric(substr(th.betas, 41, 50))
  names(th.betas) <- th.names


  return(list(psi.betas = psi.betas[-1], psi.se = psi.se[-1],
              gam.betas = gam.betas[-1], gam.se = gam.se[-1],
              eps.betas = eps.betas[-1], eps.se = eps.se[-1],
              p.betas = p.betas, th.betas = th.betas))
}
