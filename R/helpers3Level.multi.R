# This file contains functions for conducting power analysis and calculating
# sample sizes for testing HTE under three levels of randomization for a cluster
# randomized trial when the dimension of effect modifiers is greater than or
# equal to 2, as described in the manuscript: [Planning Three-Level Cluster
# Randomized Trials to Assess Treatment Effect Heterogeneity].
#
# Under each of the 3 randomization scenarios, 2 functions were created for
# sample size/power computation.

#############################################################################################
## Sample size calculations and power analysis under cluster-level (level 3) randomization ##
#############################################################################################

### Testing HTE

## Sample size (number of clusters) calculation for testing HTE under cluster-level (level 3) randomization
L3_HTE_nc.multi <- function(eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi, max.ss=10^6){
  # Argument:
  # eff: a vector of effect sizes (beta4)
  # rho0: within-subcluster covariate cross-correlation matrix (elements in [0,1))
  # rho1: between-subcluster covariate cross-correlation matrix (elements in [0,1))
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # sigma2y: outcome variance
  # power: required power level (= 1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated
  # max.ss: maximum number of clusters to be searched

  dimen <- length(eff)
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  kappa2 <- 1/lambda1-(lambda2-lambda1)/m/lambda1/lambda2-(lambda3-lambda2)/ns/m/lambda2/lambda3
  kappa1 <- -(ns-1)*(lambda3-lambda2)/ns/lambda2/lambda3
  kappa0 <- -(m-1)*(lambda2-lambda1)/m/lambda1/lambda2-(m-1)*(lambda3-lambda2)/ns/m/lambda2/lambda3
  sigma2w <- pi*(1-pi)
  sigma24 <- sigma2y/sigma2w*sqrtinvOmega %*% solve(kappa2*Upsilon+kappa1*rho1+kappa0*rho0) %*% sqrtinvOmega

  # Compute the critical chi-squared value
  chi_critical <- qchisq(1 - alpha, df = dimen)

  # Define the power equation as a function of nc
  power_function <- function(nc) {
    noncentrality <- nc * t(eff) %*% solve(sigma24) %*% eff
    pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE) - power
  }

  nc <- uniroot(power_function, interval = c(0, max.ss))$root
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing HTE under cluster-level (level 3) randomization
L3_HTE_power.multi <- function(nc, eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta4)
  # rho0: within-subcluster covariate cross-correlation matrix (elements in [0,1))
  # rho1: between-subcluster covariate cross-correlation matrix (elements in [0,1))
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  dimen <- length(eff)
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  kappa2 <- 1/lambda1-(lambda2-lambda1)/m/lambda1/lambda2-(lambda3-lambda2)/ns/m/lambda2/lambda3
  kappa1 <- -(ns-1)*(lambda3-lambda2)/ns/lambda2/lambda3
  kappa0 <- -(m-1)*(lambda2-lambda1)/m/lambda1/lambda2-(m-1)*(lambda3-lambda2)/ns/m/lambda2/lambda3
  sigma2w <- pi*(1-pi)
  sigma24 <- sigma2y/sigma2w*sqrtinvOmega %*% solve(kappa2*Upsilon+kappa1*rho1+kappa0*rho0) %*% sqrtinvOmega

  chi_critical <- qchisq(1 - alpha, df = dimen)
  noncentrality <- nc * t(eff) %*% solve(sigma24) %*% eff
  power <- pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE)
  return (power)
}




#######################################################################################
## Sample size and power calculations under subcluster-level (level 2) randomization ##
#######################################################################################

## Sample size (number of clusters) calculation for testing HTE under subcluster-level (level 2) randomization
L2_HTE_nc.multi <- function(eff, rho0, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi, max.ss=10^6){
  # Argument:
  # eff: effect size (beta4)
  # rho0: within-subcluster covariate cross-correlation matrix (elements in [0,1))
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # sigma2y: outcome variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated
  # max.ss: maximum number of clusters to be searched

  dimen <- length(eff)
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  kappa2 <- 1/lambda1-(lambda2-lambda1)/m/lambda1/lambda2
  kappa0 <- -(m-1)*(lambda2-lambda1)/m/lambda1/lambda2
  sigma2w <- pi*(1-pi)
  sigma24 <- sigma2y/sigma2w*sqrtinvOmega %*% solve(kappa2*Upsilon+kappa0*rho0) %*% sqrtinvOmega

  # Compute the critical chi-squared value
  chi_critical <- qchisq(1 - alpha, df = dimen)

  # Define the power equation as a function of nc
  power_function <- function(nc) {
    noncentrality <- nc * t(eff) %*% solve(sigma24) %*% eff
    pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE) - power
  }

  nc <- uniroot(power_function, interval = c(0, max.ss))$root
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing HTE under subcluster-level (level 2) randomization
L2_HTE_power.multi <- function(nc, eff, rho0, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta4)
  # rho0: within-subcluster covariate cross-correlation matrix (elements in [0,1))
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  dimen <- length(eff)
  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  kappa2 <- 1/lambda1-(lambda2-lambda1)/m/lambda1/lambda2
  kappa0 <- -(m-1)*(lambda2-lambda1)/m/lambda1/lambda2
  sigma2w <- pi*(1-pi)
  sigma24 <- sigma2y/sigma2w*sqrtinvOmega %*% solve(kappa2*Upsilon+kappa0*rho0) %*% sqrtinvOmega

  chi_critical <- qchisq(1 - alpha, df = dimen)
  noncentrality <- nc * t(eff) %*% solve(sigma24) %*% eff
  power <- pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE)
  return (power)
}




#######################################################################################
## Sample size and power calculations under individual-level (level 1) randomization ##
#######################################################################################

## Sample size (number of clusters) calculation for testing HTE under individual-level (level 1) randomization
L1_HTE_nc.multi <- function(eff, alpha0, sigma2x, sigma2y, power, alpha, m, ns, pi, max.ss=10^6){
  # Argument:
  # eff: effect size (beta4)
  # alpha0: within-subcluster outcome ICC
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # sigma2y: outcome variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated
  # max.ss: maximum number of clusters to be searched

  dimen <- length(eff)
  lambda1 <- 1-alpha0
  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  sigma2w <- pi*(1-pi)
  sigma24 <- sigma2y/sigma2w*lambda1*sqrtinvOmega %*% solve(Upsilon) %*% sqrtinvOmega

  # Compute the critical chi-squared value
  chi_critical <- qchisq(1 - alpha, df = dimen)

  # Define the power equation as a function of nc
  power_function <- function(nc) {
    noncentrality <- nc * t(eff) %*% solve(sigma24) %*% eff
    pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE) - power
  }

  nc <- uniroot(power_function, interval = c(0, max.ss))$root
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing HTE under individual-level (level 1) randomization
L1_HTE_power.multi <- function(nc, eff, alpha0, sigma2x, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta4)
  # alpha0: within-subcluster outcome ICC
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  dimen <- length(eff)
  lambda1 <- 1-alpha0
  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  sigma2w <- pi*(1-pi)
  sigma24 <- sigma2y/sigma2w*lambda1*sqrtinvOmega %*% solve(Upsilon) %*% sqrtinvOmega

  chi_critical <- qchisq(1 - alpha, df = dimen)
  noncentrality <- nc * t(eff) %*% solve(sigma24) %*% eff
  power <- pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE)
  return (power)
}

