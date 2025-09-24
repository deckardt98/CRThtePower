# This file contains functions for conducting power analysis and calculating
# sample sizes for testing HTE, ATE, and unadjusted ATE under three levels of
# randomization for a cluster randomized trial, as described in the manuscript:
# [Planning Three-Level Cluster Randomized Trials to Assess Treatment Effect
# Heterogeneity].
#
# Under each of the 3 randomization scenarios, 2x3 functions were created for
# sample size/power computation under different tests.

#############################################################################################
## Sample size calculations and power analysis under cluster-level (level 3) randomization ##
#############################################################################################

### Testing HTE

## Sample size (number of clusters) calculation for testing HTE under cluster-level (level 3) randomization
L3_HTE_nc <- function(eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi){
  # Argument:
  # eff: effect size (beta4)
  # rho0: within-subcluster covariate ICC
  # rho1: between-subcluster covariate ICC
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # power: required power level (= 1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  l1 <- 1-rho0                                      # zeta1
  l2 <- 1+(m-1)*rho0-m*rho1                         # zeta2
  l3 <- 1+(m-1)*rho0+m*(ns-1)*rho1                  # zeta3
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*(m-1)*lambda2*lambda3*l1+(ns-1)*lambda1*lambda3*l2+lambda1*lambda2*l3))
  nc <- (qnorm(1-alpha/2) + qnorm(power))^2*sigma24/eff^2
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing HTE under cluster-level (level 3) randomization
L3_HTE_power <- function(nc, eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta4)
  # rho0: within-subcluster covariate ICC
  # rho1: between-subcluster covariate ICC
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  l1 <- 1-rho0
  l2 <- 1+(m-1)*rho0-m*rho1
  l3 <- 1+(m-1)*rho0+m*(ns-1)*rho1
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*(m-1)*lambda2*lambda3*l1+(ns-1)*lambda1*lambda3*l2+lambda1*lambda2*l3))
  power <- pnorm( sqrt(nc*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}

### Testing ATE

## Sample size (number of clusters) calculation for testing ATE under cluster-level (level 3) randomization (based on t-approximation)
L3_ATE_nc <- function(eff, alpha0, alpha1, sigma2y, power, alpha, m, ns, pi){
  # Argument:
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2y: outcome variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  nVar <- sigma2y*lambda3/(pi*(1-pi)*ns*m)

  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <- denom
  current_power <- 0
  while (current_power < power){
    nc <- nc+denom
    current_power <- pt(qt(1-alpha/2, nc-denom), nc-denom, ncp=eff/sqrt(nVar/nc), lower.tail = F) + pt(qt(alpha/2, nc-denom), nc-denom, ncp=eff/sqrt(nVar/nc))
  }
  return(nc)
}

## Power calculation for testing ATE under cluster-level (level 3) randomization (based on t-approximation)
L3_ATE_power <- function(nc, eff, alpha0, alpha1, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  nVar <- sigma2y*lambda3/(pi*(1-pi)*ns*m)
  power <- pt(qt(1-alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc), lower.tail = F) + pt(qt(alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc))
  return(power)
}

### Testing unadjusted ATE

## Sample size (number of clusters) calculation for testing "unadjusted" ATE under cluster-level (level 3) randomization (based on t-approximation)
L3_ATE2_nc <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, power, alpha, m, ns, pi){
  # Argument:
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # rho0: within-subcluster covariate ICC
  # rho1: between-subcluster covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda3 <- 1+(m-1)*unadjusted_alpha0+m*(ns-1)*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x

  nVar <- unadjusted_sigma2y*lambda3/(pi*(1-pi)*ns*m)

  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <- denom
  current_power <- 0
  while (current_power < power){
    nc <- nc+denom
    current_power <- pt(qt(1-alpha/2, nc-denom), nc-denom, ncp=eff/sqrt(nVar/nc), lower.tail = F) + pt(qt(alpha/2, nc-denom), nc-denom, ncp=eff/sqrt(nVar/nc))
  }
  return(nc)
}

## Power calculation for testing "unadjusted" ATE under cluster-level (level 3) randomization (based on t-approximation)
L3_ATE2_power <- function(nc, beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # rho0: within-subcluster covariate ICC
  # rho1: between-subcluster covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda3 <- 1+(m-1)*unadjusted_alpha0+m*(ns-1)*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x

  nVar <- unadjusted_sigma2y*lambda3/(pi*(1-pi)*ns*m)
  power <- pt(qt(1-alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc), lower.tail = F) + pt(qt(alpha/2, nc-2), nc-2, ncp=eff/sqrt(nVar/nc))
  return(power)
}




#######################################################################################
## Sample size and power calculations under subcluster-level (level 2) randomization ##
#######################################################################################

## Sample size (number of clusters) calculation for testing HTE under subcluster-level (level 2) randomization
L2_HTE_nc <- function(eff, rho0, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi){
  # Argument:
  # eff: effect size (beta4)
  # rho0: within-subcluster covariate ICC
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  zeta1 <- 1-rho0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*m*lambda2*lambda3-ns*lambda3*(lambda2-lambda1)*(1+(m-1)*rho0)))
  nc <- (qnorm(1-alpha/2) + qnorm(power))^2*sigma24/eff^2
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing HTE under subcluster-level (level 2) randomization
L2_HTE_power <- function(nc, eff, rho0, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta4)
  # rho0: within-subcluster covariate ICC
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  lambda3 <- 1+(m-1)*alpha0+m*(ns-1)*alpha1
  zeta1 <- 1-rho0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1*lambda2*lambda3)/(sigma2w*sigma2x*(ns*m*lambda2*lambda3-ns*lambda3*(lambda2-lambda1)*(1+(m-1)*rho0)))
  power <- pnorm( sqrt(nc*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}



## Sample size (number of clusters) calculation for testing ATE under subcluster-level (level 2) randomization
L2_ATE_nc <- function(eff, alpha0, alpha1, sigma2y, power, alpha, m, ns, pi){
  # Argument:
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2y: outcome variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  nc <- (qnorm(1-alpha/2) + qnorm(power))^2*sigma2y*lambda2/(pi*(1-pi)*ns*m*eff^2)
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return(nc)
}

## Power calculation for testing ATE under subcluster-level (level 2) randomization
L2_ATE_power <- function(nc, eff, alpha0, alpha1, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda2 <- 1+(m-1)*alpha0-m*alpha1
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(sigma2y*lambda2))-qnorm(1-alpha/2) )
  return (power)
}



## Sample size (number of clusters) calculation for testing "unadjusted" ATE under subcluster-level (level 2) randomization
L2_ATE2_nc <- function(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, power, alpha, m, ns, pi){
  # Argument:
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # rho0: within-subcluster covariate ICC
  # rho1: between-subcluster covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda2 <- 1 + (m-1)*unadjusted_alpha0 - m*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  nc <- (qnorm(1-alpha/2) + qnorm(power))^2*unadjusted_sigma2y*lambda2/(pi*(1-pi)*ns*m*eff^2)
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing "unadjusted" ATE under subcluster-level (level 2) randomization
L2_ATE2_power <- function(nc, beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # alpha1: between-subcluster outcome ICC
  # rho0: within-subcluster covariate ICC
  # rho1: between-subcluster covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0
  unadjusted_alpha1 <- w*alpha1 + (1-w)*rho1
  lambda2 <- 1 + (m-1)*unadjusted_alpha0 - m*unadjusted_alpha1
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(unadjusted_sigma2y*lambda2))-qnorm(1-alpha/2) )
  return (power)
}




#######################################################################################
## Sample size and power calculations under individual-level (level 1) randomization ##
#######################################################################################

## Sample size (number of clusters) calculation for testing HTE under individual-level (level 1) randomization
L1_HTE_nc <- function(eff, alpha0, sigma2x, sigma2y, power, alpha, m, ns, pi){
  # Argument:
  # eff: effect size (beta4)
  # alpha0: within-subcluster outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1)/(sigma2w*sigma2x*ns*m)
  nc <- (qnorm(1-alpha/2) + qnorm(power))^2*sigma24/eff^2
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing HTE under individual-level (level 1) randomization
L1_HTE_power <- function(nc, eff, alpha0, sigma2x, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta4)
  # alpha0: within-subcluster outcome ICC
  # sigma2x: covariate variance
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  sigma2w <- pi*(1-pi)
  sigma24 <- (sigma2y*lambda1)/(sigma2w*sigma2x*ns*m)
  power <- pnorm( sqrt(nc*eff^2/sigma24)-qnorm(1-alpha/2) )
  return (power)
}



## Sample size (number of clusters) calculation for testing ATE under individual-level (level 1) randomization
L1_ATE_nc <- function(eff, alpha0, sigma2y, power, alpha, m, ns, pi){
  # Argument:
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # sigma2y: outcome variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  nc <- (qnorm(1-alpha/2) + qnorm(power))^2*sigma2y*lambda1/(pi*(1-pi)*ns*m*eff^2)
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return(nc)
}

## Power calculation for testing ATE under individual-level (level 1) randomization
L1_ATE_power <- function(nc, eff, alpha0, sigma2y, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # sigma2y: outcome variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  lambda1 <- 1-alpha0
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(sigma2y*lambda1))-qnorm(1-alpha/2) )
  return (power)
}



## Sample size (number of clusters) calculation for testing "unadjusted" ATE under individual-level (level 1) randomization
L1_ATE2_nc <- function(beta3, beta4, eff, alpha0, rho0, sigma2y, sigma2x, power, alpha, m, ns, pi){
  # Argument:
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # rho0: within-subcluster covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # power: required power level (1-type II error)
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0

  lambda1 <- 1 - unadjusted_alpha0
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  nc <- (qnorm(1-alpha/2) + qnorm(power))^2*unadjusted_sigma2y*lambda1/(pi*(1-pi)*ns*m*eff^2)
  frac <- fractions(pi)
  denom <-  as.numeric(strsplit(attr(frac,"fracs"),"/")[[1]][2])
  nc <-   denom*ceiling(nc/denom)
  return (nc)
}

## Power calculation for testing "unadjusted" ATE under individual-level (level 1) randomization
L1_ATE2_power <- function(nc, beta3, beta4, eff, alpha0, rho0, sigma2y, sigma2x, alpha, m, ns, pi){
  # Argument:
  # nc: number of clusters
  # beta3: true covariate effect
  # beta4: true effect of treatment effect heterogeneity
  # eff: effect size (beta2)
  # alpha0: within-subcluster outcome ICC
  # rho0: within-subcluster covariate ICC
  # sigma2y: outcome variance
  # sigma2x: covariate variance
  # alpha: type I error
  # m: number of patients under each subcluster
  # ns: number of subclusters under each cluster
  # pi: proportion of treated

  Q <- beta3^2 + beta4^2*pi + 2*beta3*beta4*pi
  w <- sigma2y/(sigma2y+Q*sigma2x)
  unadjusted_alpha0 <- w*alpha0 + (1-w)*rho0

  lambda1 <- 1 - unadjusted_alpha0
  unadjusted_sigma2y <- sigma2y + Q*sigma2x
  power <- pnorm( sqrt(nc*pi*(1-pi)*ns*m*eff^2/(unadjusted_sigma2y*lambda1))-qnorm(1-alpha/2) )
  return (power)
}
