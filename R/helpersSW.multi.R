# This file contains functions for power analysis and sample size calculations
# for stepped wedge cluster randomized trials (SW-CRTs). It supports testing
# heterogeneous treatment effects (HTE) for a multivariate effect modifier and
# average treatment effects (ATE). The methods accommodate both cross-sectional
# and closed-cohort stepped wedge designs, as detailed in the manuscript:
# [Planning Stepped Wedge Cluster Randomized Trials to Detect Treatment Effect
# Heterogeneity].
#
# For each design, the file provides four functions to perform sample size and
# power calculations tailored to different estimands.


################################################################################
######################### Heterogeneous treatment effect########################
################################################################################

####### cross-sectional design #######

## Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
## under the cross-sectional design, HTE test
ss_CS_HTE.multi <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate cross-correlation matrix (elements in [0,1))
  # rho1: between-period covariate cross-correlation matrix (elements in [0,1))
  # alpha: type I error rate
  # beta: type II error rate
  #
  # Output:
  # I: predicted number of clusters

  dimen <- length(eff)
  nca <- 0
  power <- 0
  while (power < 1-beta){
    nca <- nca + 1
    I <- (J-1)*nca

    #make the staggered table
    trtSeq <- matrix(0, ncol=J, nrow=J-1)
    trtSeq[upper.tri(trtSeq)] <- 1

    W_noN <- NULL
    for (i in 1:(J-1)){
      W_noN <- c(W_noN, rep(trtSeq[i,], nca))
    }
    W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)

    # Follow the definition in Li et al. 2018
    # W_ij, cluster-period as the unit
    U <- sum(W_noN)
    W <- sum(colSums(W_noN_matrix)^2)
    V <- sum(rowSums(W_noN_matrix)^2)

    lambda1 <- 1-alpha0
    lambda2 <- 1+(N-1)*alpha0-N*alpha1
    lambda3 <- 1+(N-1)*alpha0+(J-1)*N*alpha1

    sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
    Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
    Z1 <- Upsilon-rho0
    Z2 <- Upsilon+(N-1)*rho0-N*rho1
    Z3 <- Upsilon+(N-1)*rho0+(J-1)*N*rho1
    Thetacs <- J*(N-1)/lambda1*Z1+(J-1)/lambda2*Z2+1/lambda3*Z3

    Var4 <- sigma2y * I*J^2/(I*U-W)*sqrtinvOmega %*% solve((U^2+I*J*U-J*W-I*V)/(I*U-W)*(1/lambda2-1/lambda3)*(Z3-Z2)+J*Thetacs) %*% sqrtinvOmega
    chi_critical <- qchisq(1 - alpha, df = dimen)
    noncentrality <- t(eff) %*% solve(Var4) %*% eff

    power <- pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE)
  }
  return(I)
}


## Function to calculate the predicted power given the number of clusters
## under the cross-sectional design, HTE test
power_CS_HTE.multi <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, I){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate cross-correlation matrix (elements in [0,1))
  # rho1: between-period covariate cross-correlation matrix (elements in [0,1))
  # alpha: type I error rate
  # I: predicted number of clusters
  #
  # Output:
  # power: actual predicted power

  dimen <- length(eff)
  nca <- I/(J-1)

  #make the staggered table
  trtSeq <- matrix(0, ncol=J, nrow=J-1)
  trtSeq[upper.tri(trtSeq)] <- 1

  W_noN <- NULL
  for (i in 1:(J-1)){
    W_noN <- c(W_noN, rep(trtSeq[i,], nca))
  }
  W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)

  # Follow the definition in Li et al. 2018
  # W_ij, cluster-period as the unit
  U <- sum(W_noN)
  W <- sum(colSums(W_noN_matrix)^2)
  V <- sum(rowSums(W_noN_matrix)^2)

  lambda1 <- 1-alpha0
  lambda2 <- 1+(N-1)*alpha0-N*alpha1
  lambda3 <- 1+(N-1)*alpha0+(J-1)*N*alpha1

  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  Z1 <- Upsilon-rho0
  Z2 <- Upsilon+(N-1)*rho0-N*rho1
  Z3 <- Upsilon+(N-1)*rho0+(J-1)*N*rho1
  Thetacs <- J*(N-1)/lambda1*Z1+(J-1)/lambda2*Z2+1/lambda3*Z3

  Var4 <- sigma2y * I*J^2/(I*U-W)*sqrtinvOmega %*% solve((U^2+I*J*U-J*W-I*V)/(I*U-W)*(Z3-Z2)*(1/lambda2-1/lambda3)+J*Thetacs) %*% sqrtinvOmega
  chi_critical <- qchisq(1 - alpha, df = dimen)
  noncentrality <- t(eff) %*% solve(Var4) %*% eff

  power <- pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE)
  return(power)
}



####### closed-cohort design ######

## Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
## under the closed-cohort design, HTE test
ss_CC_HTE.multi <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate cross-correlation matrix (elements in [0,1))
  # alpha: type I error rate
  # beta: type II error rate
  #
  # Output:
  # I: predicted number of clusters

  dimen <- length(eff)
  nca <- 0
  power <- 0
  while (power < 1-beta){
    nca <- nca + 1
    I <- (J-1)*nca

    trtSeq <- matrix(0, ncol=J, nrow=J-1)
    trtSeq[upper.tri(trtSeq)] <- 1

    W_noN <- NULL
    for (i in 1:(J-1)){
      W_noN <- c(W_noN, rep(trtSeq[i,], nca))
    }
    W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)

    U <- sum(W_noN)
    W <- sum(colSums(W_noN_matrix)^2)
    V <- sum(rowSums(W_noN_matrix)^2)

    tau1 <- 1-alpha0+alpha1-alpha2
    tau2 <- 1-alpha0-(J-1)*(alpha1-alpha2)
    tau3 <- 1+(N-1)*(alpha0-alpha1)-alpha2
    tau4 <- 1+(N-1)*alpha0+(J-1)*(N-1)*alpha1+(J-1)*alpha2

    sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
    Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
    H1 <- Upsilon-rho0
    H2 <- Upsilon+(N-1)*rho0
    Thetacc <- (N-1)/tau2*H1+1/tau4*H2

    Var4 <- sigma2y * (I*J)/(I*U-W)*sqrtinvOmega %*% solve((U^2+I*J*U-J*W-I*V)/(I*U-W)*((1/tau3-1/tau4)*H2+(N-1)*(1/tau1-1/tau2)*H1)+J*Thetacc) %*% sqrtinvOmega
    chi_critical <- qchisq(1 - alpha, df = dimen)
    noncentrality <- t(eff) %*% solve(Var4) %*% eff

    power <- pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE)
  }
  return(I)
}


## Function to calculate the predicted power given the number of clusters
## under the closed-cohort design, HTE test
power_CC_HTE.multi <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, I){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: covariate covariance matrix (symmetric positive definite)
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate cross-correlation matrix (elements in [0,1))
  # alpha: type I error rate
  # I: predicted number of clusters
  #
  # Output:
  # power: actual predicted power

  dimen <- length(eff)
  nca <- I/(J-1)

  trtSeq <- matrix(0, ncol=J, nrow=J-1)
  trtSeq[upper.tri(trtSeq)] <- 1

  W_noN <- NULL
  for (i in 1:(J-1)){
    W_noN <- c(W_noN, rep(trtSeq[i,], nca))
  }
  W_noN_matrix <- matrix(W_noN, byrow=T, nrow=I)

  U <- sum(W_noN)
  W <- sum(colSums(W_noN_matrix)^2)
  V <- sum(rowSums(W_noN_matrix)^2)

  tau1 <- 1-alpha0+alpha1-alpha2
  tau2 <- 1-alpha0-(J-1)*(alpha1-alpha2)
  tau3 <- 1+(N-1)*(alpha0-alpha1)-alpha2
  tau4 <- 1+(N-1)*alpha0+(J-1)*(N-1)*alpha1+(J-1)*alpha2

  sqrtinvOmega <- diag(sqrt(diag(sigma2x))^(-1))
  Upsilon <- sqrtinvOmega %*% sigma2x %*% sqrtinvOmega
  H1 <- Upsilon-rho0
  H2 <- Upsilon+(N-1)*rho0
  Thetacc <- (N-1)/tau2*H1+1/tau4*H2

  Var4 <- sigma2y * (I*J)/(I*U-W)*sqrtinvOmega %*% solve((U^2+I*J*U-J*W-I*V)/(I*U-W)*((1/tau3-1/tau4)*H2+(N-1)*(1/tau1-1/tau2)*H1)+J*Thetacc) %*% sqrtinvOmega
  chi_critical <- qchisq(1 - alpha, df = dimen)
  noncentrality <- t(eff) %*% solve(Var4) %*% eff

  power <- pchisq(chi_critical, df = dimen, ncp = noncentrality, lower.tail = FALSE)
  return(power)
}

