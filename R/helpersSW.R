# This file contains functions for power analysis and sample size calculations
# for stepped wedge cluster randomized trials (SW-CRTs). It supports testing
# heterogeneous treatment effects (HTE) for a univariate effect modifier and
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

####### cross-sectional design ######

## Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
## under the cross-sectional design, HTE test
ss_CS_HTE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  #
  # Output:
  # I: predicted number of clusters

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
    zeta1 <- 1-rho0
    zeta2 <- 1+(N-1)*rho0-N*rho1
    zeta3 <- 1+(N-1)*rho0+(J-1)*N*rho1

    Var4 <- sigma2y/sigma2x * I*J^2/
      ( (I*U-W)*J*(J*(N-1)*zeta1/lambda1 + (J-1)*zeta2/lambda2 + zeta3/lambda3) + (U^2+I*J*U-J*W-I*V)*(1/lambda2-1/lambda3)*(zeta3-zeta2) )

    power <- pnorm( sqrt(eff^2/Var4)-qnorm(1-alpha/2) )
  }
  return(I)
}


## Function to calculate the predicted power given the number of clusters
## under the cross-sectional design, HTE test
power_CS_HTE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, I){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha: type I error rate
  # I: predicted number of clusters
  #
  # Output:
  # power: actual predicted power

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
  zeta1 <- 1-rho0
  zeta2 <- 1+(N-1)*rho0-N*rho1
  zeta3 <- 1+(N-1)*rho0+(J-1)*N*rho1

  Var4 <- sigma2y/sigma2x * I*J^2/
    ( (I*U-W)*J*(J*(N-1)*zeta1/lambda1 + (J-1)*zeta2/lambda2 + zeta3/lambda3) + (U^2+I*J*U-J*W-I*V)*(1/lambda2-1/lambda3)*(zeta3-zeta2) )

  power <- pnorm( sqrt(eff^2/Var4)-qnorm(1-alpha/2) )
  return(power)
}



####### closed-cohort design ######

## Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
## under the closed-cohort design, HTE test
ss_CC_HTE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  #
  # Output:
  # I: predicted number of clusters

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

    eta1 <- 1-rho0
    eta2 <- 1+(N-1)*rho0

    kappa1 <- (N-1)*eta1/tau2 + eta2/tau4
    kappa3 <- (1/tau3-1/tau4)*eta2 + (N-1)*(1/tau1-1/tau2)*eta1

    Var4 <- sigma2y/sigma2x * (I*J)/
      ( (I*U-W)*J*kappa1 + (U^2+I*J*U-J*W-I*V)*kappa3 )

    power <- pnorm( sqrt(eff^2/Var4)-qnorm(1-alpha/2) )
  }
  return(I)
}


## Function to calculate the predicted power given the number of clusters
## under the closed-cohort design, HTE test
power_CC_HTE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, I){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC
  # alpha: type I error rate
  # I: predicted number of clusters
  #
  # Output:
  # power: actual predicted power

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

  eta1 <- 1-rho0
  eta2 <- 1+(N-1)*rho0

  kappa1 <- (N-1)*eta1/tau2 + eta2/tau4
  kappa3 <- (1/tau3-1/tau4)*eta2 + (N-1)*(1/tau1-1/tau2)*eta1

  Var4 <- sigma2y/sigma2x * (I*J)/
    ( (I*U-W)*J*kappa1 + (U^2+I*J*U-J*W-I*V)*kappa3 )

  power <- pnorm( sqrt(eff^2/Var4)-qnorm(1-alpha/2) )
  return(power)
}

################################################################################
############################# Average treatment effect##########################
################################################################################

####### cross-sectional design ######

## Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
## under the cross-sectional design, ATE test
ss_CS_ATE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  #
  # Output:
  # I: predicted number of clusters

  nca <- 1
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

    U <- sum(W_noN)
    W <- sum(colSums(W_noN_matrix)^2)
    V <- sum(rowSums(W_noN_matrix)^2)

    lambda1 <- 1-alpha0
    lambda2 <- 1+(N-1)*alpha0-N*alpha1
    lambda3 <- 1+(N-1)*alpha0+(J-1)*N*alpha1

    Var <- sigma2y/N * (I*J*lambda2*lambda3)/
      ( (U^2+I*J*U-J*W-I*V)*lambda3-(U^2-I*V)*lambda2 )

    power <- pt(qt(1-alpha/2, I-2), I-2, ncp=eff/sqrt(Var), lower.tail = F)
  }
  return(I)
}


## Function to calculate the predicted power given the number of clusters
## under the cross-sectional design, ATE test
power_CS_ATE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, rho0, rho1, alpha=0.05, I){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha: type I error rate
  # I: predicted number of clusters
  #
  # Output:
  # power: actual predicted power

  nca <- I/(J-1)

  #make the staggered table
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

  lambda1 <- 1-alpha0
  lambda2 <- 1+(N-1)*alpha0-N*alpha1
  lambda3 <- 1+(N-1)*alpha0+(J-1)*N*alpha1

  Var <- sigma2y/N * (I*J*lambda2*lambda3)/
    ( (U^2+I*J*U-J*W-I*V)*lambda3-(U^2-I*V)*lambda2 )

  power <- pt(qt(1-alpha/2, I-2), I-2, ncp=eff/sqrt(Var), lower.tail = F)
  return(power)
}


####### closed-cohort design ######

## Function to calculate predicted number of clusters to get at least (1-beta) power, the actual power that can be reached, and the analytical variance
## under the closed-cohort design, ATE test
ss_CC_ATE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, beta=0.2){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC
  # alpha: type I error rate
  # beta: type II error rate
  #
  # Output:
  # I: predicted number of clusters

  nca <- 1
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

    tau3 <- 1+(N-1)*(alpha0-alpha1)-alpha2
    tau4 <- 1+(N-1)*alpha0+(J-1)*(N-1)*alpha1+(J-1)*alpha2

    Var <- sigma2y/N * (I*J*tau3*tau4)/( (U^2+I*J*U-J*W-I*V)*tau4-(U^2-I*V)*tau3 )

    power <- pt(qt(1-alpha/2, I-2), I-2, ncp=eff/sqrt(Var), lower.tail = F)

  }
  return(I)
}


## Function to calculate the predicted power given the number of clusters
## under the closed-cohort design, ATE test
power_CC_ATE <- function(N, J, sigma2y=1, sigma2x=1, eff, alpha0, alpha1, alpha2, rho0, alpha=0.05, I){
  # Argument:
  # N: number of individuals in each period
  # J: number of periods in each cluster
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC
  # alpha: type I error rate
  # I: predicted number of clusters
  #
  # Output:
  # power: actual predicted power

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

  tau3 <- 1+(N-1)*(alpha0-alpha1)-alpha2
  tau4 <- 1+(N-1)*alpha0+(J-1)*(N-1)*alpha1+(J-1)*alpha2

  Var <- sigma2y/N * (I*J*tau3*tau4)/( (U^2+I*J*U-J*W-I*V)*tau4-(U^2-I*V)*tau3 )

  power <- pt(qt(1-alpha/2, I-2), I-2, ncp=eff/sqrt(Var), lower.tail = F)
  return(power)
}
