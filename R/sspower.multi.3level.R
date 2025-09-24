#' @title Sample Size and Power Calculations for Three-Level Cluster Randomized Trials with Multivariate Effect Modifiers
#'
#' @description The function `sspower.multi.3level` calculates either the required number of clusters or the achieved power for three-level cluster randomized
#' trials (CRTs) with multivariate effect modifiers (dimension greater than or equal to 2). It accommodates randomization at the cluster, subcluster, or individual level
#' and supports hypothesis testing for the heterogeneous treatment effects (HTEs). Users can specify the number of clusters (level-3) using the `nc` argument.
#' If this argument is provided, the function calculates power for the specified hypothesis test and disregards any input for the power parameter. If the number of
#' clusters is not provided, the function computes the required number of clusters to achieve a specified power threshold, which defaults to 0.8.
#'
#' @usage
#' sspower.multi.3level(power = 0.8,
#'                      nc = NULL,
#'                      alpha = 0.05,
#'                      pi = 0.5,
#'                      sigma2x = diag(rep(1,3)),
#'                      sigma2y = 1,
#'                      rho0 = diag(0.02,3),
#'                      rho1 = diag(0.02,3),
#'                      alpha0 = 0.02,
#'                      alpha1 = 0.02,
#'                      m = 50,
#'                      ns = 6,
#'                      eff = rep(0.5,3),
#'                      random.level = "cluster",
#'                      max.ss = 10^6,
#'                      verbose = TRUE)
#'
#' @param power A numeric value between 0 and 1 specifying the desired power level for sample size estimation. Default is `0.8`.
#' @param nc A numeric value specifying the number of clusters provided by the user to estimate the achievable power. Default is `NULL`.
#' @param alpha A numeric value between 0 and 1 representing the type I error rate. Default is `0.05`.
#' @param pi A rational number between 0 and 1 indicating the proportion of treated units at the randomization level. Default is `0.5`, representing balanced allocation.
#' @param sigma2x A positive-definite covariance matrix for the marginal effect modifiers. Default is `diag(rep(1,3))`.
#' @param sigma2y A positive numeric value for the conditional variance of the continuous outcome. Default is `1`.
#' @param rho0 A cross-correlation matrix with elements between 0 and 1 representing the intraclass correlation coefficients (ICCs) for within-subcluster covariates. Default is `diag(0.02,3)`.
#' @param rho1 A cross-correlation matrix with elements between 0 and 1 representing the intraclass correlation coefficients (ICCs) for between-subcluster covariates. Default is `diag(0.02,3)`.
#' @param alpha0 A numeric value between 0 and 1 representing the ICC for within-subcluster outcome variability. Default is `0.02`.
#' @param alpha1 A numeric value between 0 and 1 representing the ICC for between-subcluster outcome variability. Default is `0.02`.
#' @param m A numeric value greater than 2 specifying the number of participants per subcluster. Default is `50`.
#' @param ns A numeric value greater than 2 specifying the number of subclusters per cluster. Default is `6`.
#' @param eff A nonzero vector of numeric values for the (unstandardized) effect size of the target estimand. The effect size can optionally be standardized relative to `sigma2y`, the
#' conditional variance of the outcome. Default is `rep(0.5,3)`.
#' @param random.level A character string specifying the level of randomization. Supported values are `"cluster"` (randomization at the cluster level),
#' `"subcluster"` (randomization at the subcluster level), and `"individual"` (individual-level randomization). Default is `"cluster"`.
#' @param max.ss A numeric value specifying the maximum number of clusters to be searched. Default is `10^6`.
#' @param verbose A logical value indicating whether parameter reiterations and supplementary messages should be displayed (`TRUE`) or suppressed (`FALSE`). Default is `TRUE`.
#'
#' @details
#' This method calculates the variances of the effects of interest using Generalized Least Squares estimators and large-sample approximations, based on the input parameters.
#' These variances are subsequently applied to construct classic sample size or power formulas, supporting calculations for both objectives.
#' The approach accommodates three levels of randomization—cluster, subcluster, and individual—and primarily focuses on heterogeneous treatment effects (HTE),
#' as outlined in the paper "Planning Three-Level Cluster Randomized Trials to Assess Treatment Effect Heterogeneity" by Li et al.
#' Notably, when the target estimand is the average treatment effect (ATE), the sample size and power calculations are consistent with those derived from `sspower.3level`,
#' as the intraclass correlation coefficients (ICCs) of the covariates do not affect these computations.
#'
#' @return
#' The function `sspower.multi.3level` returns one of two outputs: either an integer specifying the required number of clusters along with the actual power achieved
#' (slightly larger than the specified power), or a decimal value representing the achieved power based on the provided sample size. The power is calculated to four decimal places.
#' Informative messages summarizing key parameter choices and results are displayed by default.
#'
#' @export
#'
#' @examples
#' # Predict the achieved power for 20 clusters with default parameter values.
#' power.example <- sspower.multi.3level(power = 0.8,
#'                                       nc = 20,
#'                                       alpha = 0.05,
#'                                       pi = 0.5,
#'                                       sigma2x = diag(rep(1,3)),
#'                                       sigma2y = 1,
#'                                       rho0 = diag(0.02,3),
#'                                       rho1 = diag(0.02,3),
#'                                       alpha0 = 0.02,
#'                                       alpha1 = 0.02,
#'                                       m = 50,
#'                                       ns = 6,
#'                                       eff = rep(0.5,3),
#'                                       random.level = "cluster",
#'                                       max.ss = 10^6,
#'                                       verbose = TRUE)
#'
#' print(power.example)
#'
#' @importFrom MASS fractions
#' @importFrom stats qnorm pnorm qt pt
#'
sspower.multi.3level <- function(power = 0.8,
                                 nc = NULL,
                                 alpha = 0.05,
                                 pi = 0.5,
                                 sigma2x = diag(rep(1,3)),
                                 sigma2y = 1,
                                 rho0 = diag(0.02,3),
                                 rho1 = diag(0.02,3),
                                 alpha0 = 0.02,
                                 alpha1 = 0.02,
                                 m = 50,
                                 ns = 6,
                                 eff = rep(0.5,3),
                                 random.level = "cluster",
                                 max.ss = 10^6,
                                 verbose = TRUE) {

  ## Error messages
  if (!is.numeric(power) || power <= 0 || power >= 1 || length(power)!=1)
    stop('Target power must be a real number in (0,1).')

  if (!is.null(nc)){
    if (!is.numeric(nc) || nc <= 0 || length(nc)!=1){
      stop('Number of clusters must be a positive integer.')
    }
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha)!=1)
    stop('Type I error rate must be a single number in (0,1)')

  if (!is.numeric(pi) || pi <= 0 || pi >= 1 || length(pi)!=1)
    stop('Proportion of treated units at the randomization level must be a real number in (0,1).')

  if (!all(eigen(sigma2x)$values > 0) || !is.matrix(sigma2x))
    stop('Covariance matrix for the marginal effect modifiers must be a positive-definite matrix.')

  if (!is.numeric(sigma2y) || sigma2y <= 0 || length(sigma2y)!=1)
    stop('Conditional variance of the outcome must be a positive real number.')

  if (!all(rho0 >= 0 & rho0 < 1) || !is.matrix(rho0))
    stop('Cross-correlation matrix for the within-subcluster covariate must be a matrix with all elements within [0,1).')

  if (!all(rho1 >= 0 & rho1 < 1) || !is.matrix(rho1))
    stop('Cross-correlation matrix for the between-subcluster covariate must be a matrix with all elements within [0,1).')

  if (!is.numeric(alpha0) || alpha0 < 0 || alpha0 >= 1 || length(alpha0)!=1)
    stop('Within-subcluster outcome Intracluster Correlation Coefficient must be a real number in [0,1).')

  if (!is.numeric(alpha1) || alpha1 < 0 || alpha1 >= 1 || length(alpha1)!=1)
    stop('Between-subcluster outcome Intracluster Correlation Coefficient must be a real number in [0,1).')

  if (!is.numeric(m) || m <= 0 || length(m)!=1){
    stop('Number of participants per subcluster must be a positive integer.')
  }

  if (!is.numeric(ns) || ns <= 0 || length(ns)!=1){
    stop('Number of subclusters per cluster must be a positive integer.')
  }

  if (!is.numeric(eff) || !any(eff != 0))
    stop('Effect sizes of the target estimand must be nonzero real numbers.')

  if (length(eff)!=ncol(sigma2x) || length(eff)!=ncol(rho0) || length(eff)!=ncol(rho1))
    stop('Dimension of covariates is not compatible.')

  if ( !(random.level %in% c("cluster", "subcluster", "individual")) || length(random.level)!=1)
    stop('Randomization level should be either "cluster" or "subcluster" or "individual".')

  if (!is.numeric(max.ss) || max.ss <= 0 || length(max.ss)!=1)
    stop('Maximum number of clusters to be searched must be a positive real number.')

  if (!is.logical(verbose))
    stop('Message presentation indicator should be a logical argument.')


  ######################################################################################################

  # Re-iterate the given effect sizes and the chosen test
  if (verbose == TRUE) {

    cat('Target treatment effect estimand:\n')
    cat('Heterogeneous treatment effects')
    cat('\n')

    if (random.level == "cluster") {
      cat('Randomization at the cluster level (level-3)')
    } else if (random.level == "subcluster") {
      cat('Randomization at the subcluster level (level-2)')
    } else if (random.level == "individual") {
      cat('Individually randomized trial (level-1)')
    }

    cat(paste0('\nEffect size:\n', paste(eff, collapse = " ")))
  }

  ## Function to estimate power given the number of clusters
  pred.power <- function(nc) {

    if (random.level=="cluster"){
      pred.power <- L3_HTE_power.multi(nc, eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi)
    }

    if (random.level=="subcluster"){
      pred.power <- L2_HTE_power.multi(nc, eff, rho0, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi)
    }

    if (random.level=="individual"){
      pred.power <- L1_HTE_power.multi(nc, eff, alpha0, sigma2x, sigma2y, alpha, m, ns, pi)
    }

    return(pred.power)
  }



  ## Function to estimate number of clusters based on the required power level
  cluster.number <- function(power) {

    if (random.level=="cluster"){
      n.out <- L3_HTE_nc.multi(eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi, max.ss)
    }

    if (random.level=="subcluster"){
      n.out <- L2_HTE_nc.multi(eff, rho0, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi, max.ss)
    }

    if (random.level=="individual"){
      n.out <- L1_HTE_nc.multi(eff, alpha0, sigma2x, sigma2y, power, alpha, m, ns, pi, max.ss)
    }

    return(n.out)
  }

  ## If the user specifies nc, the provided power value will be ignored. The function will calculate the achieved power based solely on the given nc.
  if (!is.null(nc)){
    ans <- as.numeric(pred.power(nc))
    if (verbose==TRUE){
      cat(paste0("\n\nPredicted power for the provided ", nc, " clusters:\n"))
      cat(paste0(round(ans, 4), "\n"))
    }
  } else {
    ## If the user does not specify nc, the formula will give the required number of clusters to hit the provided power.
    ans <- as.numeric(cluster.number(power))
    actual.power <- round(as.numeric(pred.power(ans)), 4)
    if (verbose==TRUE){
      cat(paste0("\n\nRequired number of clusters to achieve an actual power of ", actual.power, ":\n"))
      cat(paste0(ans, "\n"))
    }
  }

  return(ans)
}

