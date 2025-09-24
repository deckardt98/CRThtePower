#' @title Sample Size and Power Calculations for Three-Level Cluster Randomized Trials with A Univariate Effect Modifier
#'
#' @description The function `sspower.3level` calculates either the required number of clusters or the achieved power for three-level cluster randomized
#' trials (CRTs) with a univariate effect modifier. It accommodates randomization at the cluster, subcluster, or individual level and supports hypothesis testing
#' for three types of effects: the heterogeneous treatment effect (HTE), the average treatment effect (ATE), and the unadjusted ATE (estimated via a linear mixed model
#' without adjusting for the effect modifier). In particular, the sample size calculations and power analysis for testing the ATE or unadjusted ATE are independent of
#' the effect modifier, meaning that the dimension of the effect modifiers does not affect the results. Users can specify the number of clusters (level-3) using the
#' `nc` argument. If this argument is provided, the function calculates power for the specified hypothesis test and disregards any input for the power parameter.
#' If the number of clusters is not provided, the function computes the required number of clusters to achieve a specified power threshold, which defaults to 0.8.
#'
#' @usage
#' sspower.3level(power = 0.8,
#'                nc = NULL,
#'                alpha = 0.05,
#'                pi = 0.5,
#'                sigma2x = 1,
#'                sigma2y = 1,
#'                rho0 = 0.02,
#'                rho1 = 0.02,
#'                alpha0 = 0.02,
#'                alpha1 = 0.02,
#'                m = 50,
#'                ns = 6,
#'                eff = 0.2,
#'                beta3 = 0.2,
#'                beta4 = 0.05,
#'                estimand = "HTE",
#'                random.level = "cluster",
#'                verbose = TRUE)
#'
#' @param power A numeric value between 0 and 1 specifying the desired power level for sample size estimation. Default is `0.8`.
#' @param nc A numeric value specifying the number of clusters provided by the user to estimate the achievable power. Default is `NULL`.
#' @param alpha A numeric value between 0 and 1 representing the type I error rate. Default is `0.05`.
#' @param pi A rational number between 0 and 1 indicating the proportion of treated units at the randomization level. Default is `0.5`, representing balanced allocation.
#' @param sigma2x A positive numeric value for the marginal variance of the univariate effect modifier. Default is `1`.
#' @param sigma2y A positive numeric value for the conditional variance of the continuous outcome. Default is `1`.
#' @param rho0 A numeric value between 0 and 1 representing the intraclass correlation coefficient (ICC) for within-subcluster covariate variability. Default is `0.02`.
#' @param rho1 A numeric value between 0 and 1 representing the ICC for between-subcluster covariate variability. Default is `0.02`.
#' @param alpha0 A numeric value between 0 and 1 representing the ICC for within-subcluster outcome variability. Default is `0.02`.
#' @param alpha1 A numeric value between 0 and 1 representing the ICC for between-subcluster outcome variability. Default is `0.02`.
#' @param m A numeric value greater than 2 specifying the number of participants per subcluster. Default is `50`.
#' @param ns A numeric value greater than 2 specifying the number of subclusters per cluster. Default is `6`.
#' @param eff A nonzero numeric value for the (unstandardized) effect size of the target estimand. The effect size can optionally be standardized relative to `sigma2y`, the
#' conditional variance of the outcome. Default is `0.1`.
#' @param beta3 A nonzero numeric value specifying the true covariate effect, corresponding to the model coefficient of the covariate.
#' Default is `0.2`. This parameter is required only when the estimand is `"unATE"`.
#' @param beta4 A nonzero numeric value specifying the true effect of treatment effect heterogeneity, corresponding to the model coefficient
#' of the covariate-treatment interaction. Default is `0.05`. This parameter is required only when the estimand is `"unATE"`.
#' @param estimand A character string indicating the type of treatment effect estimand. Supported values are `"HTE"` (heterogeneous treatment effect),
#' `"ATE"` (average treatment effect), and `"unATE"` (unadjusted average treatment effect). Default is `"HTE"`.
#' @param random.level A character string specifying the level of randomization. Supported values are `"cluster"` (randomization at the cluster level),
#' `"subcluster"` (randomization at the subcluster level), and `"individual"` (individual-level randomization). Default is `"cluster"`.
#' @param verbose A logical value indicating whether parameter reiterations and supplementary messages should be displayed (`TRUE`) or suppressed (`FALSE`). Default is `TRUE`.
#'
#' @details
#' Based on the input parameters, the method first calculates the variances of the effects of interest using Generalized Least Squares estimators and large-sample approximations.
#' These variances are then applied to construct either classic sample size formulas or power formulas, facilitating both sample size and power calculations. Specifically, normal quantiles
#' are employed for sample size calculations, except in the case of cluster-level randomization with ATE and unATE estimands, where t-quantiles are used to incorporate small-sample adjustments.
#' The methods support three estimands (HTE, ATE, and unATE) and three randomization levels (cluster, subcluster, and individual), as detailed in the paper entitled "Planning Three-Level Cluster
#' Randomized Trials to Assess Treatment Effect Heterogeneity" by Li et al.
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
#' power.example <- sspower.3level(power = 0.8,
#'                                 nc = 20,
#'                                 alpha = 0.05,
#'                                 pi = 0.5,
#'                                 sigma2x = 1,
#'                                 sigma2y = 1,
#'                                 rho0 = 0.02,
#'                                 rho1 = 0.02,
#'                                 alpha0 = 0.02,
#'                                 alpha1 = 0.02,
#'                                 m = 50,
#'                                 ns = 6,
#'                                 eff = 0.2,
#'                                 beta3 = 0.2,
#'                                 beta4 = 0.05,
#'                                 estimand = "HTE",
#'                                 random.level = "cluster",
#'                                 verbose = TRUE)
#'
#' print(power.example)
#'
#' @importFrom MASS fractions
#' @importFrom stats qnorm pnorm qt pt
#'
sspower.3level <- function(power = 0.8,
                           nc = NULL,
                           alpha = 0.05,
                           pi = 0.5,
                           sigma2x = 1,
                           sigma2y = 1,
                           rho0 = 0.02,
                           rho1 = 0.02,
                           alpha0 = 0.02,
                           alpha1 = 0.02,
                           m = 50,
                           ns = 6,
                           eff = 0.1,
                           beta3 = 0.2,
                           beta4 = 0.05,
                           estimand = "HTE",
                           random.level = "cluster",
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

  if (!is.numeric(sigma2x) || sigma2x <= 0 || length(sigma2x)!=1)
    stop('Marginal variance of the univariate effect modifier must be a positive real number.')

  if (!is.numeric(sigma2y) || sigma2y <= 0 || length(sigma2y)!=1)
    stop('Conditional variance of the outcome must be a positive real number.')

  if (!is.numeric(rho0) || rho0 < 0 || rho0 >= 1 || length(rho0)!=1)
    stop('Within-subcluster covariate Intracluster Correlation Coefficient must be a real number in [0,1).')

  if (!is.numeric(rho1) || rho1 < 0 || rho1 >= 1 || length(rho1)!=1)
    stop('Between-subcluster covariate Intracluster Correlation Coefficient must be a real number in [0,1).')

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

  if (!is.numeric(eff) || eff == 0 || length(eff)!=1)
    stop('Effect size of the target estimand must be a nonzero real number.')

  if ( !(estimand %in% c("HTE", "ATE", "unATE")) || length(estimand)!=1)
    stop('Treatment effect estimand should be either "HTE" or "ATE" or "unATE".')

  if ( !(random.level %in% c("cluster", "subcluster", "individual")) || length(random.level)!=1)
    stop('Randomization level should be either "cluster" or "subcluster" or "individual".')

  if (!is.logical(verbose))
    stop('Message presentation indicator should be a logical argument.')

  if (estimand == "unATE"){
    if (!is.numeric(beta3) || length(beta3)!=1)
      stop('True effect size of the covariate must be a real number.')
    if (!is.numeric(beta4) || length(beta4)!=1)
      stop('True effect size of the treatment-covariate interaction must be a real number.')
  }

  ######################################################################################################

  # Re-iterate the given effect sizes and the chosen test
  if (verbose == TRUE) {

    cat('Target treatment effect estimand:\n')
    if (estimand == "HTE") {
      cat('Heterogeneous treatment effect')
    } else if (estimand == "ATE") {
      cat('Average treatment effect')
    } else if (estimand == "unATE") {
      cat('Unadjusted average treatment effect')
    }

    cat('\n')

    if (random.level == "cluster") {
      cat('Randomization at the cluster level (level-3)')
    } else if (random.level == "subcluster") {
      cat('Randomization at the subcluster level (level-2)')
    } else if (random.level == "individual") {
      cat('Individually randomized trial (level-1)')
    }

    cat(paste0('\nEffect size:\n', eff))
  }

  ## Effect size might be negative
  eff <- abs(eff)

  ## Function to estimate power given the number of clusters
  pred.power <- function(nc) {

    if (estimand=="HTE"){

      if (random.level=="cluster"){
        pred.power <- L3_HTE_power(nc, eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi)
      }

      if (random.level=="subcluster"){
        pred.power <- L2_HTE_power(nc, eff, rho0, alpha0, alpha1, sigma2x, sigma2y, alpha, m, ns, pi)
      }

      if (random.level=="individual"){
        pred.power <- L1_HTE_power(nc, eff, alpha0, sigma2x, sigma2y, alpha, m, ns, pi)
      }

    } else if (estimand=="ATE"){

      if (random.level=="cluster"){
        pred.power <- L3_ATE_power(nc, eff, alpha0, alpha1, sigma2y, alpha, m, ns, pi)
      }

      if (random.level=="subcluster"){
        pred.power <- L2_ATE_power(nc, eff, alpha0, alpha1, sigma2y, alpha, m, ns, pi)
      }

      if (random.level=="individual"){
        pred.power <- L1_ATE_power(nc, eff, alpha0, sigma2y, alpha, m, ns, pi)
      }

    } else if (estimand=="unATE"){

      if (random.level=="cluster"){
        pred.power <- L3_ATE2_power(nc, beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, alpha, m, ns, pi)
      }

      if (random.level=="subcluster"){
        pred.power <- L2_ATE2_power(nc, beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, alpha, m, ns, pi)
      }

      if (random.level=="individual"){
        pred.power <- L1_ATE2_power(nc, beta3, beta4, eff, alpha0, rho0, sigma2y, sigma2x, alpha, m, ns, pi)
      }

    }
    return(pred.power)
  }



  ## Function to estimate number of clusters based on the required power level
  cluster.number <- function(power) {

    if (estimand=="HTE"){

      if (random.level=="cluster"){
        n.out <- L3_HTE_nc(eff, rho0, rho1, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi)
      }

      if (random.level=="subcluster"){
        n.out <- L2_HTE_nc(eff, rho0, alpha0, alpha1, sigma2x, sigma2y, power, alpha, m, ns, pi)
      }

      if (random.level=="individual"){
        n.out <- L1_HTE_nc(eff, alpha0, sigma2x, sigma2y, power, alpha, m, ns, pi)
      }

    } else if (estimand=="ATE"){

      if (random.level=="cluster"){
        n.out <- L3_ATE_nc(eff, alpha0, alpha1, sigma2y, power, alpha, m, ns, pi)
      }

      if (random.level=="subcluster"){
        n.out <- L2_ATE_nc(eff, alpha0, alpha1, sigma2y, power, alpha, m, ns, pi)
      }

      if (random.level=="individual"){
        n.out <- L1_ATE_nc(eff, alpha0, sigma2y, power, alpha, m, ns, pi)
      }

    } else if (estimand=="unATE"){

      if (random.level=="cluster"){
        n.out <- L3_ATE2_nc(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, power, alpha, m, ns, pi)
      }

      if (random.level=="subcluster"){
        n.out <- L2_ATE2_nc(beta3, beta4, eff, alpha0, alpha1, rho0, rho1, sigma2y, sigma2x, power, alpha, m, ns, pi)
      }

      if (random.level=="individual"){
        n.out <- L1_ATE2_nc(beta3, beta4, eff, alpha0, rho0, sigma2y, sigma2x, power, alpha, m, ns, pi)
      }

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

