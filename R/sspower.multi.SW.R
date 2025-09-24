#' @title Sample Size and Power Calculations for Stepped Wedge Cluster Randomized Trials with A Multivariate Effect Modifier
#'
#' @description
#' The `sspower.multi.SW` function calculates sample size or power for stepped wedge cluster randomized trials (SW-CRTs) with a multivariate effect modifier.
#' It supports both cross-sectional and closed-cohort stepped wedge designs and focuses on hypothesis testing for heterogeneous treatment effects (HTE).
#' Users can specify the number of clusters using the `I` argument. When `I` is provided, the function computes the power for the specified hypothesis
#' test, disregarding the `power` parameter. If `I` is not specified, the function calculates the required number of clusters to achieve the desired power
#' threshold, which defaults to 0.8.
#'
#' @usage
#' sspower.multi.SW(power = 0.8,
#'                  I = NULL,
#'                  alpha = 0.05,
#'                  sigma2x = diag(rep(1,3)),
#'                  sigma2y = 1,
#'                  rho0 = diag(rep(0.02,3)),
#'                  rho1 = diag(rep(0.02,3)),
#'                  alpha0 = 0.02,
#'                  alpha1 = 0.02,
#'                  alpha2 = 0.02,
#'                  N = 50,
#'                  J = 6,
#'                  eff = rep(0.5,3),
#'                  design = "cross-sectional",
#'                  verbose = TRUE)
#'
#' @param power A numeric value between 0 and 1 specifying the desired power level for sample size estimation. Default is `0.8`.
#' @param I A numeric value specifying the number of clusters provided by the user to estimate the achievable power. Default is `NULL`.
#' @param alpha A numeric value between 0 and 1 representing the type I error rate. Default is `0.05`.
#' @param sigma2x A positive-definite covariance matrix for the marginal effect modifiers. Default is `diag(rep(1,3))`.
#' @param sigma2y A positive numeric value for the conditional variance of the continuous outcome. Default is `1`.
#' @param rho0 A cross-correlation matrix with elements between 0 and 1 representing the intraclass correlation coefficients (ICCs) for within-period covariates. Default is `diag(0.02,3)`.
#' @param rho1 A cross-correlation matrix with elements between 0 and 1 representing the intraclass correlation coefficients (ICCs) for between-period covariates. Default is `diag(0.02,3)`.
#' @param alpha0 A numeric value between 0 and 1 representing the ICC for within-period outcome variability. Default is `0.02`.
#' @param alpha1 A numeric value between 0 and 1 representing the ICC for between-period outcome variability. Default is `0.02`.
#' @param alpha2 A numeric value between 0 and 1 representing the ICC for within-individual outcome variability, which is only required for closed-cohort design. Default is `0.02`.
#' @param N A numeric value greater than 2 specifying the number of individuals per period. Default is `50`.
#' @param J A numeric value greater than 2 specifying the number of periods per cluster. Default is `6`.
#' @param eff A nonzero vector of numeric values for the (unstandardized) effect size of the target estimand. The effect size can optionally be standardized relative to `sigma2y`, the
#' conditional variance of the outcome. Default is `rep(0.5,3)`.
#' @param design A character string indicating the type of stepped wedge design. Supported values are `"cross-sectional"` (cross-sectional design) and
#' `"closed-cohort"` (closed-cohort design). Default is `"cross-sectional"`.
#' @param verbose A logical value indicating whether parameter reiterations and supplementary messages should be displayed (`TRUE`) or suppressed (`FALSE`). Default is `TRUE`.
#'
#' @details
#' This method employs Generalized Least Squares estimators and large-sample approximations to compute the variances of the effects of interest.
#' These variances are then utilized to derive formulas for sample size and power calculations, supporting analysis for both objectives.
#' The approach accommodates two stepped wedge designs—cross-sectional and closed-cohort—with a primary focus on heterogeneous treatment effects (HTE).
#' As described in the paper "Planning Stepped Wedge Cluster Randomized Trials to Detect Treatment Effect Heterogeneity" by Li et al.,
#' the sample size and power calculations align with those produced by `sspower.SW` when the target estimand is the average treatment effect,
#' as the intraclass correlation coefficients (ICCs) of the covariates do not affect these computations.
#'
#' @return
#' The function `sspower.multi.SW` returns one of two outputs: either an integer specifying the required number of clusters along with the actual power achieved
#' (slightly larger than the specified power), or a decimal value representing the achieved power based on the provided sample size. The power is calculated to four decimal places.
#' Informative messages summarizing key parameter choices and results are displayed by default.
#'
#' @export
#'
#' @examples
#' # Predict the achieved power for 20 clusters with default parameter values.
#' power.example <- sspower.multi.SW(power = 0.8,
#'                                   I = 20,
#'                                   alpha = 0.05,
#'                                   sigma2x = diag(rep(1,3)),
#'                                   sigma2y = 1,
#'                                   rho0 = diag(0.02,3),
#'                                   rho1 = diag(0.02,3),
#'                                   alpha0 = 0.02,
#'                                   alpha1 = 0.02,
#'                                   alpha2 = 0.02,
#'                                   N = 50,
#'                                   J = 6,
#'                                   eff = rep(0.5,3),
#'                                   design = "cross-sectional",
#'                                   verbose = TRUE)
#'
#' print(power.example)
#'
#' @importFrom stats qnorm pnorm qt pt
#'
sspower.multi.SW <- function(power = 0.8,
                             I = NULL,
                             alpha = 0.05,
                             sigma2x = diag(rep(1,3)),
                             sigma2y = 1,
                             rho0 = diag(0.02,3),
                             rho1 = diag(0.02,3),
                             alpha0 = 0.02,
                             alpha1 = 0.02,
                             alpha2 = 0.02,
                             N = 50,
                             J = 6,
                             eff = rep(0.5,3),
                             design = "cross-sectional",
                             verbose = TRUE) {

  ## Error messages
  if (!is.numeric(power) || power <= 0 || power >= 1 || length(power)!=1)
    stop('Target power must be a real number in (0,1).')

  if (!is.null(I)){
    if (!is.numeric(I) || I <= 0 || length(I)!=1){
      stop('Number of clusters must be a positive integer.')
    }
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha)!=1)
    stop('Type I error rate must be a single number in (0,1)')

  if (!all(eigen(sigma2x)$values > 0) || !is.matrix(sigma2x))
    stop('Covariance matrix for the marginal effect modifiers must be a positive-definite matrix.')

  if (!is.numeric(sigma2y) || sigma2y <= 0 || length(sigma2y)!=1)
    stop('Conditional variance of the outcome must be a positive real number.')

  if (!all(rho0 >= 0 & rho0 < 1) || !is.matrix(rho0))
    stop('Cross-correlation matrix for the within-period covariate must be a matrix with all elements within [0,1).')

  if (!all(rho1 >= 0 & rho1 < 1) || !is.matrix(rho1))
    stop('Cross-correlation matrix for the between-period covariate must be a matrix with all elements within [0,1).')

  if (!is.numeric(alpha0) || alpha0 < 0 || alpha0 >= 1 || length(alpha0)!=1)
    stop('Within-period outcome Intracluster Correlation Coefficient must be a real number in [0,1).')

  if (!is.numeric(alpha1) || alpha1 < 0 || alpha1 >= 1 || length(alpha1)!=1)
    stop('Between-period outcome Intracluster Correlation Coefficient must be a real number in [0,1).')

  if (!is.numeric(alpha2) || alpha2 < 0 || alpha2 >= 1 || length(alpha2)!=1)
    stop('Within-individual outcome Intracluster Correlation Coefficient must be a real number in [0,1).')

  if (!is.numeric(N) || N <= 0 || length(N)!=1){
    stop('Number of individuals per period must be a positive integer.')
  }

  if (!is.numeric(J) || J <= 0 || length(J)!=1){
    stop('Number of periods per cluster must be a positive integer.')
  }

  if (!is.numeric(eff) || !any(eff != 0))
    stop('Effect sizes of the heterogeneous treatment effect must be nonzero real numbers.')

  if (length(eff)!=ncol(sigma2x) || length(eff)!=ncol(rho0) || length(eff)!=ncol(rho1))
    stop('Dimension of covariates is not compatible.')

  if ( !(design %in% c("cross-sectional", "closed-cohort")) || length(design)!=1)
    stop('Stepped wedge design should be either "cross-sectional" or "closed-cohort".')

  if (!is.logical(verbose))
    stop('Message presentation indicator should be a logical argument.')

  if (!is.null(I)){
    if (I %% (J - 1) != 0){
      stop('Number of clusters (I) is not divisible by number of periods minus one (J-1).')
    }
  }

  ######################################################################################################

  ## Effect size might be negative
  eff <- abs(eff)

  # Re-iterate the given effect sizes and the chosen test
  if (verbose == TRUE) {

    cat('Target treatment effect estimand:\n')
    cat('Heterogeneous treatment effect')
    cat('\n')

    if (design == "cross-sectional") {
      cat('Cross-sectional stepped wedge design')
    } else if (design == "closed-cohort") {
      cat('Closed-cohort stepped wedge design')
    }

    cat(paste0('\nEffect size:\n', paste(eff, collapse = " ")))
  }

  ## Function to estimate power given the number of clusters
  pred.power <- function(I) {

    if (design=="cross-sectional"){
      pred.power <- power_CS_HTE.multi(N, J, sigma2y, sigma2x, eff, alpha0, alpha1, rho0, rho1, alpha, I)
    }

    if (design=="closed-cohort"){
      pred.power <- power_CC_HTE.multi(N, J, sigma2y, sigma2x, eff, alpha0, alpha1, alpha2, rho0, alpha, I)
    }
    return(pred.power)
  }



  ## Function to estimate number of clusters based on the required power level
  cluster.number <- function(power) {

    if (design=="cross-sectional"){
      n.out <- ss_CS_HTE.multi(N, J, sigma2y, sigma2x, eff, alpha0, alpha1, rho0, rho1, alpha, 1-power)
    }

    if (design=="closed-cohort"){
      n.out <- ss_CC_HTE.multi(N, J, sigma2y, sigma2x, eff, alpha0, alpha1, alpha2, rho0, alpha, 1-power)
    }
    return(n.out)
  }

  ## If the user specifies I, the provided power value will be ignored. The function will calculate the achieved power based solely on the given I.
  if (!is.null(I)){
    ans <- as.numeric(pred.power(I))
    if (verbose==TRUE){
      cat(paste0("\n\nPredicted power for the provided ", I, " clusters:\n"))
      cat(paste0(round(ans, 4), "\n"))
    }
  } else {
    ## If the user does not specify I, the formula will give the required number of clusters to hit the provided power.
    ans <- as.numeric(cluster.number(power))
    actual.power <- round(as.numeric(pred.power(ans)), 4)
    if (verbose==TRUE){
      cat(paste0("\n\nRequired number of clusters to achieve an actual power of ", actual.power, ":\n"))
      cat(paste0(ans, "\n"))
    }
  }

  return(ans)
}

