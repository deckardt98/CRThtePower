
# CRThtePower

<!-- badges: start -->
<!-- badges: end -->

The `CRThtePower` package implements sample size and power calculation methods for three-level cluster randomized trials (CRTs) and stepped-wedge cluster randomized trials (SW-CRTs). 

For three-level CRTs, the package supports:
- Three types of treatment effect estimands: heterogeneous treatment effect (HTE), average treatment effect (ATE), and unadjusted average treatment effect.
- Three levels of randomization: cluster-level, subcluster-level, and individual-level randomization.

For SW-CRTs, the package supports:
- Two types of treatment effect estimands: heterogeneous treatment effect (HTE) and average treatment effect (ATE).
- Two stepped-wedge designs: cross-sectional and closed-cohort designs.

The following key functions are included:
- `sspower.3level`: Calculates the required number of clusters or the achieved power for three-level CRTs with a univariate effect modifier.
- `sspower.multi.3level`: Extends the above function to handle multivariate effect modifiers for three-level CRTs.
- `sspower.SW`: Calculates the required number of clusters or the achieved power for SW-CRTs with a univariate effect modifier.
- `sspower.multi.SW`: Extends the above function to handle multivariate effect modifiers for SW-CRTs.

The sample size and power calculation methodologies for three-level CRTs are formalized in the paper "Planning Three-Level Cluster Randomized Trials to Assess Treatment Effect Heterogeneity" by Li et al. 
Similarly, the methodologies for SW-CRTs are described in the paper "Planning Stepped Wedge Cluster Randomized Trials to Detect Treatment Effect Heterogeneity" by Li et al.

## Installation

The released version of CRThtePower can be installed from
[CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("devtools") # Run this line if you don't have devtools installed
devtools::install_github("deckardt98/CRThtePower")
```

## Example

This is an example for calculating the number of clusters for a three-level cluster randomized trial with a univariate effect modifier:

``` r
library(CRThtePower)
## basic example code
power.example <- sspower.3level(power = 0.8,
                                nc = 20,
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
                                eff = 0.2,
                                beta3 = 0.2,
                                beta4 = 0.05,
                                estimand = "HTE",
                                random.level = "cluster",
                                verbose = TRUE)
print(power.example)
```


This is an example for calculating the number of clusters for a three-level cluster randomized trial with a multivariate effect modifier:

``` r
library(CRThtePower)
## basic example code
power.example <- sspower.multi.3level(power = 0.8,
                                      nc = 20,
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
                                      verbose = TRUE)
print(power.example)
```



This is an example for calculating the number of clusters for a stepped wedge cluster randomized trial with a univariate effect modifier:

``` r
library(CRThtePower)
## basic example code
power.example <- sspower.SW(power = 0.8,
                            I = 20,
                            alpha = 0.05,
                            sigma2x = 1,
                            sigma2y = 1,
                            rho0 = 0.02,
                            rho1 = 0.02,
                            alpha0 = 0.02,
                            alpha1 = 0.02,
                            alpha2 = 0.02,
                            N = 50,
                            J = 6,
                            eff = 0.1,
                            estimand = "HTE",
                            design = "cross-sectional",
                            verbose = TRUE)
print(power.example)
```


This is an example for calculating the number of clusters for a stepped wedge cluster randomized trial with a multivariate effect modifier:

``` r
library(CRThtePower)
## basic example code
power.example <- sspower.multi.SW(power = 0.8,
                                  I = 20,
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
                                  verbose = TRUE)
print(power.example)
```

