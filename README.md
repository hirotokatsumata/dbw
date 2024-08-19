
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dbw: Doubly Robust Distribution Balancing Weighting Estimation

<!-- badges: start -->
<!-- badges: end -->

**dbw** implements the doubly robust distribution balancing weighting
proposed by Katsumata (2024) <doi:10.1017/psrm.2024.23>, which improves
the augmented inverse probability weighting (AIPW) by estimating
propensity scores with estimating equations suitable for the
pre-specified parameter of interest (e.g., the average treatment effects
or the average treatment effects on the treated) and estimating outcome
models with the estimated inverse probability weights.

It also implements the covariate balancing propensity score proposed by
Imai and Ratkovic (2014) <doi:10.1111/rssb.12027> and the entropy
balancing weighting proposed by Hainmueller (2012)
<doi:10.1093/pan/mpr025>, which use covariate balancing conditions in
propensity score estimation.

The point estimate of the parameter of interest and its uncertainty as
well as coefficients for propensity score estimation and outcome
regression are produced using the M-estimation. The same functions can
be used to estimate average outcomes in missing outcome cases.

## How to Cite

<font size="4"> Katsumata, Hiroto. 2024. “How Should We Estimate Inverse
Probability Weights with Possibly Misspecified Propensity Score Models?”
*Political Science Research and Methods*. </font>

## Installation

<!-- 
You can install the release version of **dbw** from [CRAN](https://CRAN.R-project.org/dbw):
``` r
install.packages("dbw")
```
-->

You can install the development version of **dbw** from
[GitHub](https://github.com/hirotokatsumata/dbw):

``` r
devtools::install_github("hirotokatsumata/dbw")
```

## Example

Simulation from Kang and Shafer (2007) and Imai and Ratkovic (2014)

### Make a toy data set

``` r
library(dbw)
# ATE estimation
# Make a toy data set
# True ATE is 10
tau <- 10
set.seed(12345)
n <- 1000
X <- matrix(stats::rnorm(n * 4, mean = 0, sd = 1), nrow = n, ncol = 4)
prop <- 1 / (1 + exp(X[, 1] - 0.5 * X[, 2] + 0.25 * X[, 3] + 0.1 * X[, 4]))
treat <- rbinom(n, 1, prop)
y <- 210 + 27.4 * X[, 1] + 13.7 * X[, 2] + 13.7 * X[, 3] + 13.7 * X[, 4] + 
     tau * treat + stats::rnorm(n = n, mean = 0, sd = 1)
ybinom <- (y > 210) + 0
df0 <- data.frame(X, treat, y, ybinom)
colnames(df0) <- c("x1", "x2", "x3", "x4", "treat", "y", "ybinom")

# Variables for a misspecified model
Xmis <- data.frame(x1mis = exp(X[, 1] / 2), 
                   x2mis = X[, 2] * (1 + exp(X[, 1]))^(-1) + 10,
                   x3mis = (X[, 1] * X[, 3] / 25 + 0.6)^3, 
                   x4mis = (X[, 2] + X[, 4] + 20)^2)

# Data frame and formulas for propensity score estimation
df <- data.frame(df0, Xmis)
formula_ps_c <- stats::as.formula(treat ~ x1 + x2 + x3 + x4)
formula_ps_m <- stats::as.formula(treat ~ x1mis + x2mis + x3mis + x4mis)

# Formula for a misspecified outcome model
formula_y <- stats::as.formula(y ~ x1mis + x2mis + x3mis + x4mis)
```

### Correct propensity score model

``` r
# Distribution balancing weighting with normalization and without regularization
fitdbwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c, 
               estimand = "ATE", method = "dbw",
               method_y = "wls", data = df, normalize = TRUE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
fitdbwc
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "dbw", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ATE estimate:  10.03
#> 
#> Estimate of E[Y(1)]:  220.6
#> Estimate of E[Y(0)]:  210.5
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#> ps_t_(Intercept)          ps_t_x1          ps_t_x2          ps_t_x3 
#>         -0.06034         -1.17788          0.63462         -0.25104 
#>          ps_t_x4 
#>          0.02505 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#> ps_c_(Intercept)          ps_c_x1          ps_c_x2          ps_c_x3 
#>         -0.06063         -1.01973          0.35814         -0.37044 
#>          ps_c_x4 
#>          0.24305 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#> outcome_t_(Intercept)       outcome_t_x1mis       outcome_t_x2mis 
#>                 25.80                 42.32                  2.56 
#>       outcome_t_x3mis       outcome_t_x4mis 
#>                -56.27                  0.33 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#> outcome_c_(Intercept)       outcome_c_x1mis       outcome_c_x2mis 
#>               47.6875               41.0927               -7.1336 
#>       outcome_c_x3mis       outcome_c_x4mis 
#>              132.9864                0.3921 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 279.66 
#>   For estimating E[Y(0)]:   366.01
summary(fitdbwc)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "dbw", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    10.03     0.6139   16.34        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    220.6      1.249   176.6        0 ***
#> mu_c    210.5      1.185   177.6        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -0.06034    0.09663 -0.6245 5.323e-01    
#> ps_t_x1          -1.17788    0.19477 -6.0476 1.470e-09 ***
#> ps_t_x2           0.63462    0.13761  4.6119 3.991e-06 ***
#> ps_t_x3          -0.25104    0.14051 -1.7867 7.399e-02   .
#> ps_t_x4           0.02505    0.12348  0.2029 8.392e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -0.06063    0.07894  -0.768 4.425e-01    
#> ps_c_x1          -1.01973    0.16801  -6.069 1.283e-09 ***
#> ps_c_x2           0.35814    0.10727   3.339 8.416e-04 ***
#> ps_c_x3          -0.37044    0.15321  -2.418 1.561e-02   *
#> ps_c_x4           0.24305    0.11844   2.052 4.016e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)    25.80   12.78815   2.018   0.0436   *
#> outcome_t_x1mis          42.32    1.46069  28.973   0.0000 ***
#> outcome_t_x2mis           2.56    1.67872   1.525   0.1273    
#> outcome_t_x3mis         -56.27   36.56608  -1.539   0.1239    
#> outcome_t_x4mis           0.33    0.02152  15.336   0.0000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  47.6875   17.63955   2.703 0.0068626  **
#> outcome_c_x1mis        41.0927    1.51316  27.157 0.0000000 ***
#> outcome_c_x2mis        -7.1336    1.95891  -3.642 0.0002710 ***
#> outcome_c_x3mis       132.9864   38.94135   3.415 0.0006377 ***
#> outcome_c_x4mis         0.3921    0.01616  24.256 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 279.66 
#>   For estimating E[Y(0)]:   366.01

# Covariate balancing weighting function without regularization
fitcbwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c, 
               estimand = "ATE", method = "cb",
               method_y = "wls", data = df, normalize = TRUE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitcbwc)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "cb", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    9.535     0.6666    14.3        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    220.4      1.232   178.9        0 ***
#> mu_c    210.9      1.171   180.1        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -0.06813    0.07512 -0.9069 3.644e-01    
#> ps_t_x1          -1.18187    0.13091 -9.0280 0.000e+00 ***
#> ps_t_x2           0.59690    0.09674  6.1701 6.826e-10 ***
#> ps_t_x3          -0.22633    0.10196 -2.2198 2.643e-02   *
#> ps_t_x4           0.03164    0.09635  0.3284 7.426e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -0.05808    0.07251 -0.8011 4.231e-01    
#> ps_c_x1          -1.03898    0.11108 -9.3534 0.000e+00 ***
#> ps_c_x2           0.42285    0.08744  4.8356 1.327e-06 ***
#> ps_c_x3          -0.27476    0.09883 -2.7802 5.433e-03  **
#> ps_c_x4           0.16160    0.09049  1.7860 7.411e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)  24.8545   13.01462   1.910  0.05617   .
#> outcome_t_x1mis        42.2381    1.56581  26.975  0.00000 ***
#> outcome_t_x2mis         2.6286    1.66301   1.581  0.11396    
#> outcome_t_x3mis       -55.1762   34.69947  -1.590  0.11181    
#> outcome_t_x4mis         0.3299    0.02138  15.428  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  50.7071   18.58438   2.728 0.0063628  **
#> outcome_c_x1mis        41.1159    1.59313  25.808 0.0000000 ***
#> outcome_c_x2mis        -7.5699    2.05674  -3.681 0.0002327 ***
#> outcome_c_x3mis       134.7874   38.94892   3.461 0.0005389 ***
#> outcome_c_x4mis         0.3953    0.01657  23.850 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 281.62 
#>   For estimating E[Y(0)]:   350.35

# Entropy balancing weighting function without regularization
fitebwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c, 
               estimand = "ATE", method = "eb",
               method_y = "wls", data = df, normalize = TRUE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitebwc)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "eb", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    9.532     0.7135   13.36        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    221.6      1.253   176.9        0 ***
#> mu_c    212.0      1.210   175.3        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error   z value  Pr(>|z|)    
#> ps_t_(Intercept) -0.888961    0.04269 -20.82186 0.000e+00 ***
#> ps_t_x1          -0.631701    0.06475  -9.75565 0.000e+00 ***
#> ps_t_x2           0.303084    0.04727   6.41118 1.444e-10 ***
#> ps_t_x3          -0.103030    0.04613  -2.23331 2.553e-02   *
#> ps_t_x4           0.003287    0.04676   0.07029 9.440e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept)  0.80649    0.04110  19.621 0.000e+00 ***
#> ps_c_x1          -0.52158    0.05492  -9.496 0.000e+00 ***
#> ps_c_x2           0.22086    0.04145   5.329 9.889e-08 ***
#> ps_c_x3          -0.11303    0.04068  -2.779 5.460e-03  **
#> ps_c_x4           0.05751    0.04026   1.428 1.532e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)   8.6513   13.09000  0.6609  0.50867    
#> outcome_t_x1mis        41.6189    1.73564 23.9790  0.00000 ***
#> outcome_t_x2mis         4.0877    1.61654  2.5287  0.01145   *
#> outcome_t_x3mis       -30.8163   33.57743 -0.9178  0.35874    
#> outcome_t_x4mis         0.3255    0.01982 16.4236  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  29.4010   18.23510   1.612 0.1068903    
#> outcome_c_x1mis        42.3951    1.51878  27.914 0.0000000 ***
#> outcome_c_x2mis        -4.9009    2.00448  -2.445 0.0144854   *
#> outcome_c_x3mis       132.6096   39.45882   3.361 0.0007774 ***
#> outcome_c_x4mis         0.3824    0.01617  23.650 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 330.32 
#>   For estimating E[Y(0)]:   395.69

# Standard logistic regression
fitmlec <- dbw(formula_y = formula_y, formula_ps = formula_ps_c, 
               estimand = "ATE", method = "mle",
               method_y = "wls", data = df, normalize = FALSE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlec)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, normalize = FALSE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    9.586     0.9696   9.886        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    220.2      1.273   173.0        0 ***
#> mu_c    210.6      1.318   159.8        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                  Estimate Std. Error  z value  Pr(>|z|)    
#> ps_t_(Intercept) -0.06302    0.07284  -0.8652 3.869e-01    
#> ps_t_x1          -1.19944    0.10086 -11.8924 0.000e+00 ***
#> ps_t_x2           0.53437    0.07648   6.9869 2.810e-12 ***
#> ps_t_x3          -0.19432    0.07222  -2.6907 7.129e-03  **
#> ps_t_x4           0.04193    0.07666   0.5469 5.844e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error  z value  Pr(>|z|)    
#> ps_c_(Intercept) -0.06302    0.07284  -0.8652 3.869e-01    
#> ps_c_x1          -1.19944    0.10086 -11.8924 0.000e+00 ***
#> ps_c_x2           0.53437    0.07648   6.9869 2.810e-12 ***
#> ps_c_x3          -0.19432    0.07222  -2.6907 7.129e-03  **
#> ps_c_x4           0.04193    0.07666   0.5469 5.844e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)  23.2629   13.22516   1.759  0.07858   .
#> outcome_t_x1mis        42.1074    1.94840  21.611  0.00000 ***
#> outcome_t_x2mis         2.7300    1.67016   1.635  0.10213    
#> outcome_t_x3mis       -52.7308   31.42338  -1.678  0.09333   .
#> outcome_t_x4mis         0.3299    0.02145  15.378  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  59.3617   21.09724   2.814 0.0048972  **
#> outcome_c_x1mis        42.3795    2.01014  21.083 0.0000000 ***
#> outcome_c_x2mis        -8.6234    2.26622  -3.805 0.0001417 ***
#> outcome_c_x3mis       122.0199   42.75283   2.854 0.0043162  **
#> outcome_c_x4mis         0.4025    0.01816  22.166 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 282.08 
#>   For estimating E[Y(0)]:   263.75

# Distribution balancing weighting without normalization and without regularization
fitdbwcnn <- dbw(formula_y = formula_y, formula_ps = formula_ps_c, 
                 estimand = "ATE", method = "dbw",
                 method_y = "wls", data = df, normalize = FALSE, 
                 vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> Warning in dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand =
#> "ATE", : estimated weights are not normalized
summary(fitdbwcnn)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "dbw", method_y = "wls", data = df, normalize = FALSE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    9.543     0.5952   16.03        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    220.4      1.247   176.8        0 ***
#> mu_c    210.9      1.186   177.9        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -0.07246     0.1024 -0.7076 4.792e-01    
#> ps_t_x1          -1.18369     0.1719 -6.8844 5.804e-12 ***
#> ps_t_x2           0.59795     0.1245  4.8021 1.570e-06 ***
#> ps_t_x3          -0.22692     0.1363 -1.6654 9.583e-02   .
#> ps_t_x4           0.03168     0.1224  0.2588 7.958e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -0.05808    0.09742 -0.5962 5.510e-01    
#> ps_c_x1          -1.03898    0.15141 -6.8619 6.797e-12 ***
#> ps_c_x2           0.42285    0.11229  3.7657 1.661e-04 ***
#> ps_c_x3          -0.27476    0.15474 -1.7756 7.579e-02   .
#> ps_c_x4           0.16160    0.12143  1.3309 1.832e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)   24.812   12.85168   1.931  0.05353   .
#> outcome_t_x1mis         42.213    1.43975  29.319  0.00000 ***
#> outcome_t_x2mis          2.623    1.66898   1.572  0.11605    
#> outcome_t_x3mis        -54.655   36.37724  -1.502  0.13298    
#> outcome_t_x4mis          0.330    0.02139  15.425  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  50.7071    17.5344   2.892 0.0038297  **
#> outcome_c_x1mis        41.1159     1.5239  26.981 0.0000000 ***
#> outcome_c_x2mis        -7.5699     1.9643  -3.854 0.0001163 ***
#> outcome_c_x3mis       134.7874    38.1127   3.537 0.0004054 ***
#> outcome_c_x4mis         0.3953     0.0165  23.958 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 280.46 
#>   For estimating E[Y(0)]:   350.35
```

### Misspecified propensity score model

``` r
# Distribution balancing weighting with normalization and without regularization
fitdbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
               estimand = "ATE", method = "dbw",
               method_y = "wls", data = df, normalize = TRUE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitdbwm)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "dbw", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    7.528     0.8596   8.757        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    218.4      1.258   173.6        0 ***
#> mu_c    210.9      1.234   171.0        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -4.054926   1.304617  -3.108 1.883e-03  **
#> ps_t_x1mis       -1.457935   0.255118  -5.715 1.099e-08 ***
#> ps_t_x2mis        0.349337   0.183647   1.902 5.714e-02   .
#> ps_t_x3mis       -0.585423   3.033180  -0.193 8.470e-01    
#> ps_t_x4mis        0.005433   0.002993   1.815 6.947e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -1.744982   1.858389 -0.9390 3.477e-01    
#> ps_c_x1mis       -2.953897   0.520234 -5.6780 1.363e-08 ***
#> ps_c_x2mis       -0.155648   0.298685 -0.5211 6.023e-01    
#> ps_c_x3mis       13.449387   6.430733  2.0914 3.649e-02   *
#> ps_c_x4mis        0.008381   0.003344  2.5065 1.219e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_t_(Intercept)  17.6337     11.961   1.474 0.1404136    
#> outcome_t_x1mis        42.5669      1.157  36.793 0.0000000 ***
#> outcome_t_x2mis         4.1683      1.516   2.749 0.0059748  **
#> outcome_t_x3mis       -95.8512     28.878  -3.319 0.0009029 ***
#> outcome_t_x4mis         0.3258      0.021  15.511 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  60.7096   14.69805   4.130 3.621e-05 ***
#> outcome_c_x1mis        42.8936    1.82753  23.471 0.000e+00 ***
#> outcome_c_x2mis        -7.3170    1.90618  -3.839 1.237e-04 ***
#> outcome_c_x3mis        63.1911   35.87621   1.761 7.818e-02   .
#> outcome_c_x4mis         0.3978    0.01812  21.953 0.000e+00 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 294.23 
#>   For estimating E[Y(0)]:   324.02

# Covariate balancing weighting function without regularization
fitcbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
               estimand = "ATE", method = "cb",
               method_y = "wls", data = df, normalize = TRUE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitcbwm)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "cb", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value  Pr(>|z|)    
#> ATE     6.25     0.9597   6.512 7.392e-11 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    218.2      1.274   171.3        0 ***
#> mu_c    212.0      1.244   170.4        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error  z value  Pr(>|z|)    
#> ps_t_(Intercept) -4.450192   1.291912 -3.44466 5.718e-04 ***
#> ps_t_x1mis       -1.557728   0.210025 -7.41688 1.199e-13 ***
#> ps_t_x2mis        0.420731   0.169312  2.48494 1.296e-02   *
#> ps_t_x3mis        0.080172   2.685062  0.02986 9.762e-01    
#> ps_t_x4mis        0.004619   0.002296  2.01164 4.426e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_c_(Intercept) -3.556914   1.680890  -2.116 0.034337   *
#> ps_c_x1mis       -2.404607   0.262851  -9.148 0.000000 ***
#> ps_c_x2mis        0.226348   0.211455   1.070 0.284426    
#> ps_c_x3mis        7.621652   3.921787   1.943 0.051966   .
#> ps_c_x4mis        0.005427   0.002098   2.588 0.009667  **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)  16.1437    12.1848   1.325 0.185204    
#> outcome_t_x1mis        41.9532     1.1070  37.900 0.000000 ***
#> outcome_t_x2mis         4.1782     1.5489   2.697 0.006987  **
#> outcome_t_x3mis       -86.6533    29.0835  -2.979 0.002888  **
#> outcome_t_x4mis         0.3256     0.0213  15.285 0.000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  53.2247   16.40538   3.244 0.0011772  **
#> outcome_c_x1mis        41.1675    1.72887  23.812 0.0000000 ***
#> outcome_c_x2mis        -7.1051    1.95576  -3.633 0.0002803 ***
#> outcome_c_x3mis       108.0243   38.80630   2.784 0.0053746  **
#> outcome_c_x4mis         0.3946    0.01701  23.194 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 267.85 
#>   For estimating E[Y(0)]:   379.81

# Entropy balancing weighting function without regularization
fitebwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
               estimand = "ATE", method = "eb",
               method_y = "wls", data = df, normalize = TRUE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitebwm)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "eb", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value  Pr(>|z|)    
#> ATE    4.949      1.126   4.396 1.102e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    219.1      1.338   163.8        0 ***
#> mu_c    214.1      1.268   168.8        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -2.388696   0.654434 -3.6500 2.622e-04 ***
#> ps_t_x1mis       -1.011005   0.141909 -7.1243 1.046e-12 ***
#> ps_t_x2mis        0.174877   0.084582  2.0676 3.868e-02   *
#> ps_t_x3mis       -0.612104   1.366100 -0.4481 6.541e-01    
#> ps_t_x4mis        0.002586   0.001203  2.1497 3.158e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_c_(Intercept) -2.380024   0.955776  -2.490  0.01277   *
#> ps_c_x1mis       -0.790713   0.084432  -9.365  0.00000 ***
#> ps_c_x2mis        0.296600   0.115485   2.568  0.01022   *
#> ps_c_x3mis        2.111323   1.289582   1.637  0.10159    
#> ps_c_x4mis        0.001558   0.000887   1.757  0.07892   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)   8.3954   12.97611   0.647  0.51764    
#> outcome_t_x1mis        42.3283    1.30648  32.399  0.00000 ***
#> outcome_t_x2mis         4.6286    1.58271   2.925  0.00345  **
#> outcome_t_x3mis       -69.9127   32.55074  -2.148  0.03173   *
#> outcome_t_x4mis         0.3257    0.02109  15.443  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  28.5177   16.70750   1.707 0.0878440   .
#> outcome_c_x1mis        44.2222    1.86708  23.685 0.0000000 ***
#> outcome_c_x2mis        -4.3555    1.89308  -2.301 0.0214052   *
#> outcome_c_x3mis       123.8418   36.62654   3.381 0.0007217 ***
#> outcome_c_x4mis         0.3757    0.01516  24.782 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 308.72 
#>   For estimating E[Y(0)]:   431.79

# Standard logistic regression
fitmlem <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
               estimand = "ATE", method = "mle",
               method_y = "wls", data = df, normalize = FALSE, 
               vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlem)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, normalize = FALSE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value  Pr(>|z|)    
#> ATE     5.09      1.164   4.374 1.217e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    217.9      1.296   168.1        0 ***
#> mu_c    212.8      1.307   162.8        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_t_(Intercept) -4.475967   1.450220 -3.0864 0.002026  **
#> ps_t_x1mis       -2.093636   0.213871 -9.7893 0.000000 ***
#> ps_t_x2mis        0.496130   0.178368  2.7815 0.005411  **
#> ps_t_x3mis        1.054133   1.958261  0.5383 0.590370    
#> ps_t_x4mis        0.003745   0.001803  2.0771 0.037788   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_c_(Intercept) -4.475967   1.450220 -3.0864 0.002026  **
#> ps_c_x1mis       -2.093636   0.213871 -9.7893 0.000000 ***
#> ps_c_x2mis        0.496130   0.178368  2.7815 0.005411  **
#> ps_c_x3mis        1.054133   1.958261  0.5383 0.590370    
#> ps_c_x4mis        0.003745   0.001803  2.0771 0.037788   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)   8.4033   13.44135  0.6252  0.53185    
#> outcome_t_x1mis        39.4293    1.18596 33.2467  0.00000 ***
#> outcome_t_x2mis         3.8572    1.72625  2.2345  0.02545   *
#> outcome_t_x3mis       -35.8454   29.09146 -1.2322  0.21789    
#> outcome_t_x4mis         0.3317    0.02286 14.5104  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  47.9952   18.91856   2.537 0.0111827   *
#> outcome_c_x1mis        39.5386    1.75123  22.578 0.0000000 ***
#> outcome_c_x2mis        -6.8059    2.16798  -3.139 0.0016935  **
#> outcome_c_x3mis       136.9302   35.34169   3.874 0.0001069 ***
#> outcome_c_x4mis         0.3912    0.01686  23.199 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 100.94 
#>   For estimating E[Y(0)]:   395.48

# Distribution balancing weighting without normalization and without regularization
fitdbwmnn <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
                 estimand = "ATE", method = "dbw",
                 method_y = "wls", data = df, normalize = FALSE, 
                 vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> Warning in dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand =
#> "ATE", : estimated weights are not normalized
summary(fitdbwmnn)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "dbw", method_y = "wls", data = df, normalize = FALSE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value  Pr(>|z|)    
#> ATE    6.266     0.7923   7.908 2.665e-15 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    218.2      1.250   174.5        0 ***
#> mu_c    212.0      1.246   170.2        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error  z value  Pr(>|z|)    
#> ps_t_(Intercept) -4.450192   1.334417 -3.33493 8.532e-04 ***
#> ps_t_x1mis       -1.557728   0.238994 -6.51786 7.132e-11 ***
#> ps_t_x2mis        0.420731   0.191730  2.19439 2.821e-02   *
#> ps_t_x3mis        0.080172   3.318517  0.02416 9.807e-01    
#> ps_t_x4mis        0.004619   0.003132  1.47482 1.403e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -3.582849   1.961983 -1.8261 6.783e-02   .
#> ps_c_x1mis       -2.411189   0.413959 -5.8247 5.721e-09 ***
#> ps_c_x2mis        0.230817   0.287529  0.8028 4.221e-01    
#> ps_c_x3mis        7.631834   6.079918  1.2553 2.094e-01    
#> ps_c_x4mis        0.005414   0.002753  1.9665 4.925e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)  16.1437   12.26751   1.316 0.188185    
#> outcome_t_x1mis        41.9532    1.10048  38.123 0.000000 ***
#> outcome_t_x2mis         4.1782    1.54374   2.707 0.006799  **
#> outcome_t_x3mis       -86.6533   29.25821  -2.962 0.003060  **
#> outcome_t_x4mis         0.3256    0.02124  15.333 0.000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  53.5921   15.20168   3.525 0.0004228 ***
#> outcome_c_x1mis        41.2237    1.78932  23.039 0.0000000 ***
#> outcome_c_x2mis        -7.1370    1.89505  -3.766 0.0001658 ***
#> outcome_c_x3mis       107.0804   36.27837   2.952 0.0031610  **
#> outcome_c_x4mis         0.3948    0.01721  22.944 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 267.85 
#>   For estimating E[Y(0)]:   377.78

# Distribution balancing weighting with normalization and with regularization
# Standardization
res_std_comp <- std_comp(formula_y = formula_y, 
                         formula_ps = formula_ps_m, 
                         estimand = "ATE", method_y = "wls", 
                         data = df, std = TRUE,
                         weights = NULL)
# Estimation
fitdbwmr <- dbw(formula_y = formula_y, 
                formula_ps = res_std_comp$formula_ps, 
                estimand = "ATE", method = "dbw", method_y = "wls",
                data = res_std_comp$data, normalize = TRUE, 
                vcov = TRUE, lambda = 0.01, weights = res_std_comp$weights, 
                clevel = 0.95)
summary(fitdbwmr)
#> Warning in summary.dbw(fitdbwmr): When lambda > 0, the uncertainty estimates may not be appropriate
#>   because the normal approximation may not hold.
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = res_std_comp$formula_ps, 
#>     estimand = "ATE", method = "dbw", method_y = "wls", data = res_std_comp$data, 
#>     normalize = TRUE, vcov = TRUE, lambda = 0.01, weights = res_std_comp$weights, 
#>     clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    6.917     0.8459   8.177 2.22e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    218.4      1.260   173.3        0 ***
#> mu_c    211.5      1.229   172.0        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -0.21421    0.04334  -4.943 7.695e-07 ***
#> ps_t_x1mis       -0.87989    0.15372  -5.724 1.040e-08 ***
#> ps_t_x2mis        0.18843    0.09997   1.885 5.944e-02   .
#> ps_t_x3mis       -0.02522    0.12610  -0.200 8.415e-01    
#> ps_t_x4mis        0.28949    0.16259   1.781 7.499e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error  z value  Pr(>|z|)    
#> ps_c_(Intercept) -0.28641     0.1369 -2.09266 3.638e-02   *
#> ps_c_x1mis       -1.55367     0.2689 -5.77850 7.537e-09 ***
#> ps_c_x2mis       -0.00729     0.1562 -0.04667 9.628e-01    
#> ps_c_x3mis        0.51875     0.2606  1.99043 4.654e-02   *
#> ps_c_x4mis        0.37925     0.1618  2.34329 1.911e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)  218.425     0.7296 299.364 0.000000 ***
#> outcome_t_x1mis         26.283     0.7252  36.244 0.000000 ***
#> outcome_t_x2mis          2.306     0.8361   2.758 0.005811  **
#> outcome_t_x3mis         -4.135     1.2223  -3.383 0.000717 ***
#> outcome_t_x4mis         18.085     1.1609  15.578 0.000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  211.508     0.6853 308.648 0.0000000 ***
#> outcome_c_x1mis         25.961     1.1048  23.499 0.0000000 ***
#> outcome_c_x2mis         -3.884     1.0543  -3.684 0.0002293 ***
#> outcome_c_x3mis          3.622     1.5336   2.362 0.0181903   *
#> outcome_c_x4mis         21.962     0.9754  22.515 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 303.63 
#>   For estimating E[Y(0)]:   363.1

# Covariate balancing weighting function with an estimating equation 
#  for the original covariate balancing propensity score method
fitcbwmcmb <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
                  estimand = "ATEcombined", method = "cb",
                  method_y = "wls", data = df, normalize = TRUE, 
                  vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitcbwmcmb)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATEcombined", 
#>     method = "cb", method_y = "wls", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATEcombined estimate:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATEcombined    4.766      1.113   4.283 1.848e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    218.0      1.260   173.1        0 ***
#> mu_c    213.2      1.326   160.8        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation:
#>                 Estimate Std. Error z value  Pr(>|z|)    
#> ps_(Intercept) -4.664229   1.369549 -3.4057 6.600e-04 ***
#> ps_x1mis       -1.708975   0.215824 -7.9184 2.442e-15 ***
#> ps_x2mis        0.434015   0.181683  2.3889 1.690e-02   *
#> ps_x3mis        1.999810   3.094182  0.6463 5.181e-01    
#> ps_x4mis        0.004374   0.002487  1.7589 7.860e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)  15.2506   12.45998   1.224  0.22096    
#> outcome_t_x1mis        40.8103    0.98965  41.237  0.00000 ***
#> outcome_t_x2mis         4.0300    1.60925   2.504  0.01227   *
#> outcome_t_x3mis       -71.3624   28.31457  -2.520  0.01172   *
#> outcome_t_x4mis         0.3258    0.02199  14.815  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)   41.984    21.7850   1.927 5.396e-02   .
#> outcome_c_x1mis         39.069     2.0592  18.973 0.000e+00 ***
#> outcome_c_x2mis         -6.230     2.4302  -2.564 1.036e-02   *
#> outcome_c_x3mis        148.382    37.4733   3.960 7.505e-05 ***
#> outcome_c_x4mis          0.388     0.0169  22.958 0.000e+00 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATEcombined estimation:
#>   For estimating E[Y(1)]: 232.61 
#>   For estimating E[Y(0)]:   426.14


# Formula for a misspecified outcome model for the GAM
library(mgcv)
#> Loading required package: nlme
#> This is mgcv 1.8-31. For overview type 'help("mgcv-package")'.
formula_y_gam <- stats::as.formula(y ~ s(x1mis) + s(x2mis) + 
                                       s(x3mis)+ s(x4mis))

# Distribution balancing weighting with the GAM
fitdbwmg <- dbw(formula_y = formula_y_gam, formula_ps = formula_ps_m, 
                estimand = "ATE", method = "dbw",
                method_y = "gam", data = df, normalize = TRUE, 
                vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> Warning in dbw(formula_y = formula_y_gam, formula_ps = formula_ps_m, estimand =
#> "ATE", : variance is estimated only when "method_y" is "wls" or "logit"
summary(fitdbwmg)
#> Warning in summary.dbw(fitdbwmg): "varcov" is empty
#> 
#> Call: 
#> dbw(formula_y = formula_y_gam, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "dbw", method_y = "gam", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ATE estimate:  8.688
#> 
#> Estimate of E[Y(1)]:  219.8
#> Estimate of E[Y(0)]:  211.1
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#> ps_t_(Intercept)       ps_t_x1mis       ps_t_x2mis       ps_t_x3mis 
#>        -4.054926        -1.457935         0.349337        -0.585423 
#>       ps_t_x4mis 
#>         0.005433 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#> ps_c_(Intercept)       ps_c_x1mis       ps_c_x2mis       ps_c_x3mis 
#>        -1.744982        -2.953897        -0.155648        13.449387 
#>       ps_c_x4mis 
#>         0.008381 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#> NULL
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#> NULL
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 294.23 
#>   For estimating E[Y(0)]:   324.02
```

### Binary outcome case

``` r
# Empirically correct ATE
ybinom_t <- (y + (1 - treat) * 10 > 210) + 0
ybinom_c <- (y - treat * 10 > 210) + 0
ATEbinom <- mean(ybinom_t - ybinom_c)
ATEbinom
#> [1] 0.127

# Formula for a misspecified binary outcome model
formula_y_bin <- stats::as.formula(ybinom ~ x1mis + x2mis + x3mis + x4mis)

# Distribution balancing weighting for the binary outcome
fitdbwmbin <- dbw(formula_y = formula_y_bin, formula_ps = formula_ps_m, 
                  estimand = "ATE", method = "dbw",
                  method_y = "logit", data = df, normalize = TRUE, 
                  vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitdbwmbin)
#> 
#> Call: 
#> dbw(formula_y = formula_y_bin, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "dbw", method_y = "logit", data = df, normalize = TRUE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value  Pr(>|z|)    
#> ATE   0.1266    0.02028   6.244 4.276e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t   0.6233    0.01973   31.59        0 ***
#> mu_c   0.4967    0.01810   27.44        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -4.054926   1.304617  -3.108 1.883e-03  **
#> ps_t_x1mis       -1.457935   0.255118  -5.715 1.099e-08 ***
#> ps_t_x2mis        0.349337   0.183647   1.902 5.714e-02   .
#> ps_t_x3mis       -0.585423   3.033180  -0.193 8.470e-01    
#> ps_t_x4mis        0.005433   0.002993   1.815 6.947e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -1.744982   1.858389 -0.9390 3.477e-01    
#> ps_c_x1mis       -2.953897   0.520234 -5.6780 1.363e-08 ***
#> ps_c_x2mis       -0.155648   0.298685 -0.5211 6.023e-01    
#> ps_c_x3mis       13.449387   6.430733  2.0914 3.649e-02   *
#> ps_c_x4mis        0.008381   0.003344  2.5065 1.219e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                        Estimate Std. Error z value  Pr(>|z|)    
#> outcome_t_(Intercept) -21.80313   3.339793 -6.5283 6.653e-11 ***
#> outcome_t_x1mis         8.25570   1.272638  6.4871 8.752e-11 ***
#> outcome_t_x2mis         0.25006   0.380657  0.6569 5.112e-01    
#> outcome_t_x3mis       -25.56250   8.732696 -2.9272 3.420e-03  **
#> outcome_t_x4mis         0.04372   0.006128  7.1344 9.723e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept) -33.5006   6.032885  -5.553 2.808e-08 ***
#> outcome_c_x1mis         8.7056   1.213130   7.176 7.170e-13 ***
#> outcome_c_x2mis        -0.5904   0.569327  -1.037 2.997e-01    
#> outcome_c_x3mis        18.6834   9.153927   2.041 4.125e-02   *
#> outcome_c_x4mis         0.0647   0.009179   7.048 1.810e-12 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 294.23 
#>   For estimating E[Y(0)]:   324.02
```

### Horvitz-Thompson and Hajek estimators

``` r
# Standard logistic regression with the Horvitz-Thompson estimator
fitmlem_ht <- dbw(formula_y = y ~ 0, formula_ps = formula_ps_m, 
                  estimand = "ATE", method = "mle",
                  method_y = "wls", data = df, normalize = FALSE, 
                  vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlem_ht)
#> 
#> Call: 
#> dbw(formula_y = y ~ 0, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, normalize = FALSE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)  
#> ATE    51.47         24   2.145  0.03197 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    261.3     23.152   11.29        0 ***
#> mu_c    209.9      2.269   92.50        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_t_(Intercept) -4.475967   1.450220 -3.0864 0.002026  **
#> ps_t_x1mis       -2.093636   0.213871 -9.7893 0.000000 ***
#> ps_t_x2mis        0.496130   0.178368  2.7815 0.005411  **
#> ps_t_x3mis        1.054133   1.958261  0.5383 0.590370    
#> ps_t_x4mis        0.003745   0.001803  2.0771 0.037788   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_c_(Intercept) -4.475967   1.450220 -3.0864 0.002026  **
#> ps_c_x1mis       -2.093636   0.213871 -9.7893 0.000000 ***
#> ps_c_x2mis        0.496130   0.178368  2.7815 0.005411  **
#> ps_c_x3mis        1.054133   1.958261  0.5383 0.590370    
#> ps_c_x4mis        0.003745   0.001803  2.0771 0.037788   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 100.94 
#>   For estimating E[Y(0)]:   395.48

# Standard logistic regression with the Hajek estimator
fitmlem_hj <- dbw(formula_y = y ~ 1, formula_ps = formula_ps_m, 
                  estimand = "ATE", method = "mle",
                  method_y = "wls", data = df, normalize = FALSE, 
                  vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlem_hj)
#> 
#> Call: 
#> dbw(formula_y = y ~ 1, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, normalize = FALSE, 
#>     vcov = TRUE, lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value  Pr(>|z|)    
#> ATE    14.57      4.174    3.49 0.0004822 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    227.6      4.294    53.0        0 ***
#> mu_c    213.0      1.498   142.2        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_t_(Intercept) -4.475967   1.450220 -3.0864 0.002026  **
#> ps_t_x1mis       -2.093636   0.213871 -9.7893 0.000000 ***
#> ps_t_x2mis        0.496130   0.178368  2.7815 0.005411  **
#> ps_t_x3mis        1.054133   1.958261  0.5383 0.590370    
#> ps_t_x4mis        0.003745   0.001803  2.0771 0.037788   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> ps_c_(Intercept) -4.475967   1.450220 -3.0864 0.002026  **
#> ps_c_x1mis       -2.093636   0.213871 -9.7893 0.000000 ***
#> ps_c_x2mis        0.496130   0.178368  2.7815 0.005411  **
#> ps_c_x3mis        1.054133   1.958261  0.5383 0.590370    
#> ps_c_x4mis        0.003745   0.001803  2.0771 0.037788   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)    227.6      4.294      53        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_c_(Intercept)      213      1.498   142.2        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 100.94 
#>   For estimating E[Y(0)]:   395.48
```

## Covariate balance plot

``` r
# Distribution balancing weighting without regularization
oldpar <- par(no.readonly = TRUE) # Just for adjusting plot margins
par(mar = c(5.1, 5.1, 4.1, 2.1)) # Just for adjusting plot margins
plot(fitdbwm, addcov = ~ x1 + x2 + x3 + x4, threshold = c(0.1, 0.2))
```

<img src="man/figures/README-covariate balance-1.png" width="100%" />

``` r
par(oldpar)

# Display the covariate balance matrix
cb <- plot(fitdbwm, addcov = ~ x1 + x2 + x3 + x4, plot = FALSE)
cb
#>    covariates       type     diff.adj     diff.un
#> 1 (Intercept) continuous  0.000000000  0.00000000
#> 2       x1mis continuous -0.007119096 -0.80804880
#> 3       x2mis continuous  0.039028380  0.30242521
#> 4       x3mis continuous -0.039931997  0.05771619
#> 5       x4mis continuous -0.010003270  0.30923044
#> 6          x1 continuous -0.060863139 -0.88090056
#> 7          x2 continuous  0.135848558  0.37398901
#> 8          x3 continuous -0.100778245 -0.13194982
#> 9          x4 continuous -0.158338606  0.05416632
```
