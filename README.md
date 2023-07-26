
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dbw: Doubly Robust Distribution Balancing Weighting Estimation

<!-- badges: start -->
<!-- badges: end -->

**dbw** implements the doubly robust distribution balancing weighting
proposed by Katsumata (2023), which improves the augmented inverse
probability weighting (AIPW) by estimating propensity scores with
estimating equations suitable for the pre-specified parameter of
interest (e.g., the average treatment effects or the average treatment
effects on the treated) and estimating outcome models with the estimated
inverse probability weights.

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

<font size="4"> Katsumata, Hiroto. 2023. “How Should We Estimate Inverse
Probability Weights with Possibly Misspecified Propensity Score Models?”
</font>

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
# Distribution balancing weighting without regularization
fitdbwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c, 
               estimand = "ATE", method = "dbw",
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
fitdbwc
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "dbw", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
#>     method = "dbw", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    10.03     0.6065   16.54        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    220.6      1.250   176.5        0 ***
#> mu_c    210.5      1.184   177.7        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -0.06034     0.1059 -0.5700 5.687e-01    
#> ps_t_x1          -1.17788     0.1727 -6.8204 9.077e-12 ***
#> ps_t_x2           0.63462     0.1274  4.9828 6.269e-07 ***
#> ps_t_x3          -0.25104     0.1390 -1.8054 7.101e-02   .
#> ps_t_x4           0.02505     0.1235  0.2029 8.393e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -0.06063     0.1120 -0.5412 5.884e-01    
#> ps_c_x1          -1.01973     0.1518 -6.7197 1.821e-11 ***
#> ps_c_x2           0.35814     0.1057  3.3873 7.059e-04 ***
#> ps_c_x3          -0.37044     0.1547 -2.3953 1.661e-02   *
#> ps_c_x4           0.24305     0.1190  2.0430 4.105e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)    25.80   12.75626   2.023  0.04308   *
#> outcome_t_x1mis          42.32    1.46038  28.979  0.00000 ***
#> outcome_t_x2mis           2.56    1.67496   1.528  0.12648    
#> outcome_t_x3mis         -56.27   36.56700  -1.539  0.12387    
#> outcome_t_x4mis           0.33    0.02152  15.335  0.00000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  47.6875   17.56453   2.715 0.0066278  **
#> outcome_c_x1mis        41.0927    1.51049  27.205 0.0000000 ***
#> outcome_c_x2mis        -7.1336    1.95902  -3.641 0.0002712 ***
#> outcome_c_x3mis       132.9864   38.77449   3.430 0.0006042 ***
#> outcome_c_x4mis         0.3921    0.01617  24.247 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 279.66 
#>   For estimating E[Y(0)]:   366.01

# Covariate balancing weighting function without regularization
fitcbwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c, 
               estimand = "ATE", method = "cb",
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
summary(fitcbwc)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "cb", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
summary(fitebwc)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "eb", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlec)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_c, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
```

### Misspecified propensity score model

``` r
# Distribution balancing weighting without regularization
fitdbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
               estimand = "ATE", method = "dbw",
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
summary(fitdbwm)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "dbw", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    7.528     0.8337   9.029        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    218.4      1.261   173.2        0 ***
#> mu_c    210.9      1.240   170.1        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -4.054927   1.279210 -3.1699 1.525e-03  **
#> ps_t_x1mis       -1.457935   0.231495 -6.2979 3.017e-10 ***
#> ps_t_x2mis        0.349337   0.181716  1.9224 5.455e-02   .
#> ps_t_x3mis       -0.585423   3.037680 -0.1927 8.472e-01    
#> ps_t_x4mis        0.005433   0.002941  1.8470 6.474e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -1.744982   1.913758 -0.9118 3.619e-01    
#> ps_c_x1mis       -2.953897   0.478347 -6.1752 6.607e-10 ***
#> ps_c_x2mis       -0.155648   0.316227 -0.4922 6.226e-01    
#> ps_c_x3mis       13.449387   6.485158  2.0739 3.809e-02   *
#> ps_c_x4mis        0.008381   0.003411  2.4573 1.400e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_t_(Intercept)  17.6337     11.973   1.473 0.1407976    
#> outcome_t_x1mis        42.5669      1.165  36.525 0.0000000 ***
#> outcome_t_x2mis         4.1683      1.516   2.749 0.0059777  **
#> outcome_t_x3mis       -95.8512     28.887  -3.318 0.0009062 ***
#> outcome_t_x4mis         0.3258      0.021  15.513 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  60.7096   14.69386   4.132 3.602e-05 ***
#> outcome_c_x1mis        42.8936    1.82903  23.452 0.000e+00 ***
#> outcome_c_x2mis        -7.3170    1.91333  -3.824 1.312e-04 ***
#> outcome_c_x3mis        63.1911   35.57815   1.776 7.571e-02   .
#> outcome_c_x4mis         0.3978    0.01815  21.915 0.000e+00 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 294.23 
#>   For estimating E[Y(0)]:   324.02

# Covariate balancing weighting function without regularization
fitcbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
               estimand = "ATE", method = "cb",
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
summary(fitcbwm)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "cb", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
summary(fitebwm)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "eb", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
               method_y = "wls", data = df, vcov = TRUE, 
               lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlem)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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

# Distribution balancing weighting with regularization
# Standardization
res_std_comp <- std_comp(formula_y = formula_y, 
                         formula_ps = formula_ps_m, 
                         estimand = "ATE", method_y = "wls", 
                         data = df, std = TRUE,
                         weights = NULL)
fitdbwmr <- dbw(formula_y = formula_y, 
                formula_ps = res_std_comp$formula_ps, 
                estimand = "ATE", method = "dbw", method_y = "wls",
                data = res_std_comp$data, vcov = TRUE, 
                lambda = 0.01, weights = res_std_comp$weights, 
                clevel = 0.95)
summary(fitdbwmr)
#> Warning in summary.dbw(fitdbwmr): When lambda > 0, the uncertainty estimates may not be appropriate
#>   because the normal approximation may not hold.
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = res_std_comp$formula_ps, 
#>     estimand = "ATE", method = "dbw", method_y = "wls", data = res_std_comp$data, 
#>     vcov = TRUE, lambda = 0.01, weights = res_std_comp$weights, 
#>     clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value Pr(>|z|)    
#> ATE    6.924     0.8382    8.26 2.22e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t    218.4      1.266   172.5        0 ***
#> mu_c    211.5      1.271   166.4        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                  Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept)  -0.2467    0.07443 -3.3148 9.171e-04 ***
#> ps_t_x1mis        -0.8770    0.13832 -6.3402 2.294e-10 ***
#> ps_t_x2mis         0.1906    0.09606  1.9839 4.726e-02   *
#> ps_t_x3mis        -0.0254    0.12319 -0.2062 8.367e-01    
#> ps_t_x4mis         0.2921    0.15354  1.9024 5.712e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                  Estimate Std. Error  z value  Pr(>|z|)    
#> ps_c_(Intercept) -0.35981     0.1884 -1.91018 5.611e-02   .
#> ps_c_x1mis       -1.56758     0.2805 -5.58798 2.297e-08 ***
#> ps_c_x2mis       -0.01004     0.1845 -0.05443 9.566e-01    
#> ps_c_x3mis        0.52870     0.2873  1.84018 6.574e-02   .
#> ps_c_x4mis        0.39297     0.1858  2.11459 3.446e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> outcome_t_(Intercept)  218.429     0.7348 297.271 0.000000 ***
#> outcome_t_x1mis         26.260     0.7236  36.291 0.000000 ***
#> outcome_t_x2mis          2.319     0.8374   2.770 0.005608  **
#> outcome_t_x3mis         -4.106     1.2229  -3.357 0.000787 ***
#> outcome_t_x4mis         18.087     1.1634  15.546 0.000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept)  211.505     0.7457 283.622 0.0000000 ***
#> outcome_c_x1mis         25.823     1.1242  22.969 0.0000000 ***
#> outcome_c_x2mis         -3.899     1.0598  -3.679 0.0002341 ***
#> outcome_c_x3mis          3.721     1.5223   2.444 0.0145149   *
#> outcome_c_x4mis         21.946     0.9753  22.502 0.0000000 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Effective sample size for the ATE estimation:
#>   For estimating E[Y(1)]: 301.7 
#>   For estimating E[Y(0)]:   366.39

# Covariate balancing weighting function with an estimating equation 
#  for the original covariate balancing propensity score method
fitcbwmcmb <- dbw(formula_y = formula_y, formula_ps = formula_ps_m, 
                  estimand = "ATEcombined", method = "cb",
                  method_y = "wls", data = df, vcov = TRUE, 
                  lambda = 0, weights = NULL, clevel = 0.95)
summary(fitcbwmcmb)
#> 
#> Call: 
#> dbw(formula_y = formula_y, formula_ps = formula_ps_m, estimand = "ATEcombined", 
#>     method = "cb", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
                method_y = "gam", data = df, vcov = TRUE, 
                lambda = 0, weights = NULL, clevel = 0.95)
#> Warning in dbw(formula_y = formula_y_gam, formula_ps = formula_ps_m, estimand =
#> "ATE", : variance is estimated only when method_y is "wls" or "logit"
summary(fitdbwmg)
#> Warning in summary.dbw(fitdbwmg): "varcov" is empty
#> 
#> Call: 
#> dbw(formula_y = formula_y_gam, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "dbw", method_y = "gam", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ATE estimate:  8.688
#> 
#> Estimate of E[Y(1)]:  219.8
#> Estimate of E[Y(0)]:  211.1
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#> ps_t_(Intercept)       ps_t_x1mis       ps_t_x2mis       ps_t_x3mis 
#>        -4.054927        -1.457935         0.349337        -0.585423 
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
                  method_y = "logit", data = df, vcov = TRUE, 
                  lambda = 0, weights = NULL, clevel = 0.95)
summary(fitdbwmbin)
#> 
#> Call: 
#> dbw(formula_y = formula_y_bin, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "dbw", method_y = "logit", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
#> 
#> ############################################################
#> ATE estimate:
#>     Estimate Std. Error z value  Pr(>|z|)    
#> ATE   0.1266    0.02024   6.255 3.983e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> ############################################################
#> 
#> Estimate of E[Y(1)] and E[Y(0)]:
#>      Estimate Std. Error z value Pr(>|z|)    
#> mu_t   0.6233    0.01973   31.59        0 ***
#> mu_c   0.4967    0.01812   27.42        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(1)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_t_(Intercept) -4.054927   1.279210 -3.1699 1.525e-03  **
#> ps_t_x1mis       -1.457935   0.231495 -6.2979 3.017e-10 ***
#> ps_t_x2mis        0.349337   0.181716  1.9224 5.455e-02   .
#> ps_t_x3mis       -0.585423   3.037680 -0.1927 8.472e-01    
#> ps_t_x4mis        0.005433   0.002941  1.8470 6.474e-02   .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the propensity score model estimation for estimating E[Y(0)]:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> ps_c_(Intercept) -1.744982   1.913758 -0.9118 3.619e-01    
#> ps_c_x1mis       -2.953897   0.478347 -6.1752 6.607e-10 ***
#> ps_c_x2mis       -0.155648   0.316227 -0.4922 6.226e-01    
#> ps_c_x3mis       13.449387   6.485158  2.0739 3.809e-02   *
#> ps_c_x4mis        0.008381   0.003411  2.4573 1.400e-02   *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(1)]:
#>                        Estimate Std. Error z value  Pr(>|z|)    
#> outcome_t_(Intercept) -21.80313   3.339835 -6.5282 6.656e-11 ***
#> outcome_t_x1mis         8.25570   1.272546  6.4875 8.725e-11 ***
#> outcome_t_x2mis         0.25006   0.380672  0.6569 5.113e-01    
#> outcome_t_x3mis       -25.56250   8.726163 -2.9294 3.396e-03  **
#> outcome_t_x4mis         0.04372   0.006127  7.1357 9.632e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Coefficients for the outcome model estimation for estimating E[Y(0)]:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> outcome_c_(Intercept) -33.5006   6.030989  -5.555 2.780e-08 ***
#> outcome_c_x1mis         8.7056   1.213498   7.174 7.283e-13 ***
#> outcome_c_x2mis        -0.5904   0.569420  -1.037 2.998e-01    
#> outcome_c_x3mis        18.6834   9.122549   2.048 4.056e-02   *
#> outcome_c_x4mis         0.0647   0.009183   7.045 1.849e-12 ***
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
                  method_y = "wls", data = df, vcov = TRUE, 
                  lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlem_ht)
#> 
#> Call: 
#> dbw(formula_y = y ~ 0, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
                  method_y = "wls", data = df, vcov = TRUE, 
                  lambda = 0, weights = NULL, clevel = 0.95)
summary(fitmlem_hj)
#> 
#> Call: 
#> dbw(formula_y = y ~ 1, formula_ps = formula_ps_m, estimand = "ATE", 
#>     method = "mle", method_y = "wls", data = df, vcov = TRUE, 
#>     lambda = 0, weights = NULL, clevel = 0.95)
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
#>    covariates       type      diff.adj     diff.un
#> 1 (Intercept) continuous  4.012344e-08  0.00000000
#> 2       x1mis continuous -7.119011e-03 -0.80804880
#> 3       x2mis continuous  3.902910e-02  0.30242521
#> 4       x3mis continuous -3.993179e-02  0.05771619
#> 5       x4mis continuous -1.000299e-02  0.30923044
#> 6          x1 continuous -6.086313e-02 -0.88090056
#> 7          x2 continuous  1.358486e-01  0.37398901
#> 8          x3 continuous -1.007782e-01 -0.13194982
#> 9          x4 continuous -1.583386e-01  0.05416632
```
