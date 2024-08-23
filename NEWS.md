# dbw 1.1.4
### Minor changes
* replace dontrun{} with donttest{} in an example code in dbw() function

# dbw 1.1.3
### Minor changes
* Add a citation of the article proposing the distribution balancing weighting (Katsumata 2024) published in Political Science Research and Methods

# dbw 1.1.2
### Bug fixes
* Fix an issue in assigning column names for std_comp() function when there is only one covariate for the propensity score model.

# dbw 1.1.1
### Bug fixes
* Update examples for plot.dbw() function.

# dbw 1.1.0
### Minor changes
* Add "normalize" argument to dbw() function for "method = "dbw"". If set "FALSE", it estimates the non-normalized distribution balancing weights.
* Add internal functions for estimating the non-normalized distribution balancing weights.
* Add "normalize" component to dbw class.
* Add some examples for estimating non-normalized distribution balancing weights.
* Update the README files
### Bug fixes
* Revise error and warning messages.

# dbw 1.0.4
### Minor changes
* Delete unnecessary calculation in DC algorithm.
* Add explanation for "formula_y" not to include the treatment variable.
* Add an error message when "formula_y" includes the treatment variable.
* Fix typos in a warning message in dbw().
### Bug fixes
* Fix issues in calculating the derivatives of the regularization terms.
* Fix an issue in std_comp() function when the number of covariates is one.

# dbw 1.0.3
### Bug fixes
* Fix typos in DESCRIPTION file and an error in dbw().

# dbw 1.0.2
### Bug fixes
* Fix typos and delete unnecessary spaces in the DESCRIPTION file.

# dbw 1.0.1
### Minor changes
* Use dontrun() in the example for the time constraint.

# dbw 1.0.0
### Major changes
* Improve the DC algorithm for "method = "dbw"".
### Minor changes
* Revise the variance-covariance matrix estimation for "method = "dbw"".
* Add "tol" option, which is the tolerance parameter for "method = "dbw"".
* Add "init_lambda" option, which is a parameter for "method = "dbw"" to set the lambda value for the initial values estimation.

# dbw 0.5.0
### Minor changes
* Revert the internal calculation of eta when lambda > 0.

# dbw 0.4.1
### Minor changes
* Normalize the estimated weights in "plot.dbw()" function.

# dbw 0.4.0
### Minor changes
* Slightly change the internal calculation of eta when lambda > 0.

# dbw 0.3.0
### Minor changes
* Add "std_comp()" function, which generates a complete-case data frame with standardized covariates for propensity score estimation. This function is meant to be used before dbw() with regularization (lambda > 0).
### Bug fixes
* incomplete case adjustment for weights.

# dbw 0.2.2
### Bug fixes
* Scale checks when lambda > 0.

# dbw 0.2.1
### Bug fixes
* Outcome name extraction for "method = "AO"".

# dbw 0.2.0
### Minor changes
* Add "vcov" option, which indicates whether to estimate the variance.
* Delete unused "std" function.

# dbw 0.1.0
### Minor changes
* Add references.

# dbw 0.0.1.9000
### Major changes
* This is the first release of dbw.
