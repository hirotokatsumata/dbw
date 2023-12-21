#' Generate standardized complete cases
#'
#' Standardizes covariates for propensity score estimation before the 
#' distribution balancing weighting \code{\link{dbw}} with the regularization.
#'
#' \code{std_comp} first extracts complete cases for both propensity score estimation
#' and outcome model estimation. Then it standardizes covariates for propensity 
#' score estimation by takeing the provided weights into account. The returned
#' data frame is the "design matrix", which contains a set of dummy variables 
#' (depending on the contrasts) for factors and similarly expanded interaction
#' terms.
#'
#' For the AO estimation, NA values for the outcome variable for missing cases
#' (the response variable taking "0") are not deleted. For this processing, the 
#' outcome variable name must not contain spaces.
#' 
#' @export
#'
#' @param formula_y an object of class \code{\link[stats]{formula}} (or one that
#'   can be coerced to that class): a symbolic description of the potential 
#'   outcome model to be fitted. When you want to use non-DR type estimators,
#'   only include "1" in the right hand side of the formula for the Hajek
#'   estimator and only include "0" for the Horvitz-Thompson estimator. 
#'   See example of \code{\link{dbw}} for more details.
#' @param formula_ps an object of class \code{\link[stats]{formula}} (or one that
#'   can be coerced to that class): a symbolic description of the propensity score
#'   model to be fitted. For the entropy balancing weighting (\code{method = "eb"}), 
#'   variables in the right hand side of the formula will be mean balanced.
#' @param estimand a character string specifying a parameter of interest. Choose
#'   "ATT" for the average treatment effects on the treated estimation, "ATE"
#'   for the average treatment effects estimation, "ATC" for the average 
#    treatment effects on the controlled estimation, or "AO" for the average
#'   outcomes estimation in missing outcome cases. You can choose "ATEcombined" 
#'   for the combined estimation for the average treatment effects estimation 
#'   when using the covariate balancing weighting (\code{method = "cb"}).
#' @param method_y a character string specifying a method for potential outcome 
#'   prediction. Choose "wls" for the linear model, "logit" for the logistic
#'   regression, "gam" for the generalized additive model for the continuous 
#'   outcome, and "gambinom" for the generalized additive model for the binary 
#'   outcome.
#' @param data a data frame (or one that can be coerced to that class) 
#'   containing the outcomes and the variables in the model.
#' @param std a logical parameter indicating whether to standardize the 
#'   covariates for propensity score estimation.
#' @param weights an optional vector of ‘prior weights’ (e.g. sampling weights)
#'   to be used in the fitting process. Should be NULL or a numeric vector.
#'
#' @return
#' \item{data}{the complete-case data frame containing the outcome variable, 
#'   the response (treatment) variable, and the standardized covariates for 
#'   propensity score estimation.}
#' \item{weights}{the initially-supplied weights for the complete cases, a 
#'	 vector of 1s if none were.}
#' \item{formula_ps}{the automatically-created formula for propensity score 
#'	 estimation, which includes expanded factors and interactions.}
#' \item{mean_x}{the weighted mean of each covariates.}
#' \item{sd_x}{the weighted standard deviation of each covariates.}
#'
#' @author Hiroto Katsumata
#'
#' @seealso \code{\link{dbw}}
#'
#' @examples
#' # Simulation from Kang and Shafer (2007) and Imai and Ratkovic (2014)
#' # ATE estimation
#' # True ATE is 10
#' tau <- 10
#' set.seed(12345)
#' n <- 1000
#' X <- matrix(stats::rnorm(n * 4, mean = 0, sd = 1), nrow = n, ncol = 4)
#' prop <- 1 / (1 + exp(X[, 1] - 0.5 * X[, 2] + 
#'                      0.25 * X[, 3] + 0.1 * X[, 4]))
#' treat <- rbinom(n, 1, prop)
#' y <- 210 + 
#'      27.4 * X[, 1] + 13.7 * X[, 2] + 13.7 * X[, 3] + 13.7 * X[, 4] + 
#'      tau * treat + stats::rnorm(n = n, mean = 0, sd = 1)
#' ybinom <- (y > 210) + 0
#' df0 <- data.frame(X, treat, y, ybinom)
#' colnames(df0) <- c("x1", "x2", "x3", "x4", "treat", "y", "ybinom")
#'
#' # Variables for a misspecified model
#' Xmis <- data.frame(x1mis = exp(X[, 1] / 2), 
#'                    x2mis = X[, 2] * (1 + exp(X[, 1]))^(-1) + 10,
#'                    x3mis = (X[, 1] * X[, 3] / 25 + 0.6)^3, 
#'                    x4mis = (X[, 2] + X[, 4] + 20)^2)
#'
#' # Data frame and formulas for propensity score estimation
#' df <- data.frame(df0, Xmis)
#' formula_ps_c <- stats::as.formula(treat ~ x1 + x2 + x3 + x4)
#' formula_ps_m <- stats::as.formula(treat ~ x1mis + x2mis + 
#'                                           x3mis + x4mis)
#'
#' # Formula for a misspecified outcome model
#' formula_y <- stats::as.formula(y ~ x1mis + x2mis + x3mis + x4mis)
#'
#'
#' # Distribution balancing weighting with regularization
#' # Standardization
#' res_std_comp <- std_comp(formula_y = formula_y, 
#'                          formula_ps = formula_ps_m, 
#'                          estimand = "ATE", method_y = "wls", 
#'                          data = df, std = TRUE,
#'                          weights = NULL)
#' fitdbwmr <- dbw(formula_y = formula_y, 
#'                 formula_ps = res_std_comp$formula_ps, 
#'                 estimand = "ATE", method = "dbw", method_y = "wls",
#'                 data = res_std_comp$data, vcov = TRUE, 
#'                 lambda = 0.01, weights = res_std_comp$weights, 
#'                 clevel = 0.95)
#' summary(fitdbwmr)
std_comp <- function (formula_y, formula_ps, estimand = "ATE", method_y = "wls",
	                    data, std = TRUE, weights = NULL) {
  # Check
  if (estimand %in% c("ATE", "ATT", "ATC", "AO", "ATEcombined") == 0) {
    stop("estimand must be \"ATE\", \"ATT\", \"ATC\", 
         \"AO\", or \"ATEcombined\"")
  }
  if (method_y %in% c("wls", "logit", "gam", "gambinom") == 0) {
    stop("method must be \"wls\", \"logit\", \"gam\", or \"gambinom\"")
  }
  data <- data.frame(data)
  formula_y <- stats::as.formula(formula_y)
  formula_ps <- stats::as.formula(formula_ps)
  if (method_y %in% c("gam", "gambinom")) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop(
        "Package \"mgcv\" must be installed to use \"gam\" or \"gambinom\".",
        call. = FALSE
      )
    }
    formula_y0 <- formula_y
    formula_y_chr <- deparse1(formula_y)
    names_y <- substr(x = formula_y_chr, 
                      start = 1, 
                      stop = regexpr(pattern = "~", text = formula_y_chr) - 1)
    res_gam <- mgcv::gam(formula = formula_y, data = data, fit = FALSE)
    formula_y <- paste(names_y, deparse1(res_gam$pred.formula))
    formula_y <- stats::as.formula(formula_y)
  }
  if (estimand == "AO") {
    formula_y_chr <- deparse1(formula_y)
    names_y <- substr(x = formula_y_chr, 
                      start = 1, 
                      stop = regexpr(pattern = "~", text = formula_y_chr) - 1)
    names_y <- gsub(pattern = " ", replacement = "", x = names_y)
  }
  model_ps <- stats::model.frame(formula_ps, data = data)
  incomplete_model_ps <- c(attr(model_ps, which = "na.action"))
  complete_model_ps <- setdiff(c(1:nrow(data)), incomplete_model_ps)
  data_ps <- data[complete_model_ps, ]
  if (estimand == "AO") {
    data_ps2 <- data_ps
    response_temp <- c(stats::model.extract(model_ps, "response"))
    data_ps[response_temp == 0, names_y] <- 0
  }
  model_y <- stats::model.frame(formula_y, data = data_ps)
  incomplete_model_y <- c(attr(model_y, which = "na.action"))
  complete_model_y <- setdiff(c(1:nrow(data_ps)), incomplete_model_y)
  data <- data_ps[complete_model_y, ]
  if (estimand == "AO") {
    data <- data_ps2[complete_model_y, ]
    model_y <- stats::model.frame(formula_y, data = data, 
    	                            na.action = stats::na.pass)
  }
  names_y <- names(model_y)[1]
  y <- c(stats::model.extract(model_y, "response"))
  attr(y, which = "names") <- NULL
  model_ps <- stats::model.frame(formula_ps, data = data, 
                                 na.action = stats::na.pass)
  names_response <- names(model_ps)[1]
  response <- c(stats::model.extract(model_ps, "response"))
  attr(response, which = "names") <- NULL
  x_ps <- as.matrix(stats::model.matrix(formula_ps, model_ps))
  N <- nrow(data)
  if (setequal(response, c(0, 1)) == FALSE) {
    stop("treatment (response) variable must be binary (0, 1)")
  }
  if (is.null(weights) == 1) {
    weights <- rep(x = 1, times = N)
  } else {
    weights <- (weights[complete_model_ps])[complete_model_y]
  }
  if (length(weights) != N) {
    stop("length of weights must be the same as the number of rows of data")
  }
  model_mat <- stats::model.matrix(formula_ps, data = data)
  model_mat <- data.frame(model_mat)
  names_x <- as.character(names(model_mat))
  names_x <- names_x[(names_x != "X.Intercept.")]
  covariates <- as.matrix(model_mat[, names_x])
  mean_x <- numeric(ncol(covariates))
  sd_x <- numeric(ncol(covariates))
  m <- sum(weights > 0)
  for (i in 1:ncol(covariates)) {
    mean_x[i] <- sum(covariates[, i] * weights) / sum(weights)
    sd_x[i] <- sqrt(sum((covariates[, i] - mean_x[i])^2 * weights) / 
                    (sum(weights) * (m - 1) /  m))
  }
  if (std == TRUE) {
  	mean_mat <- matrix(rep(x = mean_x, each = nrow(covariates)), 
	  	                 nrow = nrow(covariates), 
		                   ncol = ncol(covariates))
  	sd_mat <- matrix(rep(x = sd_x, each = nrow(covariates)), 
	  	               nrow = nrow(covariates), 
		                 ncol = ncol(covariates))
  	covariates <- (covariates - mean_mat) / sd_mat
  }
  data <- data.frame(y, response)
  colnames(data) <- c(names_y, names_response)
  data <- data.frame(data, covariates)
  colnames(data)[-c(1, 2)] <- names_x
  names(mean_x) <- names_x
  names(sd_x) <- names_x
  formula_new <- paste(names_response, " ~ ", paste(names_x, collapse = " + "))
  formula_new <- stats::as.formula(formula_new)
  list(data = data,
  	   weights = weights,
  	   formula_ps = formula_new,
  	   mean_x = mean_x,
  	   sd_x = sd_x)
}
