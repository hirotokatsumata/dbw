#' Doubly robust distribution balancing weighting (DBW) estimation
#'
#' \code{dbw} estimates a pre-specified parameter of interest (e.g.,
#' the average treatment effects (ATE) or the average treatment effects
#' on the treated (ATT)) with the augmented inverse probability weighting
#' (AIPW), where propensity scores are estimated using estimating
#' equations suitable for the parameter of interest and outcome models are
#' estimated using inverse probability weights. \code{dbw} can
#' also be used to estimate average outcomes (AO) in missing outcome cases.
#'
#' The treatment variable (or, response variable in missing outcome cases)
#' must be binary and coded as 0 (for controlled or missing observations)
#' or 1 (for treated or non-missing observations).
#'
#' When the data frame has incomplete cases, which have NAs for either of
#' the treatment variable, the outcome variable, or explanatory variables
#' either for propensity score or outcome model estimation, \code{dbw} conducts
#' listwise deletion. Returned values (e.g., \code{est_weights}, \code{ps}, \code{data})
#' do not contain values for these deleted cases.
#'
#' For propensity score estimation, \code{dbw} can utilize the distribution
#' balancing weighting (\code{method = "dbw"}), covariate balancing weighting
#' (\code{method = "cb"}), entropy balancing weighting (\code{method = "eb"}), or
#' standard maximum likelihood estimation (\code{method = "mle"}). For the
#' covariate balancing weighting and entropy balancing weighting, \code{dbw} runs
#' much faster than the original functions (\code{\link[CBPS]{CBPS}} and \code{\link[ebal]{ebalance}})
#' by using loss-function-based algorithms, which also results in more
#' accurate covariate balance. For the ATT and ATC estimation, the distribution
#' balancing weighting, covariate balancing weighting, and entropy balancing
#' weighting are theoretically equivalent and \code{dbw} implements accordingly.
#'
#' The parameter of interest is estimated by the AIPW estimator, where inverse
#' probability weights are standardized within each treatment group by being
#' devided by the size of the group after being calculated as
#' \eqn{t_i / \pi_i - (1 - t_i) / (1 - \pi_i)} for the ATE estimation,
#' \eqn{(t_i - \pi_i) / (1 - \pi_i)} for the ATT estimation,
#' \eqn{(t_i - \pi_i) / \pi_i} for the ATC estimation, and
#' \eqn{t_i / \pi_i} for the missing outcome cases. The resulting inverse probability
#' weights sum to 1 for the distribution balancing weighting, covariate
#' balancing weighting, and entropy balancing weighting estimators without
#' regularization.
#'
#' The variance-covariance matrix for the parameter of interest and ancillary
#' parameters is calculated using the sandwich variance formula obtained in the
#' M-estimation framework.
#'
#' When using regularization for propensity score estimation (\code{lambda > 0}),
#' you should standardize the covariates for propensity score estimation
#' by \code{\link{std_comp}} before using \code{dbw}. See example below for more details.
#'
#' For the ATE estimation, it is recommended to specify the \code{estimand} as
#' \code{"ATE"}, but you may specify it as \code{"ATEcombined"} when using the
#' covariate balancing weighting. The former utilizes the separated propensity
#' score estimation whereas the latter utilizes the combined estimation, and
#' the former should produce smaller biases and variances. Note that the
#' former estimates two propensity scores for each observation by estimating
#' two propensity score functions with different estimating equations.
#'
#' For the AO estimation, NA values for the outcome variable for missing cases
#' (the response variable taking "0") are not deleted. For this processing, the
#' outcome variable name must not contain spaces.
#'
#' @export
#'
#' @param formula_y an object of class \code{\link[stats]{formula}} (or one that
#'   can be coerced to that class): a symbolic description of the potential
#'   outcome model to be fitted. This is a model for potential outcomes, so do not
#'   include treatment variable in the model. When you want to use non-DR type 
#'   estimators, only include "1" in the right hand side of the formula for the 
#'   Hajek estimator and only include "0" for the Horvitz-Thompson estimator.
#'   See example below for more details.
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
#' @param method a character string specifying a method for propensity score
#'   estimation. Choose "dbw" for the distribution balancing weighting, "cb"
#'   for the covariate balancing weighting, "eb" for the entropy balancing
#'   weighting, and "mle" for the logistic regression with the maximum
#'   likelihood estimation.
#' @param method_y a character string specifying a method for potential outcome
#'   prediction. Choose "wls" for the linear model, "logit" for the logistic
#'   regression, "gam" for the generalized additive model for the continuous
#'   outcome, and "gambinom" for the generalized additive model for the binary
#'   outcome. Note that variance-covariance matrix is calculated only when
#'   \code{method_y = "wls"} or \code{method_y = "logit"}.
#' @param data a data frame (or one that can be coerced to that class)
#'   containing the outcomes and the variables in the model.
#' @param normalize a logical parameter indicating whether to normalize the estimated
#'   weights to sum up to one for each treatment group for \code{method = "dbw"}. Default is TRUE.
#' @param vcov a logical parameter indicating whether to estimate the variance.
#'   Default is TRUE.
#' @param lambda a parameter taking 0 or larger specifying the degree of the
#'   L2-regularization for propensity score estimation. \code{lambda = 0} (default)
#'   means no regularization.
#' @param weights an optional vector of ‘prior weights’ (e.g. sampling weights)
#'   to be used in the fitting process. Should be NULL or a numeric vector.
#' @param clevel confidence levels. Default is 0.95.
#' @param tol a tolerance parameter for \code{method = "dbw"}. Default is 1e-10.
#' @param init_lambda a parameter for \code{method = "dbw"} to set the lambda value
#'   for the initial values estimation, where \eqn{lambda_init = lambda * init_lambda}.
#'   Default is 0.01.
#'
#' @return \code{dbw} returns an object of "dbw" class.
#'
#' The function summary (i.e., \code{\link{summary.dbw}}) can be used to obtain or print a
#'   summary of the results.
#'
#' An object of class "dbw" is a list containing the following components:
#'
#' \item{est}{the point estimate of the parameter of interest.}
#' \item{coef_ps}{a named vector of coefficients for propensity score estimation.
#'   A list of two sets of coefficients for two sets of propensity scores
#'   (one for estimating \eqn{E[Y(1)]} and the other for estimating \eqn{E[Y(0)]}) is
#'   returned when \code{estimand = "ATE"}.}
#' \item{coef_y}{a named vector of coefficients for outcome model estimation.
#'   A list of two sets of coefficients for two sets of outcome models
#'   (one for estimating \eqn{E[Y(1)]} and the other for estimating \eqn{E[Y(0)]}) is
#'   returned when \code{estimand = "ATE"} and \code{estimand = "ATEcombined"}.}
#' \item{varcov}{the variance-covariance matrix of the coefficients and the
#'   parameter of interest.}
#' \item{est_weights}{the estimated inverse probability weights.}
#' \item{ps}{the estimated propensity scores. A list of two sets of the
#'   estimated propensity scores (one for estimating \eqn{E[Y(1)]} and the other
#'   for estimating \eqn{E[Y(0)])} is returned when \code{estimand = "ATE"}.}
#' \item{predicted_y}{the predicted outcomes. A list of two sets of the
#'   predicted outcomes (one for estimating \eqn{E[Y(1)]} and the other for
#'   estimating \eqn{E[Y(0)]}) is returned when \code{estimand = "ATE"} and
#'   \code{estimand = "ATEcombined"}.}
#' \item{converged}{logical. Were the propensity score estimation algorithms
#'   judged to have converged?}
#' \item{effn}{the effective sample size for the parameter of interest estimation.}
#' \item{effn_original}{the effective sample size with the initial weights.}
#' \item{estimand}{the parameter of interest specified.}
#' \item{method}{the method for propensity score estimation specified.}
#' \item{method_y}{the method for outcome model estimation specified.}
#' \item{response}{the treatment vector. The response (non-missingness) vector
#' when the missing outcome cases.}
#' \item{outcome}{the outcome vector.}
#' \item{original_weights}{the weights initially supplied, a vector of 1s if
#'   none were.}
#' \item{ci}{a matrix of the confidence intervals for the parameter of interest.}
#' \item{formula_y}{the outcome model formula specified.}
#' \item{formula_ps}{the propensity score model formula specified.}
#' \item{call}{the matched call.}
#' \item{data}{the data argument.}
#' \item{normalize}{a logical argument indicating whether to normalize the
#'	 estimated weights for each treatment group for \code{method = "dbw"}.}
#' \item{lambda}{the parameter specifying the degree of the L2-regularization.}
#'
#' @author Hiroto Katsumata
#'
#' @references Katsumata, Hiroto. 2024. "How Should We Estimate Inverse
#'   Probability Weights with Possibly Misspecified Propensity Score Models?"
#'   Political Science Research and Methods.
#'
#' Imai, Kosuke and Marc Ratkovic. 2014. "Covariate Balancing Propensity Score."
#'   Journal of the Royal Statistical Society, Series B (Statistical Methodology)
#'   76 (1): 243--63.
#'
#' Hainmueller, Jens. 2012. "Entropy Balancing for Causal Effects: A
#'   Multivariate Reweighting Method to Produce Balanced Samples in
#'   Observational Studies." Political Analysis 20 (1): 25--46.
#'
#' @seealso \code{\link{summary.dbw}}, \code{\link{std_comp}}, \code{\link[mgcv]{gam}}
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
#' # Correct propensity score model
#'
#' # Distribution balancing weighting with normalization and
#' #   without regularization
#' fitdbwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c,
#'                estimand = "ATE", method = "dbw",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' fitdbwc
#' summary(fitdbwc)
#'
#' # Covariate balancing weighting function without regularization
#' fitcbwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c,
#'                estimand = "ATE", method = "cb",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitcbwc)
#'
#' # Entropy balancing weighting function without regularization
#' fitebwc <- dbw(formula_y = formula_y, formula_ps = formula_ps_c,
#'                estimand = "ATE", method = "eb",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitebwc)
#'
#' # Standard logistic regression
#' fitmlec <- dbw(formula_y = formula_y, formula_ps = formula_ps_c,
#'                estimand = "ATE", method = "mle",
#'                method_y = "wls", data = df, normalize = FALSE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitmlec)
#'
#' # Distribution balancing weighting without normalization and
#' #   without regularization
#' fitdbwcnn <- dbw(formula_y = formula_y, formula_ps = formula_ps_c,
#'                  estimand = "ATE", method = "dbw",
#'                  method_y = "wls", data = df, normalize = FALSE,
#'                  vcov = TRUE, lambda = 0, weights = NULL,
#'                  clevel = 0.95)
#' summary(fitdbwcnn)
#'
#'
#' # Misspecified propensity score model
#'
#' # Distribution balancing weighting with normalization and 
#' #   without regularization
#' fitdbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "dbw",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitdbwm)
#'
#' # Covariate balancing weighting function without regularization
#' fitcbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "cb",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitcbwm)
#'
#' # Entropy balancing weighting function without regularization
#' fitebwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "eb",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitebwm)
#'
#' # Standard logistic regression
#' fitmlem <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "mle",
#'                method_y = "wls", data = df, normalize = FALSE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitmlem)
#'
#' # Distribution balancing weighting without normalization and 
#' #   without regularization
#' fitdbwmnn <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "dbw",
#'                method_y = "wls", data = df, normalize = FALSE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' summary(fitdbwmnn)
#'
#' # Distribution balancing weighting with normalization and 
#' #   with regularization
#' # Standardization
#' res_std_comp <- std_comp(formula_y = formula_y,
#'                          formula_ps = formula_ps_m,
#'                          estimand = "ATE", method_y = "wls",
#'                          data = df, std = TRUE,
#'                          weights = NULL)
#' # Estimation
#' fitdbwmr <- dbw(formula_y = formula_y,
#'                 formula_ps = res_std_comp$formula_ps,
#'                 estimand = "ATE", method = "dbw", method_y = "wls",
#'                 data = res_std_comp$data, normalize = TRUE,
#'                 vcov = TRUE, lambda = 0.01,
#'                 weights = res_std_comp$weights, clevel = 0.95)
#' summary(fitdbwmr)
#'
#' # Covariate balancing weighting function with an estimating equation
#' #  for the original covariate balancing propensity score method
#' fitcbwmcmb <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                   estimand = "ATEcombined", method = "cb",
#'                   method_y = "wls", data = df, normalize = TRUE,
#'                   vcov = TRUE, lambda = 0, weights = NULL,
#'                   clevel = 0.95)
#' summary(fitcbwmcmb)
#'
#'
#' \donttest{
#' # Formula for a misspecified outcome model for the GAM
#' library(mgcv)
#' formula_y_gam <- stats::as.formula(y ~ s(x1mis) + s(x2mis) +
#'                                        s(x3mis)+ s(x4mis))
#'
#' # Distribution balancing weighting with the GAM
#' fitdbwmg <- dbw(formula_y = formula_y_gam, formula_ps = formula_ps_m,
#'                 estimand = "ATE", method = "dbw",
#'                 method_y = "gam", data = df, normalize = TRUE,
#'                 vcov = TRUE, lambda = 0, weights = NULL,
#'                 clevel = 0.95)
#' summary(fitdbwmg)
#'
#'
#' # Binary outcome case
#' # Empirically correct ATE
#' ybinom_t <- (y + (1 - treat) * 10 > 210) + 0
#' ybinom_c <- (y - treat * 10 > 210) + 0
#' ATEbinom <- mean(ybinom_t - ybinom_c)
#' ATEbinom
#'
#' # Formula for a misspecified binary outcome model
#' formula_y_bin <- stats::as.formula(ybinom ~ x1mis + x2mis +
#'                                             x3mis + x4mis)
#'
#' # Distribution balancing weighting for the binary outcome
#' fitdbwmbin <- dbw(formula_y = formula_y_bin,
#'                   formula_ps = formula_ps_m,
#'                   estimand = "ATE", method = "dbw",
#'                   method_y = "logit", data = df, normalize = TRUE,
#'                   vcov = TRUE, lambda = 0, weights = NULL,
#'                   clevel = 0.95)
#' summary(fitdbwmbin)
#' }
#'
#' # Standard logistic regression with the Horvitz-Thompson estimator
#' fitmlem_ht <- dbw(formula_y = y ~ 0, formula_ps = formula_ps_m,
#'                   estimand = "ATE", method = "mle",
#'                   method_y = "wls", data = df, normalize = FALSE,
#'                   vcov = TRUE, lambda = 0, weights = NULL,
#'                   clevel = 0.95)
#' summary(fitmlem_ht)
#'
#' # Standard logistic regression with the Hajek estimator
#' fitmlem_hj <- dbw(formula_y = y ~ 1, formula_ps = formula_ps_m,
#'                   estimand = "ATE", method = "mle",
#'                   method_y = "wls", data = df, normalize = FALSE,
#'                   vcov = TRUE, lambda = 0, weights = NULL,
#'                   clevel = 0.95)
#' summary(fitmlem_hj)
dbw <- function (formula_y, formula_ps, estimand = "ATE", method = "dbw",
                 method_y = "wls", data, normalize = TRUE, vcov = TRUE, lambda = 0,
                 weights = NULL, clevel = 0.95, tol = 1e-10, init_lambda = 0.01) {
  # Check
  if (estimand %in% c("ATE", "ATT", "ATC", "AO", "ATEcombined") == 0) {
    stop("\"estimand\" must be \"ATE\", \"ATT\", \"ATC\",
         \"AO\", or \"ATEcombined\"")
  }
  if (estimand == "ATEcombined") {
    if (method != "cb") {
      stop("\"ATEcombined\" estimand can only be used with \"method == \"cb\"\"")
    }
  }
  if (method %in% c("dbw", "cb", "eb", "mle") == 0) {
    stop("\"method\" must be \"dbw\", \"cb\", \"eb\", or \"mle\"")
  }
  if (method_y %in% c("wls", "logit", "gam", "gambinom") == 0) {
    stop("\"method_y\" must be \"wls\", \"logit\", \"gam\", or \"gambinom\"")
  }
  if (method_y %in% c("wls", "logit") == 0) {
    warning("variance is estimated only when \"method_y\" is \"wls\" or \"logit\"")
    vcov <- FALSE
  }
  if (is.logical(normalize) == FALSE) {
    stop("\"normalize\" must be logical")
  }
  if (normalize == FALSE & method == "dbw") {
    warning("estimated weights are not normalized")
  }
  if (normalize == TRUE & method == "mle") {
    warning("estimated weights are not normalized when \"method == \"mle\"\"")
  }
  if (is.logical(vcov) == FALSE) {
    stop("\"vcov\" must be logical")
  }
  if (length(lambda) != 1) {
    stop("length of \"lambda\" must be 1")
  }
  if (min(clevel) <= 0 | max(clevel) >= 1) {
    stop("\"clevel\" must be between 0 and 1")
  }
  if (method == "dbw" & tol > 1e-6) {
    warning("\"tol\" may be too large. The algorithm might not converge even when judged to do.")
  }
  call <- match.call()
  data <- data.frame(data)
  formula_y <- stats::as.formula(formula_y)
  formula_ps <- stats::as.formula(formula_ps)
  if (method_y %in% c("gam", "gambinom")) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop(
        "package \"mgcv\" must be installed to use \"gam\" or \"gambinom\".",
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
    data2 <- data_ps2[complete_model_y, ]
  }
  y <- c(stats::model.extract(model_y, "response"))
  attr(y, which = "names") <- NULL
  z <- as.matrix(stats::model.matrix(formula_y, model_y))
  if (method_y %in% c("gam", "gambinom")) {
    formula_y <- formula_y0
  }
  model_ps <- stats::model.frame(formula_ps, data = data,
                                 na.action = stats::na.pass)
  response <- c(stats::model.extract(model_ps, "response"))
  attr(response, which = "names") <- NULL
  x_ps <- as.matrix(stats::model.matrix(formula_ps, model_ps))
  response_name <- attr(model_ps, "names")[1]
  z_response <- which(colnames(z) == response_name)
  if(length(z_response) > 0) {
    stop("\"formula_y\" must not include the treatment variable")
  }
  N <- nrow(data)
  if (lambda > 0 & abs(1 - max(apply(x_ps, 2, stats::sd))) > 1e-3) {
    warning("use std_comp() function before dbw() function when \"lambda\" > 0")
  }
  if (setequal(response, c(0, 1)) == FALSE) {
    stop("treatment (response) variable must be binary (0, 1)")
  }
  if (is.null(weights) == 1) {
    weights <- rep(x = 1, times = N)
  } else {
    weights <- (weights[complete_model_ps])[complete_model_y]
  }
  if (length(weights) != N) {
    stop("length of \"weights\" must be the same as the number of rows of \"data\"")
  }
  # Estimation
  result <- dbw_internal(estimand = estimand,
                         data = data,
                         lambda = lambda,
                         method = method,
                         method_y = method_y,
                         response = response,
                         x_ps = x_ps,
                         outcome = y,
                         z = z,
                         weights = weights,
                         vcov = vcov,
                         formula_y = formula_y,
                         tol = tol,
                         init_lambda = init_lambda,
                         normalize = normalize)
  if (sum(result$converged == FALSE) > 0) {
    warning("estimation algorithm did not converge. Estimates are not reliable")
  }
  if (vcov == TRUE) {
    cilength <- sqrt(diag(result$varcov)[nrow(result$varcov)]) *
                 stats::qnorm((1 + clevel) / 2)
    result$ci <- cbind(lower = result$est - cilength,
                       upper = result$est + cilength)
    rownames(result$ci) <- paste0(clevel * 100, "%")
  } else {
    result$ci <- NULL
  }
  if (estimand == "AO") {
    data <- data2
  }
  result$formula_y <- formula_y
  result$formula_ps <- formula_ps
  result$call <- call
  result$data <- data
  result$normalize <- normalize
  result$lambda <- lambda
  class(result) <- "dbw"
  result
}
