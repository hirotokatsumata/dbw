#' Summarize and plot covariate balance
#'
#' Summarizes and plots covariate balance between treatment and control groups 
#'   before and after the distribution balancing weighting.
#'
#' Estimated weights are normalized to sum to one. Position of the legend is 
#' determined internally. Set \code{absolute = TRUE} and \code{sort = TRUE} for choosing 
#' a better location automatically.
#'
#' @export
#'
#' @param x an object of class “dbw”, usually, a result of a call to \code{\link{dbw}}.
#' @param addcov a one-sided formula specifying additional covariates whose 
#'   balance is checked. Covariates containing NAs are automatically dropped.
#' @param standardize a logical value indicating whether weighted mean 
#'   differences are standardized or not.
#' @param plot a logical value indicating whether a covariate balance plot is
#'   displayed. 
#' @param absolute a logical value indicating whether the absolute values of 
#'   differences in weighted means are used in the covariate balance plot.
#' @param threshold an optional numeric vector used as threshold markers in the 
#'   covariate balance plot.
#' @param sort a logical value indicating whether covariates in the covariate 
#'   balance plot are sorted by the values of differences in the weighted means 
#'   before the distribution balancing weighting.
#' @param ... other graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @return A matrix whose rows are the covariates and columns are the 
#'   differences in the (un)standardized weighted mean between the treatment and
#'   control groups before (\code{diff.un}) and after (\code{diff.adj}) the 
#'   distribution balancing weighting. The standardized weighted mean is the 
#'   weighted mean divided by the standard deviation of the covariate for the 
#'   target population (the treatment group for the average treatment effects on
#'   the treated estimation, the control group for the average treatment effects
#'   on the controlled estimation and the whole population for the other 
#'   quantity of interest). The differences in the categorical variables are not
#'   standardized.
#'
#' @author Hiroto Katsumata
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
#' # Set plot margins
#' oldpar <- par(no.readonly = TRUE) # Just for adjusting plot margins
#' par(mar = c(5.1, 5.1, 4.1, 2.1)) # Just for adjusting plot margins
#'
#' # Misspecified propensity score model
#'
#' # Distribution balancing weighting without regularization
#' fitdbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "dbw",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' plot(fitdbwm, addcov = ~ x1 + x2 + x3 + x4, threshold = c(0.1, 0.2))
#' 
#' # Non-absolute values
#' plot(fitdbwm, addcov = ~ x1 + x2 + x3 + x4, absolute = FALSE, 
#'      threshold = c(0.1, 0.2))
#' 
#' # No sorting
#' plot(fitdbwm, addcov = ~ x1 + x2 + x3 + x4, 
#'      threshold = c(0.1, 0.2), sort = FALSE)
#' 
#' 
#' # Covariate balancing weighting function without regularization
#' fitcbwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "cb",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' plot(fitcbwm, addcov = ~ x1 + x2 + x3 + x4, threshold = c(0.1, 0.2))
#' 
#' # Entropy balancing weighting function without regularization
#' fitebwm <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "eb",
#'                method_y = "wls", data = df, normalize = TRUE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' plot(fitebwm, addcov = ~ x1 + x2 + x3 + x4, threshold = c(0.1, 0.2))
#' 
#' # Standard logistic regression
#' fitmlem <- dbw(formula_y = formula_y, formula_ps = formula_ps_m,
#'                estimand = "ATE", method = "mle",
#'                method_y = "wls", data = df, normalize = FALSE,
#'                vcov = TRUE, lambda = 0, weights = NULL,
#'                clevel = 0.95)
#' plot(fitmlem, addcov = ~ x1 + x2 + x3 + x4, threshold = c(0.1, 0.2))
#'
#' par(oldpar)
#'
#' # Display the covariate balance matrix
#' cb <- plot(fitdbwm, addcov = ~ x1 + x2 + x3 + x4, plot = FALSE)
#' cb
plot.dbw <- function (x, addcov = NULL, standardize = TRUE, plot = TRUE, 
                      absolute = TRUE, threshold = 0, sort = TRUE, ...) {
	data <- x$data
  formula_ps <- stats::as.formula(x$formula_ps)
  model_ps <- stats::model.frame(formula_ps, data = data)
  response <- c(stats::model.extract(model_ps, "response"))
  x_ps <- as.matrix(stats::model.matrix(formula_ps, model_ps))
  formula_y <- stats::as.formula(x$formula_y)
  model_y <- stats::model.frame(formula_y, data = data)
  z <- as.matrix(stats::model.matrix(formula_y, model_y))
  names_z <- setdiff(x = colnames(z), y = colnames(x_ps))
  z <- z[, names_z]
  cov <- as.matrix(cbind(x_ps, z))
  colnames(cov) <- c(colnames(x_ps), names_z)
  original_weights <- x$original_weights
  est_weights <- x$est_weights
  if (is.null(addcov) == 0) {
  	attr(data, which = "na.action") <- stats::na.pass
    formula2 <- stats::as.formula(addcov)
    model2 <- stats::model.frame(formula2, data = data)
    cov2 <- as.matrix(stats::model.matrix(formula2, model2))
    incomplete_cov2 <- which(stats::complete.cases(t(cov2)) == 0)
    if (length(incomplete_cov2) > 0) {
      cov2 <- cov2[, -incomplete_cov2]
      warning("Additional covariates which contain NA values are dropped")
    }
    colnamescov <- colnames(cov)
    colnamescov2 <- colnames(cov2)[-1]
    cov <- as.matrix(cbind(cov, cov2[, -1]))
    colnames(cov) <- c(colnamescov, colnamescov2)
  }
  type <- apply(cov, 2, function(cov) bincont(cov))
  type <- factor(type, levels = c("continuous", "binary"))
  if (x$estimand == "AO") {
    diff.adj <- apply(original_weights * cov, 2, sum) / sum(original_weights) - 
                 apply(est_weights * cov * response, 2, sum) / 
                        sum(est_weights * response)
    diff.un <- apply(original_weights * cov, 2, sum) / sum(original_weights) - 
                apply(original_weights * cov * response, 2, sum) / 
                       sum(original_weights * response)
  } else { # ATT, ATC, ATE, and ATEcombined
    diff.adj <- apply(est_weights * cov * response, 2, sum) / 
                      sum(est_weights * response) - 
                 apply(est_weights * cov * (1 - response), 2, sum) / 
                        sum(est_weights * (1 - response))
    diff.un <- apply(original_weights * cov * response, 2, sum) / 
                      sum(original_weights * response) - 
                apply(original_weights * cov * (1 - response), 2, sum) / 
                       sum(original_weights * (1 - response))
  }
  covariates <- factor(colnames(cov), levels = colnames(cov))
  res <- data.frame(covariates = covariates, 
                    type = type,
                    diff.adj = diff.adj,
                    diff.un = diff.un)
  if (standardize == TRUE) {
    if (x$estimand == "ATT") {
      std <- 
        apply(cov[response == 1, ], 2, 
              function(cov) weighted_sd(cov, 
                                        weights = original_weights[response == 1]))
    } else if (x$estimand == "ATC") {
      std <- 
        apply(cov[response == 0, ], 2, 
              function(cov) weighted_sd(cov, 
                                        weights = original_weights[response == 0]))
    } else { # ATE, ATEcombined and AO
      std <- apply(cov, 2, function(cov) weighted_sd(cov, weights = original_weights))
    }
    std[type == "binary"] <- 1
    std[1] <- 1
    res$diff.adj <- res$diff.adj / std
    res$diff.un <- res$diff.un / std
  }
  rownames(res) <- NULL
  if (plot == TRUE) {
    plot_balance(result = res, standardize = standardize, 
                 absolute = absolute, threshold = threshold, sort = sort)
  }
  invisible(res)
}
