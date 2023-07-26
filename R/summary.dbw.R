#' Summarize the distribution balancing weighting estimation results
#'
#' Prints a summary of a fitted \code{dbw} object.
#'
#' Prints a summary of a \code{dbw} object, in a format similar to glm.
#'
#' @export
#'
#' @param object an object of class “dbw”, usually, a result of a call to \code{\link{dbw}}.
#' @param ... additional arguments to be passed to summary.
#'
#' @return 
#' \item{call}{the matched call.}
#' \item{est}{the point estimate of the parameter of interest.}
#' \item{coefficients}{a table including coefficients, standard errors, 
#'   z-values, and two-sided p-values.}
#' \item{effn}{the effective sample size for the parameter of interest
#'   estimation.}
#'
#' @author Hiroto Katsumata
#' 
#' @seealso \code{\link{dbw}}, \code{\link[base]{summary}}
#'
#' @examples # For examples see example(dbw)
summary.dbw <- function (object, ...) {
  if (is.null(object$varcov) == 1) {
    warning("\"varcov\" is empty")
    return(print(object))
  }
  if (object$lambda > 0) {
    warning("When lambda > 0, the uncertainty estimates may not be appropriate\n", 
             "  because the normal approximation may not hold.")
  }
  estimand <- object$estimand
  est <- signif(object$est, digits = getOption("digits") - 3)
  coef_ps <- object$coef_ps
  coef_y <- object$coef_y
  if (estimand == "ATE") {
    est_sub <- signif(object$est_sub, digits = getOption("digits") - 3)
    coef <- c(coef_ps$coef_ps_t, coef_y$coef_y_t, est_sub["mu_t"], 
              coef_ps$coef_ps_c, coef_y$coef_y_c, est_sub["mu_c"])
    k_ps <- length(coef_ps$coef_ps_t)
    k_y <- length(coef_y$coef_y_t)
  } else if (estimand == "ATEcombined") {
    est_sub <- signif(object$est_sub, digits = getOption("digits") - 3)
    coef <- c(coef_ps, coef_y$coef_y_t, coef_y$coef_y_c, 
              est_sub["mu_t"], est_sub["mu_c"])
    k_ps <- length(coef_ps)
    k_y <- length(coef_y$coef_y_t)
  } else { # estimand %in% c("AO", "ATT", "ATC")
    coef <- c(coef_ps, coef_y)
    k_ps <- length(coef_ps)
    k_y <- length(coef_y)
  }
  se <- c(sqrt(diag(object$varcov)))
  pval <- as.vector(2 - 2 * stats::pnorm(abs(c(est, coef) / se)))
  coef.table <- cbind(as.vector(c(est, coef)),
                      as.vector(se),
                      as.vector(c(est, coef) / se),
                      pval)
  colnames(coef.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(coef.table) <- c(names(est), names(coef))
  symp <- stats::symnum(pval, corr = FALSE,
                        cutpoints = c(0, .001, .01, .05, .1, 1),
                        symbols = c("***", "**", "*", ".", " "))
  coef.print <- data.frame(coef.table, as.vector(symp))
  colnames(coef.print) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "")
  cat("\nCall: \n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("\n", rep(x = "#", times = 60), "\n", paste(estimand, "estimate:"), "\n", sep = "")
  print(coef.print[1, ], digits = getOption("digits") - 3)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n", 
      rep(x = "#", times = 60), "\n", sep = "")
  if (estimand == "ATE") {
    cat("\nEstimate of E[Y(1)] and E[Y(0)]:\n")
    print(coef.print[c((2 + k_ps + k_y), (3 + 2 * k_ps + 2 * k_y)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    cat("\nCoefficients for the propensity score model estimation for estimating E[Y(1)]:\n")
    print(coef.print[c(2:(1 + k_ps)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    cat("\nCoefficients for the propensity score model estimation for estimating E[Y(0)]:\n")
    print(coef.print[c((3 + k_ps + k_y):(2 + 2 * k_ps + k_y)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    if (k_y != 0) {
      cat("\nCoefficients for the outcome model estimation for estimating E[Y(1)]:\n")
      print(coef.print[c((2 + k_ps):(1 + k_ps + k_y)), ], 
            digits = getOption("digits") - 3)
      cat("---\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
      cat("\nCoefficients for the outcome model estimation for estimating E[Y(0)]:\n")
      print(coef.print[c((3 + 2 * k_ps + k_y):(2 + 2 * k_ps + 2 * k_y)), ], 
            digits = getOption("digits") - 3)
      cat("---\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    }
  } else if (estimand == "ATEcombined") {
    cat("\nEstimate of E[Y(1)] and E[Y(0)]:\n")
    print(coef.print[c((2 + k_ps + 2 * k_y), (3 + k_ps + 2 * k_y)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    cat("\nCoefficients for the propensity score model estimation:\n")
    print(coef.print[c(2:(1 + k_ps)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    cat("\nCoefficients for the outcome model estimation for estimating E[Y(1)]:\n")
    print(coef.print[c((2 + k_ps):(1 + k_ps + k_y)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    cat("\nCoefficients for the outcome model estimation for estimating E[Y(0)]:\n")
    print(coef.print[c((2 + k_ps + k_y):(1 + k_ps + 2 * k_y)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  } else { # estimand %in% c("AO", "ATT", "ATC")
    cat("\nCoefficients for the propensity score model estimation:\n")
    print(coef.print[c(2:(1 + k_ps)), ], 
          digits = getOption("digits") - 3)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    if (k_y != 0) {
      cat("\nCoefficients for the outcome model estimation:\n")
      print(coef.print[c((2 + k_ps):(1 + k_ps + k_y)), ], 
            digits = getOption("digits") - 3)
      cat("---\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    }
  }
  if (estimand == "AO") {
    cat("\nEffective sample size for the", estimand, "estimation: ", 
        round(object$effn, digits = 2), "\n")
  } else if (estimand == "ATT") {
    cat("\nEffective sample size for the", estimand, "estimation:\n")
    cat("  For estimating E[Y(1)|D_i=1]:", round(object$effn[1], digits = 2), 
        "\n  For estimating E[Y(0)|D_i=1]:  ", round(object$effn[2], digits = 2), "\n")
  } else if (estimand == "ATC") {
    cat("\nEffective sample size for the", estimand, "estimation:\n")
    cat("  For estimating E[Y(1)|D_i=0]:", round(object$effn[1], digits = 2), 
        "\n  For estimating E[Y(0)|D_i=0]:  ", round(object$effn[2], digits = 2), "\n")
  } else { # object$estimand %in% c("ATE", "ATEcombined")
    cat("\nEffective sample size for the", estimand, "estimation:\n")
    cat("  For estimating E[Y(1)]:", round(object$effn[1], digits = 2), 
        "\n  For estimating E[Y(0)]:  ", round(object$effn[2], digits = 2), "\n")
  }
  out <- list("call" = object$call, 
              "est" = est,
              "coefficients" = coef.table,
              "effn" = object$effn)
  invisible(out)
}
