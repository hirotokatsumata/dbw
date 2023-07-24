#' Print the distribution balancing weighting estimation results
#'
#' Prints a fitted \code{dbw} object.
#'
#' @export
#'
#' @param x an object of class “dbw”, usually, a result of a call to \code{\link{dbw}}.
#' @param ... additional arguments to be passed to print.
#'
#' @return No retrun value, called for side effects.
#'
#' @author Hiroto Katsumata
#' 
#' @seealso \code{\link{dbw}}, \code{\link[base]{print}}
print.dbw <- function (x, ...) {
  est <- signif(x$est, digits = getOption("digits") - 3)
  cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("\n", paste(x$estimand, "estimate: ", est), "\n", sep = "")
  if (x$estimand %in% c("ATE", "ATEcombined")) {
    cat("\n", paste("Estimate of E[Y(1)]: ", 
        signif(x$est_sub[1], digits = getOption("digits") - 3)), "\n", sep = "")
    cat(paste("Estimate of E[Y(0)]: ", 
        signif(x$est_sub[2], digits = getOption("digits") - 3)), "\n", sep = "")
  }
  if (x$estimand == "ATE") {
    cat("\nCoefficients for the propensity score model estimation for estimating E[Y(1)]:\n")
    print(x$coef_ps$coef_ps_t, digits = getOption("digits") - 3)
    cat("\nCoefficients for the propensity score model estimation for estimating E[Y(0)]:\n")
    print(x$coef_ps$coef_ps_c, digits = getOption("digits") - 3)
    cat("\nCoefficients for the outcome model estimation for estimating E[Y(1)]:\n")
    print(x$coef_y$coef_y_t, digits = getOption("digits") - 3)
    cat("\nCoefficients for the outcome model estimation for estimating E[Y(1)]:\n")
    print(x$coef_y$coef_y_c, digits = getOption("digits") - 3)
  } else if (x$estimand == "ATEcombined") {
    cat("\nCoefficients for the propensity score model estimation:\n")
    print(x$coef_ps, digits = getOption("digits") - 3)
    cat("\nCoefficients for the outcome model estimation for estimating E[Y(1)]:\n")
    print(x$coef_y$coef_y_t, digits = getOption("digits") - 3)
    cat("\nCoefficients for the outcome model estimation for estimating E[Y(1)]:\n")
    print(x$coef_y$coef_y_c, digits = getOption("digits") - 3)
  } else { # x$estimand %in% c("AO", "ATT", "ATC"))
    cat("\nCoefficients for the propensity score model estimation:\n")
    print(x$coef_ps, digits = getOption("digits") - 3)
    cat("\nCoefficients for the outcome model estimation:\n")
    print(x$coef_y, digits = getOption("digits") - 3)
  }
  if (x$estimand == "AO") {
    cat("\nEffective sample size for the", x$estimand, "estimation: ", 
        round(x$effn, digits = 2), "\n")
  } else if (x$estimand == "ATT") {
    cat("\nEffective sample size for the", x$estimand, "estimation:\n")
    cat("  For estimating E[Y(1)|D_i=1]:", round(x$effn[1], digits = 2), 
        "\n  For estimating E[Y(0)|D_i=1]:  ", round(x$effn[2], digits = 2), "\n")
  } else if (x$estimand == "ATC") {
    cat("\nEffective sample size for the", x$estimand, "estimation:\n")
    cat("  For estimating E[Y(1)|D_i=0]:", round(x$effn[1], digits = 2), 
        "\n  For estimating E[Y(0)|D_i=0]:  ", round(x$effn[2], digits = 2), "\n")
  } else { # x$estimand %in% c("ATE", "ATEcombined")
    cat("\nEffective sample size for the", x$estimand, "estimation:\n")
    cat("  For estimating E[Y(1)]:", round(x$effn[1], digits = 2), 
        "\n  For estimating E[Y(0)]:  ", round(x$effn[2], digits = 2), "\n")
  }
}
