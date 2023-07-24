#' Standardize covariates for propensity score estimation
#'
#' Standardizes covariates for propensity score estimation before the 
#' distribution balancing weighting \code{dbw} when using the regularization.
#'
#' @export
#'
#' @param x a data frame (or one that can be coerced to that class) 
#'   containing the covariates for the propensity score estimation.
#'
#' @return A data frame of the standardized covariates, where the intercept term
#'   is coerced to be a vector of 1s if included.
#'
#' @author Hiroto Katsumata
#' 
#' @seealso \code{\link{dbw}}
#'
#' @examples # For examples see example(dbw)
std <- function (x) {
	x <- data.frame(x)
	mu <- apply(x, 2, mean)
	sigma <- apply(x, 2, stats::sd)
	intercept <- which(sigma == 0)
	mumat <- matrix(rep(mu, each = nrow(x)), nrow = nrow(x), ncol = ncol(x))
	sigmamat <- matrix(rep(sigma, each = nrow(x)), nrow = nrow(x), ncol = ncol(x))
	stdx <- (x - mumat) / sigmamat
	stdx[, intercept] <- 1
	stdx
}
