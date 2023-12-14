## Transform results based on sigular value decomposition: coefficients
svdtranscoef <- function (coef, svdx) {
  c(solve(diag(svdx$d) %*% t(svdx$v)) %*% coef)
}

## Derivatives of the SVDed beta
deriv_svdtranscoef <- function (svdx) {
  solve(diag(svdx$d) %*% t(svdx$v))
}
