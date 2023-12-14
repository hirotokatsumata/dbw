## Covariate balancing weighting loss function
### Note: this is a loss function to be MINIMIZED
cbloss <- function (beta, estimand, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  if (estimand == "ATEcombined") {
    sum((response / ps + (1 - response) / (1 - ps) + (1 - 2 * response) * x_ps %*% beta) * weights) + 
     sum(lambda * beta_trans^2 / 2) * length(response)
  } else {
    sum((response / ps + (1 - response) * x_ps %*% beta) * weights) + 
     sum(lambda * beta_trans^2 / 2) * length(response)
  }
}

## First derivative of the covariate balancing weighting loss function
cblossgradient <- function (beta, estimand, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  if (estimand == "ATEcombined") {
    fde <- -(response / ps - (1 - response) / (1 - ps)) * x_ps * weights
  } else {
    fde <- (-(response - ps) / ps) * x_ps * weights
  }
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  beta_trans_deriv <- deriv_svdtranscoef(svdx = svdx_ps)
  regularization <- c(t(beta_trans_deriv) %*% (lambda * beta_trans)) * length(response)
  apply(fde, 2, sum) + regularization
}


## Entropy balancing weighting loss function
### Note: this is a loss function to be MINIMIZED
ebloss <- function (beta, lambda, response, x_ps, weights, svdx_ps) {
  ps <- c(exp(x_ps %*% beta))
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  sum((response / ps + x_ps %*% beta) * weights) + 
   sum(lambda * beta_trans^2 / 2) * length(response)
}

## First derivative of the entropy balancing weighting loss function
eblossgradient <- function (beta, lambda, response, x_ps, weights, svdx_ps) {
  ps <- c(exp(x_ps %*% beta))
  fde <- -(response / ps - 1) * x_ps * weights
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  beta_trans_deriv <- deriv_svdtranscoef(svdx = svdx_ps)
  regularization <- c(t(beta_trans_deriv) %*% (lambda * beta_trans)) * length(response)
  apply(fde, 2, sum) + regularization
}

## MLE loss function
### Note: this is a loss function to be MINIMIZED
mleloss <- function (beta, estimand, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  -sum((response * log(ps) + (1 - response) * log(1 - ps)) * weights) + 
    sum(lambda * beta_trans^2 / 2) * length(response)
}

## First derivative of the MLE loss function
mlelossgradient <- function (beta, estimand, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  fde <- -(response - ps) * x_ps * weights
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  beta_trans_deriv <- deriv_svdtranscoef(svdx = svdx_ps)
  regularization <- c(t(beta_trans_deriv) %*% (lambda * beta_trans)) * length(response)
  apply(fde, 2, sum) + regularization
}
