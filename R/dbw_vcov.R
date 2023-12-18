## Distribution balancing weighting variance-covariance matrix estimation function
varcov_dbw <- function (ps, beta_trans, lambda, response, x_ps_trans, weights, normalize) {
  if (normalize == TRUE) {
    n <- length(response)
    n0 <- sum(1 - response)
    sum1 <- sum(1 - ps) / n0
    sum2 <- sum(response / ps - 1)
    H1 <- t(x_ps_trans) %*% (x_ps_trans * (-ps + sum1 * response / ps + sum2 * (2 * ps - 1) * ps) * (1 - ps) * weights) / n
    H2 <- 2 * apply(ps * (1 - ps) / n0 * x_ps_trans, 2, sum) %*% 
              t(apply(response * (1 - ps) / ps * x_ps_trans, 2, sum)) / n
    H <- H1 + H2
    H <- H + diag(lambda) / n
    score <- ((1 - sum1 * response / ps - sum2 * ps / n0) * (1 - ps) * x_ps_trans + 
              beta_trans * lambda) / n
  } else { # normalize == FALSE
    n <- length(response)
    H <- t(x_ps_trans) %*% (-x_ps_trans * (response / ps - ps) * (1 - ps) * weights) / n
    H <- H + diag(lambda) / n
    score <- ((response / ps - 1) * (1 - ps) * x_ps_trans + 
              beta_trans * lambda) / n
  }
  list(H = H, score = score)
}

## Covariate balancing weighting variance-covariance matrix estimation function
varcov_cbw <- function (ps, beta_trans, lambda, response, x_ps_trans, weights) {
  n <- length(response)
  H <- t(x_ps_trans) %*% (-x_ps_trans * response * (1 - ps) / ps * weights) / n
  H <- H + diag(lambda) / n
  score <- ((response / ps - 1) * x_ps_trans + 
            beta_trans * lambda) / n
  list(H = H, score = score)
}

## Covariate balancing weighting variance-covariance matrix estimation function for ATEcombined
varcov_cbw_atecombined <- function (beta_trans, lambda, response, x_ps_trans, weights) {
  ps <- logistic(x_ps_trans %*% beta_trans)
  n <- length(response)
  H <- t(x_ps_trans) %*% (-x_ps_trans * (response - ps)^2 / (ps * (1 - ps)) * weights) / n
  H <- H + diag(lambda) / n
  score <- ((response / ps - (1 - response) / (1 - ps)) * x_ps_trans + 
            beta_trans * lambda) / n
  list(H = H, score = score)
}

## Entropy balancing weighting variance-covariance matrix estimation function
varcov_ebw <- function (ps, beta_trans, lambda, response, x_ps_trans, weights) {
  n <- length(response)
  H <- t(x_ps_trans) %*% (-x_ps_trans * response / ps * weights) / n
  H <- H + diag(lambda) / n
  score <- ((response / ps - 1) * x_ps_trans + 
            beta_trans * lambda) / n
  list(H = H, score = score)
}

## MLE inverse probability weighting variance-covariance matrix estimation function
varcov_mlew <- function (ps, beta_trans, lambda, response, x_ps_trans, weights) {
  n <- length(response)
  H <- t(x_ps_trans) %*% (-x_ps_trans * ps * (1 - ps) * weights) / n
  H <- H + diag(lambda) / n
  score <- ((response - ps) * x_ps_trans + 
            beta_trans * lambda) / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the outcome regression for the PS coefficients
varcov_beta <- function (ps, response, x_ps_trans, weights, resid_y, z_trans) {
  n <- length(response)
  H <- t(z_trans) %*% (x_ps_trans * response * resid_y * (1 - ps) / ps * weights) / n
  score <- -response * resid_y / ps * z_trans / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the outcome regression for the PS coefficients for the entropy balancing
varcov_beta_eb <- function (ps, response, x_ps_trans, weights, resid_y, z_trans) {
  n <- length(response)
  H <- t(z_trans) %*% (x_ps_trans * response * resid_y / ps * weights) / n
  score <- -response * resid_y / ps * z_trans / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the outcome regression for the outcome regression coefficients
varcov_wls_gamma <- function (ps, response, weights, resid_y, z_trans) {
  n <- length(response)
  H <- t(z_trans) %*% (z_trans * response / ps * weights) / n
  score <- -response * resid_y / ps * z_trans / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the outcome regression for the binary outcome regression coefficients
varcov_logit_gamma <- function (ps, response, weights, y_hat, resid_y, z_trans) {
  n <- length(response)
  H <- t(z_trans) %*% (z_trans * response * y_hat * (1 - y_hat) / ps * weights) / n
  score <- -response * resid_y / ps * z_trans / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the AO estimation for the PS coefficients
varcov_ao_beta <- function (ps, response, x_ps_trans, weights, y_hat, resid_y, est) {
  n <- length(response)
  H <- apply(-x_ps_trans * response * resid_y * (1 - ps) / ps * weights / n, 2, sum)
  score <- (response * resid_y / ps + y_hat - est) / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the AO estimation for the PS coefficients for the entropy balancing
varcov_ao_beta_eb <- function (ps, response, x_ps_trans, weights, y_hat, resid_y, est) {
  n <- length(response)
  H <- apply(-x_ps_trans * response * resid_y / ps * weights / n, 2, sum)
  score <- (response * resid_y / ps + y_hat - est) / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the AO estimation for the outcome regression coefficients
varcov_ao_wls_gamma <- function (ps, response, weights, y_hat, resid_y, z_trans, est) {
  n <- length(response)
  H <- apply(z_trans * (1 - response / ps) * weights / n, 2, sum)
  score <- (response * resid_y / ps + y_hat - est) / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the AO estimation for the binary outcome regression coefficients
varcov_ao_logit_gamma <- function (ps, response, weights, y_hat, resid_y, z_trans, est) {
  n <- length(response)
  H <- apply(z_trans * (1 - response / ps) * y_hat * (1 - y_hat) * weights / n, 2, sum)
  score <- (response * resid_y / ps + y_hat - est) / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the ATC estimation for the PS coefficients
varcov_atc_beta <- function (ps, response, x_ps_trans, weights, resid_y, est) {
  n <- length(response)
  n0 <- sum(1 - response)
  H <- apply(-x_ps_trans * response * resid_y * (1 - ps) / ps * weights / n0, 2, sum)
  score <- (response * resid_y * (1 - ps) / ps - (1 - response) * resid_y) / n0 - est / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the ATC estimation for the outcome regression coefficients
varcov_atc_wls_gamma <- function (ps, response, weights, resid_y, z_trans, est) {
  n <- length(response)
  n0 <- sum(1 - response)
  H <- apply(-z_trans * ((1 - response) - response * (1 - ps) / ps) * weights / n0, 2, sum)
  score <- (response * resid_y * (1 - ps) / ps - (1 - response) * resid_y) / n0 - est / n
  list(H = H, score = score)
}

## Variance-covariance matrix of the ATC estimation for the binary outcome regression coefficients
varcov_atc_logit_gamma <- function (ps, response, weights, y_hat, resid_y, z_trans, est) {
  n <- length(response)
  n0 <- sum(1 - response)
  H <- apply(-z_trans * ((1 - response) - response * (1 - ps) / ps) * y_hat * (1 - y_hat) * weights / n0, 2, sum)
  score <- (response * resid_y * (1 - ps) / ps - (1 - response) * resid_y) / n0 - est / n
  list(H = H, score = score)
}
