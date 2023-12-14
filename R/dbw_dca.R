## Surrogate function
surrogate <- function (beta, beta_old, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  n0 <- sum(1 - response)
  sum((response / ps + (1 - response) * x_ps %*% beta - sum(response / ps - 1) * log(ps) / n0 -
        ((ps_old - response / ps_old) + 
         sum(1 - ps_old + log(ps_old)) * response * (1 - ps_old) / ps_old / n0 -
         sum(response / ps_old - 1) * (1 - ps_old)^2 / n0) * x_ps %*% beta) * weights) + 
   sum(lambda * beta_trans^2 / 2) * length(response)
}

## First derivative of the surrogate function
u_gradient <- function (beta, beta_old, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  beta_trans_deriv <- deriv_svdtranscoef(svdx = svdx_ps)
  n0 <- sum(1 - response)
  fde1 <- (((1 - response / ps) + 
            sum(log(ps) / n0) * response / ps * (1 - ps) + sum(1 - response / ps) * (1 - ps) / n0) * x_ps) * weights
  fde2 <- (((ps_old - response / ps_old) + 
            sum(1 - ps_old + log(ps_old)) * response * (1 - ps_old) / ps_old / n0 -
            sum(response / ps_old - 1) * (1 - ps_old)^2 / n0) * x_ps) * weights
  fde <- fde1 - fde2
  regularization <- c(t(beta_trans_deriv) %*% (lambda * beta_trans)) * length(response)
  apply(fde, 2, sum) + regularization
}

## DC algorithm for distribution balancing weighting
dca <- function (init, lambda, response, x_ps, weights, svdx_ps, 
                 tol, maxiter = 1000, min_r = 50) {
  ## Initialization
  converged <- 0
  beta <- init
  ps <- logistic(x = x_ps %*% beta)
  ## Start estimation
  for (r in 1:maxiter) {
    beta_old <- beta
    res <- stats::optim(par = beta, 
                        fn = surrogate, 
                        gr = u_gradient,
                        method = "BFGS", 
                        control = list(fnscale = 1, trace = FALSE), 
                        hessian = FALSE,
                        beta_old = beta_old, lambda = lambda, 
                        response = response, x_ps = x_ps, weights = weights, 
                        svdx_ps = svdx_ps)
    beta <- c(res$par)
    beta_trans_dif <- max(abs(svdtranscoef(coef = beta, svdx = svdx_ps) - 
                              svdtranscoef(coef = beta_old, svdx = svdx_ps)))
    if (beta_trans_dif < tol & r > min_r) {
      converged <- r
      break
    }
  }
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  converged <- as.logical(converged != 0)
  list(beta_trans = beta_trans, converged = converged)
}
