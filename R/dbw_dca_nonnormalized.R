## Surrogate function without normalization
surrogatenn <- function (beta, beta_old, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  sum((response / ps + (1 - response) * (log(ps) - log(1 - ps)) + 
       (response - ps_old) * x_ps %*% beta) * weights) + 
   sum(lambda * beta_trans^2 / 2) * length(response)
}

## First derivative of the surrogate function without normalization
u_gradientnn <- function (beta, beta_old, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  beta_trans_deriv <- deriv_svdtranscoef(svdx = svdx_ps)
  fde <- ((response / ps + 1 - 2 * response) * x_ps + response - ps_old) * weights
  regularization <- c(t(beta_trans_deriv) %*% (lambda * beta_trans)) * length(response)
  apply(fde, 2, sum) + regularization
}

## DC algorithm for distribution balancing weighting without normalization
dcann <- function (init, lambda, response, x_ps, weights, svdx_ps, 
                   tol, maxiter = 1000, min_r = 50) {
  ## Initialization
  converged <- 0
  beta <- init
  ps <- logistic(x = x_ps %*% beta)
  ## Start estimation
  for (r in 1:maxiter) {
    beta_old <- beta
    res <- stats::optim(par = beta, 
                        fn = surrogatenn, 
                        gr = u_gradientnn,
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
  x_ps_trans <- x_ps %*% diag(svdx_ps$d) %*% t(svdx_ps$v)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  converged <- as.logical(converged != 0)
  list(beta_trans = beta_trans, converged = converged)
}
