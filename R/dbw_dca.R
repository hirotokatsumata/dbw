## Surrogate function for eta > 1
surrogateA <- function (beta, beta_old, eta, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  sum((eta * (response / ps + (1 - response) * x_ps %*% beta +
        ((response - ps_old) * x_ps %*% (beta - beta_old))) + (1 - eta) * log(ps)) * weights) + 
   sum(lambda * beta_trans^2 / 2) * length(response)
}

## Surrogate function for eta < 1
surrogateB <- function (beta, beta_old, eta, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  sum((eta * (response / ps + (1 - response) * x_ps %*% beta) +
        ((eta * (response - ps_old) + (1 - eta) * (1 - ps_old)) * x_ps %*% (beta - beta_old))) * weights) + 
    sum(lambda * beta_trans^2 / 2) * length(response)
}

## First derivative of the surrogate function for eta > 1
u_gradientA <- function (beta, beta_old, eta, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  beta_trans_deriv <- deriv_svdtranscoef(svdx = svdx_ps)
  fde <- ((-eta * (response / ps - 1) + eta * (response - ps_old) + (1 - eta) * (1 - ps)) * x_ps) * weights
  regularization <- c(lambda * t(beta_trans_deriv) %*% beta_trans) * length(response)
  apply(fde, 2, sum) + regularization
}

## First derivative of the surrogate function for eta < 1
u_gradientB <- function (beta, beta_old, eta, lambda, response, x_ps, weights, svdx_ps) {
  ps <- logistic(x = x_ps %*% beta)
  ps_old <- logistic(x = x_ps %*% beta_old)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  beta_trans_deriv <- deriv_svdtranscoef(svdx = svdx_ps)
  fde <- ((-eta * (response / ps - 1) + eta * (response - ps_old) + (1 - eta) * (1 - ps_old)) * x_ps) * weights
  regularization <- c(lambda * t(beta_trans_deriv) %*% beta_trans) * length(response)
  apply(fde, 2, sum) + regularization
}

## DC algorithm for distribution balancing weighting
dca <- function (init, lambda, response, x_ps, weights, svdx_ps, 
                 maxiter1 = 1000, maxiter2 = 400, tol1 = 1e-7, tol2 = 1e-8) {
  ## Initialization
  converged <- 0
  beta <- init
  beta_old <- beta
  ps <- logistic(x = x_ps %*% beta)
  tol3 <- tol2 * 10
  tol3_old <- tol3
  tol3_old2 <- tol3_old
  ## Start estimation
  for (r in 1:maxiter1) {
    maxiter2_2 <- round(maxiter2 * min(c(0.5 + 10 * r / maxiter1), 1))
    beta_r <- beta
    eta <- sum((1 - ps) + (response / ps - 1) * 0.1) / sum(response * (1 - ps) / ps)
    if (eta > 1) {
      for (t in 1:maxiter2_2) {
        beta_old <- beta
        res <- stats::optim(par = beta, 
                            fn = surrogateA, 
                            gr = u_gradientA,
                            method = "BFGS", 
                            control = list(fnscale = 1, trace = FALSE), 
                            hessian = FALSE,
                            beta_old = beta_old, eta = eta, lambda = lambda, 
                            response = response, x_ps = x_ps, weights = weights, 
                            svdx_ps = svdx_ps)
        beta <- c(res$par)
        beta_trans_dif <- max(abs(svdtranscoef(coef = beta, svdx = svdx_ps) - 
                                    svdtranscoef(coef = beta_old, svdx = svdx_ps)))
        if (beta_trans_dif < tol3 & t > 40) {
          break
        }
      }
      beta_trans_dif_r <- max(abs(svdtranscoef(coef = beta, svdx = svdx_ps) - 
                                    svdtranscoef(coef = beta_r, svdx = svdx_ps)))
      tol1_ratio <- beta_trans_dif_r / tol1
      tol3_temp <- sqrt(sqrt(tol3_old2 * tol3_old) * sqrt(tol3 * max(beta_trans_dif_r * 1e-3, tol2)))
      tol3_old2 <- tol3_old
      tol3_old <- tol3
      tol3 <- tol3_temp
      ps <- logistic(x = x_ps %*% beta)
      psavoidzero <- ps
      psavoidzero[response == 0] <- 1
      if (tol1_ratio < 1 & 
          abs(1 - sum(response / psavoidzero) / length(response)) < 1e-4 &
          r > 50) {
        converged <- r
        break
      }
    } else {  # eta <= 1
      for (t in 1:maxiter2_2) {
        beta_old <- beta
        res <- stats::optim(par = beta, 
                            fn = surrogateB, 
                            gr = u_gradientB,
                            method = "BFGS", 
                            control = list(fnscale = 1, trace = FALSE), 
                            hessian = FALSE,
                            beta_old = beta_old, eta = eta, lambda = lambda, 
                            response = response, x_ps = x_ps, weights = weights, 
                            svdx_ps = svdx_ps)
        beta <- c(res$par)
        beta_trans_dif <- max(abs(svdtranscoef(coef = beta, svdx = svdx_ps) - 
                                    svdtranscoef(coef = beta_old, svdx = svdx_ps)))
        if (beta_trans_dif < tol3 & t > 40) {
          break
        }
      }
      beta_trans_dif_r <- max(abs(svdtranscoef(coef = beta, svdx = svdx_ps) - 
                                    svdtranscoef(coef = beta_r, svdx = svdx_ps)))
      tol1_ratio <- beta_trans_dif_r / tol1
      tol3_temp <- sqrt(sqrt(tol3_old2 * tol3_old) * sqrt(tol3 * max(beta_trans_dif_r * 1e-3, tol2)))
      tol3_old2 <- tol3_old
      tol3_old <- tol3
      tol3 <- tol3_temp
      ps <- logistic(x = x_ps %*% beta)
      psavoidzero <- ps
      psavoidzero[response == 0] <- 1
      if (tol1_ratio < 1 & 
          abs(1 - sum(response / psavoidzero) / length(response)) < 1e-4 &
          r > 50) {
        converged <- r
        break
      }
    }
  }
  x_ps_trans <- x_ps %*% diag(svdx_ps$d) %*% t(svdx_ps$v)
  beta_trans <- svdtranscoef(coef = beta, svdx = svdx_ps)
  converged <- as.logical(converged < 1000)
  list(beta_trans = beta_trans, converged = converged)
}
