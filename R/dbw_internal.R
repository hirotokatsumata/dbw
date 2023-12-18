## Distribution balancing weighting estimation function
dbw0 <- function (initval, lambda, response, x_ps, weights, svdx_ps, tol, init_lambda, normalize) {
  ## Borrowing initial values from cbw
  init <- stats::optim(par = initval, 
                       fn = cbloss, 
                       gr = cblossgradient,
                       method = "BFGS", 
                       control = list(fnscale = 1, trace = FALSE), 
                       hessian = FALSE,
                       estimand = "AO", lambda = lambda * init_lambda, response = response, x_ps = x_ps,
                       weights = weights, svdx_ps = svdx_ps)$par

  ## DC algorithm
  if (normalize == TRUE) {
    ## With normalization
    result <- dca(init = init, lambda = lambda, response = response, x_ps = x_ps, 
                  weights = weights, svdx_ps = svdx_ps, 
                  tol = tol, maxiter = 1000, min_r = 50)
  } else { # normalize == FALSE
    ## Without normalization
    result <- dcann(init = init, lambda = lambda, response = response, x_ps = x_ps, 
                    weights = weights, svdx_ps = svdx_ps, 
                    tol = tol, maxiter = 1000, min_r = 50)
  }
  result
}

## Covariate balancing weighting estimation function
cbw <- function (initval, estimand, lambda, response, x_ps, weights, svdx_ps) {
  result <- stats::optim(par = initval, 
                         fn = cbloss, 
                         gr = cblossgradient,
                         method = "BFGS", 
                         control = list(fnscale = 1, trace = FALSE), 
                         hessian = FALSE,
                         estimand = estimand, lambda = lambda, response = response, x_ps = x_ps,
                         weights = weights, svdx_ps = svdx_ps)

  result$beta_trans <- svdtranscoef(coef = c(result$par), svdx = svdx_ps)
  result$converged <- as.logical(1 - result$convergence)
  result
}

## Entropy balancing weighting estimation function
ebw <- function (initval, lambda, response, x_ps, weights, svdx_ps) {
  result <- stats::optim(par = initval, 
                         fn = ebloss, 
                         gr = eblossgradient,
                         method = "BFGS", 
                         control = list(fnscale = 1, trace = FALSE), 
                         hessian = FALSE,
                         lambda = lambda, response = response, x_ps = x_ps,
                         weights = weights, svdx_ps = svdx_ps)

  result$beta_trans <- svdtranscoef(coef = c(result$par), svdx = svdx_ps)
  result$converged <- as.logical(1 - result$convergence)
  result
}

## Inverse probability weighting with maximum likelihood estimation function
mlew <- function (initval, lambda, response, x_ps, weights, svdx_ps) {
  result <- stats::optim(par = initval, 
                         fn = mleloss, 
                         gr = mlelossgradient,
                         method = "BFGS", 
                         control = list(fnscale = 1, trace = FALSE), 
                         hessian = FALSE,
                         lambda = lambda, response = response, x_ps = x_ps,
                         weights = weights, svdx_ps = svdx_ps)

  result$beta_trans <- svdtranscoef(coef = c(result$par), svdx = svdx_ps)
  result$converged <- as.logical(1 - result$convergence)
  result
}

## Distribution balancing weighting internal function
dbw_internal <- function (estimand, data, lambda, method, method_y, 
                          response, x_ps, outcome, z, weights, 
                          vcov, formula_y, tol, init_lambda, normalize) {
  names_z <- colnames(z)
  names_z[apply(z, 2, stats::sd) == 0] <- "(Intercept)"
  names_x_ps <- colnames(x_ps)
  names_x_ps[apply(x_ps, 2, stats::sd) == 0] <- "(Intercept)"
  intercept <- which(apply(x_ps, 2, stats::sd) == 0)
  x_ps0 <- x_ps
  svdx_ps <- svd(x = x_ps)
  x_ps <- svdx_ps$u
  lambda <- rep(x = lambda, times = ncol(x_ps))
  lambda[intercept] <- 0
  estimand0 <- estimand
  if (estimand == "ATT") {
    estimand <- "ATC"
    response <- 1 - response
  }
  if (estimand == "AO") {
    weights[response == 1] <- weights[response == 1] / mean(weights[response == 1])
  } else { # estimand %in% c("ATE", "ATC", "ATEcombined")
    weights[response == 1] <- weights[response == 1] / mean(weights[response == 1])
    weights[response == 0] <- weights[response == 0] / mean(weights[response == 0])
  }
  initval <- rep(x = 0, times = ncol(x_ps))
  est_sub <- NULL
  ##########
  ## Weight estimation
  if (method == "dbw") {
    ## Distribution balancing weighting
    if (estimand == "ATE") {
      result_w_t <- dbw0(initval = initval, lambda = lambda, response = response, 
                         x_ps = x_ps, weights = weights, svdx_ps = svdx_ps, 
                         tol = tol, init_lambda = init_lambda, normalize = normalize)
      result_w_c <- dbw0(initval = initval, lambda = lambda, response = 1 - response, 
                         x_ps = x_ps, weights = weights, svdx_ps = svdx_ps, 
                         tol = tol, init_lambda = init_lambda, normalize = normalize)
      coef_t <- result_w_t$beta_trans
      coef_c <- -result_w_c$beta_trans
      ps_t <- logistic(x = x_ps0 %*% coef_t)
      ps_c <- logistic(x = x_ps0 %*% coef_c)
      est_weights_t <- invprob(ps = ps_t, response = response, estimand = "AO", weights = weights)
      est_weights_c <- invprob(ps = 1 - ps_c, response = 1 - response, estimand = "AO", weights = weights)
      est_weights <- est_weights_t + est_weights_c
      converged <- c(result_w_t$converged, result_w_c$converged)
    } else if (estimand == "AO") {
      result_w <- dbw0(initval = initval, lambda = lambda, response = response, 
                       x_ps = x_ps, weights = weights, svdx_ps = svdx_ps, 
                       tol = tol, init_lambda = init_lambda, normalize = normalize)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights <- invprob(ps = ps, response = response, estimand = "AO", weights = weights)
      converged <- result_w$converged
    } else { # estimand == ATC: equivalent to cbw
      result_w <- cbw(initval = initval, estimand = "ATC", lambda = lambda, response = response, 
                      x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights <- invprob(ps = ps, response = response, estimand = "ATC", weights = weights)
      converged <- result_w$converged
    }
  } else if (method == "cb") {
    ## Covariate balancing weighting
    if (estimand == "ATE") {
      result_w_t <- cbw(initval = initval, estimand = "AO", lambda = lambda, response = response, 
                        x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      result_w_c <- cbw(initval = initval, estimand = "AO", lambda = lambda, response = 1 - response, 
                        x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef_t <- result_w_t$beta_trans
      coef_c <- -result_w_c$beta_trans
      ps_t <- logistic(x = x_ps0 %*% coef_t)
      ps_c <- logistic(x = x_ps0 %*% coef_c)
      est_weights_t <- invprob(ps = ps_t, response = response, estimand = "AO", weights = weights)
      est_weights_c <- invprob(ps = 1 - ps_c, response = 1 - response, estimand = "AO", weights = weights)
      est_weights <- est_weights_t + est_weights_c
      converged <- c(result_w_t$converged, result_w_c$converged)
    } else if (estimand == "AO") {
      result_w <- cbw(initval = initval, estimand = "AO", lambda = lambda, response = response, 
                      x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights <- invprob(ps = ps, response = response, estimand = "AO", weights = weights)
      converged <- result_w$converged
    } else if (estimand == "ATC") {
      result_w <- cbw(initval = initval, estimand = "ATC", lambda = lambda, response = response, 
                      x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights <- invprob(ps = ps, response = response, estimand = "ATC", weights = weights)
      converged <- result_w$converged
    } else { # estimand == ATEcombined
      result_w <- cbw(initval = initval, estimand = "ATEcombined", lambda = lambda, response = response, 
                      x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights_t <- invprob(ps = ps, response = response, estimand = "AO", weights = weights)
      est_weights_c <- invprob(ps = 1 - ps, response = 1 - response, estimand = "AO", weights = weights)
      est_weights <- est_weights_t + est_weights_c
      converged <- result_w$converged
    }
  } else if (method == "eb") {
    ## Entropy balancing weighting
    if (estimand == "ATE") {
      result_w_t <- ebw(initval = initval, lambda = lambda, response = response, 
                        x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      result_w_c <- ebw(initval = initval, lambda = lambda, response = 1 - response, 
                        x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef_t <- result_w_t$beta_trans
      coef_c <- -result_w_c$beta_trans
      ps_t <- c(exp(x = x_ps0 %*% coef_t))
      ps_c <- c(exp(x = x_ps0 %*% coef_c))
      est_weights_t <- invprob(ps = ps_t, response = response, estimand = "AO", weights = weights)
      est_weights_c <- invprob(ps = 1 / ps_c, response = 1 - response, estimand = "AO", weights = weights)
      est_weights <- est_weights_t + est_weights_c
      converged <- c(result_w_t$converged, result_w_c$converged)
    } else if (estimand == "AO") {
      result_w <- ebw(initval = initval, lambda = lambda, response = response, 
                      x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- c(exp(x = x_ps0 %*% coef))
      est_weights <- invprob(ps = ps, response = response, estimand = "AO", weights = weights)
      converged <- result_w$converged
    } else { # estimand == ATC: equivalent to cbw
      result_w <- cbw(initval = initval, estimand = "ATC", lambda = lambda, response = response, 
                      x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights <- invprob(ps = ps, response = response, estimand = "ATC", weights = weights)
      converged <- result_w$converged
    }
  } else { # method == mle
    ## Maximum likelihood estimation
    if (estimand == "ATE") {
      result_w_t <- mlew(initval = initval, lambda = lambda, response = response, 
                         x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      result_w_c <- mlew(initval = initval, lambda = lambda, response = 1 - response, 
                         x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef_t <- result_w_t$beta_trans
      coef_c <- -result_w_c$beta_trans
      ps_t <- logistic(x = x_ps0 %*% coef_t)
      ps_c <- logistic(x = x_ps0 %*% coef_c)
      est_weights_t <- invprob(ps = ps_t, response = response, estimand = "AO", weights = weights)
      est_weights_c <- invprob(ps = 1 - ps_c, response = 1 - response, estimand = "AO", weights = weights)
      est_weights <- est_weights_t + est_weights_c
      converged <- c(result_w_t$converged, result_w_c$converged)
    } else if (estimand == "AO") {
      result_w <- mlew(initval = initval, lambda = lambda, response = response, 
                       x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights <- invprob(ps = ps, response = response, estimand = "AO", weights = weights)
      est <- sum(est_weights * outcome * response)
      converged <- result_w$converged
    } else { # estimand == "ATC"
      result_w <- mlew(initval = initval, lambda = lambda, response = response, 
                       x_ps = x_ps, weights = weights, svdx_ps = svdx_ps)
      coef <- result_w$beta_trans
      ps <- logistic(x = x_ps0 %*% coef)
      est_weights <- invprob(ps = ps, response = response, estimand = "ATC", weights = weights)
      converged <- result_w$converged
    }
  }
  ## End of weight estimation
  ##########
  ## Outcome model estimation and variance estimation
  if (estimand == "AO") {
    ## Outcome model estimation
    result_y <- dbw_y(formula_y = formula_y,
                      df = data[response == 1, ], 
                      fulldf = data, 
                      method_y = method_y,
                      est_weights = est_weights[response == 1])
    y_hat <- result_y$pred
    resid_y <- rep(x = 0, times = length(response))
    resid_y[response == 1] <- outcome[response == 1] - result_y$pred[response == 1]
    est <- sum(est_weights * resid_y + weights * y_hat / sum(weights))
    gamma <- result_y$gamma
    ## Variance estimation
    if (vcov == FALSE) {
      varcov <- NULL
    } else {
    if (method == "eb") {
      hs1 <- varcov_ebw(ps = ps, beta_trans = coef, lambda = lambda, 
                        response = response, x_ps_trans = x_ps0, 
                        weights = weights)
      hs2 <- varcov_beta_eb(ps = ps, response = response, 
                            x_ps_trans = x_ps0, weights = weights, 
                            resid_y = resid_y, z_trans = z)
      hs4 <- varcov_ao_beta_eb(ps = ps, response = response, 
                               x_ps_trans = x_ps0, weights = weights, 
                               y_hat = y_hat, resid_y = resid_y, 
                               est = est)
    } else if (method == "dbw") {
      hs1 <- varcov_dbw(ps = ps, beta_trans = coef, lambda = lambda, 
                        response = response, x_ps_trans = x_ps0, 
                        weights = weights, normalize = normalize)
    } else if (method == "cb") {
      hs1 <- varcov_cbw(ps = ps, beta_trans = coef, lambda = lambda, 
                        response = response, x_ps_trans = x_ps0, 
                        weights = weights)
    } else { # method == "mlew"
      hs1 <- varcov_mlew(ps = ps, beta_trans = coef, lambda = lambda, 
                         response = response, x_ps_trans = x_ps0, 
                         weights = weights)      
    }
    if (method != "eb") {
      hs2 <- varcov_beta(ps = ps, response = response, 
                         x_ps_trans = x_ps0, weights = weights, 
                         resid_y = resid_y, z_trans = z)
      hs4 <- varcov_ao_beta(ps = ps, response = response, 
                            x_ps_trans = x_ps0, weights = weights, 
                            y_hat = y_hat, resid_y = resid_y, 
                            est = est)
    }
    if (method_y == "wls") {
      hs3 <- varcov_wls_gamma(ps = ps, response = response, 
                              weights = weights, resid_y = resid_y, 
                              z_trans = z)
      hs5 <- varcov_ao_wls_gamma(ps = ps, response = response, 
                                 weights = weights, y_hat = y_hat, 
                                 resid_y = resid_y, z_trans = z, 
                                 est = est)
    } else { # method_y == "logit"
      hs3 <- varcov_logit_gamma(ps = ps, response = response, 
                                weights = weights, y_hat = y_hat, 
                                resid_y = resid_y, z_trans = z)
      hs5 <- varcov_ao_logit_gamma(ps = ps, response = response, 
                                   weights = weights, y_hat = y_hat, 
                                   resid_y = resid_y, z_trans = z, 
                                   est = est)
    }
    if (ncol(z) == 0) {
      H <- rbind(cbind(hs1$H, 0),
                 c(hs4$H, -1))
    } else { 
      zero_matrix <- matrix(0, nrow = ncol(x_ps0), ncol = ncol(z))
      H <- rbind(cbind(hs1$H, zero_matrix, 0),
                 cbind(hs2$H, hs3$H, 0),
                 c(hs4$H, hs5$H, -1))
    }
    H <- as.matrix(H)
    score <- as.matrix(cbind(hs1$score, hs2$score, hs4$score))
    Sigma <- t(score) %*% (score * weights)
    varcov <- solve(H) %*% Sigma %*% solve(t(H))
    }
  } else if (estimand == "ATC") {
    ## Outcome model estimation
    result_y <- dbw_y(formula_y = formula_y,
                      df = data[response == 1, ], 
                      fulldf = data, 
                      method_y = method_y,
                      est_weights = est_weights[response == 1])
    y_hat <- result_y$pred
    resid_y <- outcome - result_y$pred
    est <- sum(response * est_weights * resid_y -
               (1 - response) * weights * resid_y / sum((1 - response) * weights))
    gamma <- result_y$gamma
    ## Variance estimation
    if (vcov == FALSE) {
      varcov <- NULL
    } else {
    if (method == "mlew") {
      hs1 <- varcov_mlew(ps = ps, beta_trans = coef, lambda = lambda, 
                         response = response, x_ps_trans = x_ps0, 
                         weights = weights)      
    } else { # method %in% c("dbw", "cbw", "ebw")
      hs1 <- varcov_cbw(ps = ps, beta_trans = coef, lambda = lambda, 
                        response = response, x_ps_trans = x_ps0, 
                        weights = weights)
    }
    hs2 <- varcov_beta(ps = ps, response = response, 
                       x_ps_trans = x_ps0, weights = weights, 
                       resid_y = resid_y, z_trans = z)
    hs4 <- varcov_atc_beta(ps = ps, response = response, 
                           x_ps_trans = x_ps0, weights = weights, 
                           resid_y = resid_y, est = est)
    if (method_y == "wls") {
      hs3 <- varcov_wls_gamma(ps = ps, response = response, 
                              weights = weights, resid_y = resid_y, 
                              z_trans = z)
      hs5 <- varcov_atc_wls_gamma(ps = ps, response = response, 
                                  weights = weights, resid_y = resid_y, 
                                  z_trans = z, est = est)
    } else { # method_y == "logit"
      hs3 <- varcov_logit_gamma(ps = ps, response = response, 
                                weights = weights, y_hat = y_hat, 
                                resid_y = resid_y, z_trans = z)
      hs5 <- varcov_atc_logit_gamma(ps = ps, response = response, 
                                    weights = weights, y_hat = y_hat, 
                                    resid_y = resid_y, z_trans = z, 
                                    est = est)
    }
    if (ncol(z) == 0) {
      H <- rbind(cbind(hs1$H, 0),
                 c(hs4$H, -1))
    } else { 
      zero_matrix <- matrix(0, nrow = ncol(x_ps0), ncol = ncol(z))
      H <- rbind(cbind(hs1$H, zero_matrix, 0),
                 cbind(hs2$H, hs3$H, 0),
                 c(hs4$H, hs5$H, -1))
    }
    H <- as.matrix(H)
    score <- as.matrix(cbind(hs1$score, hs2$score, hs4$score))
    Sigma <- t(score) %*% (score * weights)
    varcov <- solve(H) %*% Sigma %*% solve(t(H))
    }
  } else if (estimand == "ATE") {
    ## Outcome model estimation
    result_y_t <- dbw_y(formula_y = formula_y,
                        df = data[response == 1, ], 
                        fulldf = data, 
                        method_y = method_y,
                        est_weights = est_weights[response == 1])
    y_hat_t <- result_y_t$pred
    resid_y_t <- rep(x = 0, times = length(response))
    resid_y_t[response == 1] <- outcome[response == 1] - result_y_t$pred[response == 1]
    result_y_c <- dbw_y(formula_y = formula_y,
                        df = data[response == 0, ], 
                        fulldf = data, 
                        method_y = method_y,
                        est_weights = est_weights[response == 0])
    y_hat_c <- result_y_c$pred
    resid_y_c <- rep(x = 0, times = length(response))
    resid_y_c[response == 0] <- outcome[response == 0] - result_y_c$pred[response == 0]
    ao_t <- sum(response * est_weights * resid_y_t + 
                weights * y_hat_t / sum(weights))
    ao_c <- sum((1 - response) * est_weights * resid_y_c + 
                weights * y_hat_c / sum(weights))
    est <- ao_t - ao_c
    est_sub <- c(ao_t, ao_c)
    names(est_sub) <- c("mu_t", "mu_c")
    gamma_t <- result_y_t$gamma
    gamma_c <- result_y_c$gamma
    ## Variance estimation
    if (vcov == FALSE) {
      varcov <- NULL
    } else {
    if (method == "eb") {
      hs11 <- varcov_ebw(ps = ps_t, beta_trans = coef_t, lambda = lambda, 
                         response = response, x_ps_trans = x_ps0, 
                         weights = weights)
      hs21 <- varcov_ebw(ps = 1 / ps_c, beta_trans = coef_c, lambda = lambda, 
                         response = 1 - response, x_ps_trans = x_ps0, 
                         weights = weights)
      hs12 <- varcov_beta_eb(ps = ps_t, response = response, 
                             x_ps_trans = x_ps0, weights = weights, 
                             resid_y = resid_y_t, z_trans = z)
      hs22 <- varcov_beta_eb(ps = 1 / ps_c, response = 1 - response, 
                             x_ps_trans = x_ps0, weights = weights, 
                             resid_y = resid_y_c, z_trans = z)
      hs14 <- varcov_ao_beta_eb(ps = ps_t, response = response, 
                                x_ps_trans = x_ps0, weights = weights, 
                                y_hat = y_hat_t, resid_y = resid_y_t, 
                                est = ao_t)
      hs24 <- varcov_ao_beta_eb(ps = 1 / ps_c, response = 1 - response, 
                                x_ps_trans = x_ps0, weights = weights, 
                                y_hat = y_hat_c, resid_y = resid_y_c, 
                                est = ao_c)
    } else if (method == "dbw") {
      hs11 <- varcov_dbw(ps = ps_t, beta_trans = coef_t, lambda = lambda, 
                         response = response, x_ps_trans = x_ps0, 
                         weights = weights, normalize = normalize)
      hs21 <- varcov_dbw(ps = 1 - ps_c, beta_trans = coef_c, lambda = lambda, 
                         response = 1 - response, x_ps_trans = x_ps0, 
                         weights = weights, normalize = normalize)
    } else if (method == "cb") {
      hs11 <- varcov_cbw(ps = ps_t, beta_trans = coef_t, lambda = lambda, 
                         response = response, x_ps_trans = x_ps0, 
                         weights = weights)
      hs21 <- varcov_cbw(ps = 1 - ps_c, beta_trans = coef_c, lambda = lambda, 
                         response = 1 - response, x_ps_trans = x_ps0, 
                         weights = weights)
    } else { # method == "mlew"
      hs11 <- varcov_mlew(ps = ps_t, beta_trans = coef_t, lambda = lambda, 
                          response = response, x_ps_trans = x_ps0, 
                          weights = weights)      
      hs21 <- varcov_mlew(ps = 1 - ps_c, beta_trans = coef_c, lambda = lambda, 
                          response = 1 - response, x_ps_trans = x_ps0, 
                          weights = weights)      
    }
    if (method != "eb") {
      hs12 <- varcov_beta(ps = ps_t, response = response, 
                          x_ps_trans = x_ps0, weights = weights, 
                          resid_y = resid_y_t, z_trans = z)
      hs22 <- varcov_beta(ps = 1 - ps_c, response = 1 - response, 
                          x_ps_trans = x_ps0, weights = weights, 
                          resid_y = resid_y_c, z_trans = z)
      hs14 <- varcov_ao_beta(ps = ps_t, response = response, 
                             x_ps_trans = x_ps0, weights = weights, 
                             y_hat = y_hat_t, resid_y = resid_y_t, 
                             est = ao_t)
      hs24 <- varcov_ao_beta(ps = 1 - ps_c, response = 1 - response, 
                             x_ps_trans = x_ps0, weights = weights, 
                             y_hat = y_hat_c, resid_y = resid_y_c, 
                             est = ao_c)
    }
    if (method_y == "wls") {
      hs13 <- varcov_wls_gamma(ps = ps_t, response = response, 
                               weights = weights, resid_y = resid_y_t, 
                               z_trans = z)
      hs15 <- varcov_ao_wls_gamma(ps = ps_t, response = response, 
                                  weights = weights, y_hat = y_hat_t, 
                                  resid_y = resid_y_t, z_trans = z, 
                                  est = ao_t)
      if (method == "eb") {
        hs23 <- varcov_wls_gamma(ps = 1 / ps_c, response = 1 - response, 
                                 weights = weights, resid_y = resid_y_c, 
                                  z_trans = z)
        hs25 <- varcov_ao_wls_gamma(ps = 1 / ps_c, response = 1 - response, 
                                    weights = weights, y_hat = y_hat_c, 
                                    resid_y = resid_y_c, z_trans = z, 
                                    est = ao_c)
      } else {
        hs23 <- varcov_wls_gamma(ps = 1 - ps_c, response = 1 - response, 
                                 weights = weights, resid_y = resid_y_c, 
                                  z_trans = z)
        hs25 <- varcov_ao_wls_gamma(ps = 1 - ps_c, response = 1 - response, 
                                    weights = weights, y_hat = y_hat_c, 
                                    resid_y = resid_y_c, z_trans = z, 
                                    est = ao_c)
      }
    } else { # method_y == "logit"
      hs13 <- varcov_logit_gamma(ps = ps_t, response = response, 
                                 weights = weights, y_hat = y_hat_t, 
                                 resid_y = resid_y_t, z_trans = z)
      hs15 <- varcov_ao_logit_gamma(ps = ps_t, response = response, 
                                    weights = weights, y_hat = y_hat_t, 
                                    resid_y = resid_y_t, z_trans = z, 
                                    est = ao_t)
      if (method == "eb") {
        hs23 <- varcov_logit_gamma(ps = 1 / ps_c, response = 1 - response, 
                                   weights = weights, y_hat = y_hat_c, 
                                   resid_y = resid_y_c, z_trans = z)
        hs25 <- varcov_ao_logit_gamma(ps = 1 / ps_c, response = 1 - response, 
                                      weights = weights, y_hat = y_hat_c, 
                                      resid_y = resid_y_c, z_trans = z, 
                                      est = ao_c)
      } else {
        hs23 <- varcov_logit_gamma(ps = 1 - ps_c, response = 1 - response, 
                                   weights = weights, y_hat = y_hat_c, 
                                   resid_y = resid_y_c, z_trans = z)
        hs25 <- varcov_ao_logit_gamma(ps = 1 - ps_c, response = 1 - response, 
                                      weights = weights, y_hat = y_hat_c, 
                                      resid_y = resid_y_c, z_trans = z, 
                                      est = ao_c)
      }
    }
    if (ncol(z) == 0) {
      large_zero <- matrix(0, nrow = ncol(x_ps0) + 1, 
                              ncol = ncol(x_ps0) + 1)
      zero_vec_x <- rep(x = 0, times = ncol(x_ps0))
      H1 <- rbind(cbind(hs11$H, 0),
                  c(hs14$H, -1))
      H2 <- rbind(cbind(hs21$H, 0),
                  c(hs24$H, -1))
      H <- rbind(cbind(H1, large_zero, 0),
                 cbind(large_zero, H2, 0),
                 c(zero_vec_x, 1, zero_vec_x, -1, -1))
    } else {
      zero_matrix <- matrix(0, nrow = ncol(x_ps0), ncol = ncol(z))
      large_zero <- matrix(0, nrow = ncol(x_ps0) + ncol(z) + 1, 
                              ncol = ncol(x_ps0) + ncol(z) + 1)
      zero_vec_x <- rep(x = 0, times = ncol(x_ps0))
      zero_vec_z <- rep(x = 0, times = ncol(z))
      H1 <- rbind(cbind(hs11$H, zero_matrix, 0),
                  cbind(hs12$H, hs13$H, 0),
                  c(hs14$H, hs15$H, -1))
      H2 <- rbind(cbind(hs21$H, zero_matrix, 0),
                  cbind(hs22$H, hs23$H, 0),
                  c(hs24$H, hs25$H, -1))
      H <- rbind(cbind(H1, large_zero, 0),
                 cbind(large_zero, H2, 0),
                 c(zero_vec_x, zero_vec_z, 1, 
                   zero_vec_x, zero_vec_z, -1, -1))
    }
    H <- as.matrix(H)
    score <- as.matrix(cbind(hs11$score, hs12$score, hs14$score,
                             hs21$score, hs22$score, hs24$score,
                             0))
    Sigma <- t(score) %*% (score * weights)
    varcov <- solve(H) %*% Sigma %*% solve(t(H))
    }
  } else { # estimand == "ATEcombined"
    ## Outcome model estimation
    result_y_t <- dbw_y(formula_y = formula_y,
                        df = data[response == 1, ], 
                        fulldf = data, 
                        method_y = method_y,
                        est_weights = est_weights[response == 1])
    y_hat_t <- result_y_t$pred
    resid_y_t <- rep(x = 0, times = length(response))
    resid_y_t[response == 1] <- outcome[response == 1] - result_y_t$pred[response == 1]
    result_y_c <- dbw_y(formula_y = formula_y,
                        df = data[response == 0, ], 
                        fulldf = data, 
                        method_y = method_y,
                        est_weights = est_weights[response == 0])
    y_hat_c <- result_y_c$pred
    resid_y_c <- rep(x = 0, times = length(response))
    resid_y_c[response == 0] <- outcome[response == 0] - result_y_c$pred[response == 0]
    ao_t <- sum(response * est_weights * resid_y_t + 
                weights * y_hat_t / sum(weights))
    ao_c <- sum((1 - response) * est_weights * resid_y_c + 
                weights * y_hat_c / sum(weights))
    est <- ao_t - ao_c
    est_sub <- c(ao_t, ao_c)
    names(est_sub) <- c("mu_t", "mu_c")
    gamma_t <- result_y_t$gamma
    gamma_c <- result_y_c$gamma
    ## Variance estimation
    if (vcov == FALSE) {
      varcov <- NULL
    } else {
    hs11 <- varcov_cbw(ps = ps, beta_trans = coef, lambda = lambda, 
                       response = response, x_ps_trans = x_ps0, 
                       weights = weights)
    hs12 <- varcov_beta(ps = ps, response = response, 
                        x_ps_trans = x_ps0, weights = weights, 
                        resid_y = resid_y_t, z_trans = z)
    hs22 <- varcov_beta(ps = 1 - ps, response = 1 - response, 
                        x_ps_trans = x_ps0, weights = weights, 
                        resid_y = resid_y_c, z_trans = z)
    hs14 <- varcov_ao_beta(ps = ps, response = response, 
                           x_ps_trans = x_ps0, weights = weights, 
                           y_hat = y_hat_t, resid_y = resid_y_t, 
                           est = ao_t)
    hs24 <- varcov_ao_beta(ps = 1 - ps, response = 1 - response, 
                           x_ps_trans = x_ps0, weights = weights, 
                           y_hat = y_hat_c, resid_y = resid_y_c, 
                           est = ao_c)
    if (method_y == "wls") {
      hs13 <- varcov_wls_gamma(ps = ps, response = response, 
                               weights = weights, resid_y = resid_y_t, 
                               z_trans = z)
      hs15 <- varcov_ao_wls_gamma(ps = ps, response = response, 
                                  weights = weights, y_hat = y_hat_t, 
                                  resid_y = resid_y_t, z_trans = z, 
                                  est = ao_t)
      hs23 <- varcov_wls_gamma(ps = 1 - ps, response = 1 - response, 
                               weights = weights, resid_y = resid_y_c, 
                               z_trans = z)
      hs25 <- varcov_ao_wls_gamma(ps = 1 - ps, response = 1 - response, 
                                  weights = weights, y_hat = y_hat_c, 
                                  resid_y = resid_y_c, z_trans = z, 
                                  est = ao_c)
    } else { # method_y == "logit"
      hs13 <- varcov_logit_gamma(ps = ps, response = response, 
                                 weights = weights, y_hat = y_hat_t, 
                                 resid_y = resid_y_t, z_trans = z)
      hs15 <- varcov_ao_logit_gamma(ps = ps, response = response, 
                                    weights = weights, y_hat = y_hat_t, 
                                    resid_y = resid_y_t, z_trans = z, 
                                    est = ao_t)
      hs23 <- varcov_logit_gamma(ps = 1 - ps, response = 1 - response, 
                                 weights = weights, y_hat = y_hat_c, 
                                 resid_y = resid_y_c, z_trans = z)
      hs25 <- varcov_ao_logit_gamma(ps = 1 - ps, response = 1 - response, 
                                    weights = weights, y_hat = y_hat_c, 
                                    resid_y = resid_y_c, z_trans = z, 
                                    est = ao_c)
    }
    if (ncol(z) == 0) {
      zero_vec_x <- rep(x = 0, times = ncol(x_ps0))
      H <- rbind(cbind(hs11$H, 0, 0, 0),
                 c(hs14$H, -1, 0, 0), 
                 c(hs24$H, 0, -1, 0), 
                 c(zero_vec_x, 1, -1, -1))
    } else {
      zero_matrix <- matrix(0, nrow = ncol(x_ps0), ncol = ncol(z))
      zero_matrix2 <- matrix(0, nrow = ncol(z), ncol = ncol(z))
      zero_vec_x <- rep(x = 0, times = ncol(x_ps0))
      zero_vec_z <- rep(x = 0, times = ncol(z))
      H <- rbind(cbind(hs11$H, zero_matrix, zero_matrix, 0, 0, 0),
                 cbind(hs12$H, hs13$H, zero_matrix2, 0, 0, 0),
                 cbind(hs22$H, zero_matrix2, hs23$H, 0, 0, 0),
                 c(hs14$H, hs15$H, zero_vec_z, -1, 0, 0), 
                 c(hs24$H, zero_vec_z, hs25$H, 0, -1, 0), 
                 c(zero_vec_x, zero_vec_z, zero_vec_z, 1, -1, -1))
    }
    H <- as.matrix(H)
    score <- as.matrix(cbind(hs11$score, hs12$score, hs22$score,
                             hs14$score, hs24$score, 
                             0))
    Sigma <- t(score) %*% (score * weights)
    varcov <- solve(H) %*% Sigma %*% solve(t(H))
    }
  }
  ## Move varcov for QoI to the top row and column
  if (vcov == TRUE) {
    varcov1 <- varcov[nrow(varcov), -ncol(varcov)]
    varcov2 <- varcov[-nrow(varcov), ncol(varcov)]
    varcov <- rbind(c(varcov[nrow(varcov), ncol(varcov)], varcov1),
                    cbind(varcov2, varcov[-nrow(varcov), -ncol(varcov)]))
    varcov <- as.matrix(varcov)
  }
  ## End of outcome model estimation and variance estimation
  ##########
  ## Effective sample size
  if (estimand != "AO") {
    effn <- c(sum(est_weights[response == 1])^2 / 
               sum(est_weights[response == 1]^2),
              sum(est_weights[response == 0])^2 / 
               sum(est_weights[response == 0]^2))
  } else { # estimand == "AO"
    effn <- c(sum(est_weights[response == 1])^2 / 
              sum(est_weights[response == 1]^2))
  }
  effn_original <- sum(weights)^2 / sum((weights)^2)
  ## Names
  names(est) <- estimand0
  if (estimand == "ATEcombined") {
    if (ncol(z) != 0 & method_y %in% c("wls", "logit")) {
      names_gamma_t <- paste0("outcome_t_", names_z)
      names_gamma_c <- paste0("outcome_c_", names_z)
      names(gamma_t) <- names_gamma_t
      names(gamma_c) <- names_gamma_c
    } else {
      names_gamma_t <- NULL
      names_gamma_c <- NULL
    }
    names(coef) <- paste0("ps_", names_x_ps)
    gamma <- list(coef_y_t = gamma_t, coef_y_c = gamma_c)
    y_hat <- list(predicted_y_t = y_hat_t, predicted_y_c = y_hat_c)
    if (vcov == TRUE) {
      names_varcov <- c(estimand0, names(coef), names_gamma_t, 
                        names_gamma_c, "mu_t", "mu_c")
      colnames(varcov) <- names_varcov
      rownames(varcov) <- names_varcov
    }
  } else if (estimand %in% c("AO", "ATC")) {
    if (ncol(z) != 0 & method_y %in% c("wls", "logit")) {
      names(gamma) <- paste0("outcome_", names_z)
    }
    names(coef) <- paste0("ps_", names_x_ps)
    if (vcov == TRUE) {
      colnames(varcov) <- c(estimand0, names(coef), names(gamma))
      rownames(varcov) <- c(estimand0, names(coef), names(gamma))
    }
  } else { # estimand == "ATE"
    if (ncol(z) != 0 & method_y %in% c("wls", "logit")) {
      names_gamma_t <- paste0("outcome_t_", names_z)
      names_gamma_c <- paste0("outcome_c_", names_z)
      names(gamma_t) <- names_gamma_t
      names(gamma_c) <- names_gamma_c
    } else {
      names_gamma_t <- NULL
      names_gamma_c <- NULL
    }
    names_coef_t <- paste0("ps_t_", names_x_ps)
    names_coef_c <- paste0("ps_c_", names_x_ps)
    names(coef_t) <- names_coef_t
    names(coef_c) <- names_coef_c
    coef <- list(coef_ps_t = coef_t, coef_ps_c = coef_c)
    ps <- list(ps_t = ps_t, ps_c = ps_c)
    gamma <- list(coef_y_t = gamma_t, coef_y_c = gamma_c)
    y_hat <- list(predicted_y_t = y_hat_t, predicted_y_c = y_hat_c)
    if (vcov == TRUE) {
      names_varcov <- c(estimand0, names_coef_t, names_gamma_t, "mu_t", 
                        names_coef_c, names_gamma_c, "mu_c")
      colnames(varcov) <- names_varcov
      rownames(varcov) <- names_varcov
    }
  }
  ## Some triks for estimand0 == ATT
  if (estimand0 == "ATT") {
    if (method == "eb") {
      ps <- 1 / ps
    } else {
      ps <- 1 - ps
    }
    est <- -est
    coef <- -coef
    response <- 1 - response
    effn <- c(effn[2], effn[1])
  }
  if (estimand != "AO") {
    names(effn) <- c("mu_t", "mu_c")
  }
  list(est = est, 
       coef_ps = coef, 
       coef_y = gamma, 
       varcov = varcov, 
       est_weights = est_weights, 
       ps = ps, 
       predicted_y = y_hat,
       est_sub = est_sub,
       converged = converged,
       effn = effn,
       effn_original = effn_original,
       estimand = estimand0,
       method = method,
       method_y = method_y,
       response = response,
       outcome = outcome,
       original_weights = weights)
}
