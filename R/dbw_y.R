dbw_y <- function (formula_y, df, fulldf, 
                   method_y, est_weights) {
  df <- data.frame(df, est_weights = est_weights)
  family <- stats::gaussian
  if (method_y %in% c("logit", "gambinom")) {
    family <- stats::quasibinomial
  }
  if (method_y %in% c("wls", "logit")) {
    result <- stats::glm(formula = formula_y, family = family, data = df, 
                         weights = est_weights)
    pred <- stats::predict(object = result, newdata = fulldf, type = "response")
    return(list(pred = pred, gamma = result$coefficients))
  } else { # method_y %in% c("gam", "gambinom"))
    result <- mgcv::gam(formula = formula_y, family = family, data = df, 
                        weights = est_weights)
    pred <- c(stats::predict(object = result, newdata = fulldf, type = "response"))
  }
  list(pred = pred, gamma = NULL)
}
