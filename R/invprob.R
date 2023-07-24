## inverse probability weights
invprob <- function (ps, response, estimand, weights) {
  invp <- numeric(length(response))
  if (estimand == "ATC") {
    invp[response == 1] <- (1 - ps[response == 1]) / ps[response == 1] * 
                             weights[response == 1] / sum(weights[response == 0])
    invp[response == 0] <- weights[response == 0] / sum(weights[response == 0])
  } else { # estimand == "AO")
    invp[response == 1] <- 1 / ps[response == 1] * weights[response == 1] / sum(weights)
    invp[response == 0] <- 0
  }
  invp
}
