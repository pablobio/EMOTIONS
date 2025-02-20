#' Estimate normalized model`s weights based on a  Expectation–Maximization (EM) algorithm using a gamma distribution
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x  data frame containing the daily milking records for the individual
#' @param trait The name of the column containing daily milking records
#' @return A vector containing the model normalized weight based on the Bayeasian Model Average (BMA) appaoch proposed by Duan et al. (2006) with an adaptation for to use a gamma distribution
BMAweight_gamma <- function(converged_models, x, trait) {

  # Step 1: Generate predictions for each model
  predictions <- lapply(converged_models, predict, newdata = x)
  predictions <- as.data.frame(predictions)

  k <- ncol(predictions)  # Number of models
  t <- nrow(predictions)  # Number of data points

  observed <- x[, trait]

  # Initial weights and variance estimates
  weights <- rep(1 / k, k)
  variances <- colMeans((observed - predictions)^2)  # Initial variance estimate

  max_iter <- 100
  tol <- 1e-6
  log_likelihood <- -Inf

  for (iter in 1:max_iter) {
    # E-step
    z <- matrix(0, nrow = t, ncol = k)
    for (c in 1:k) {
      pred_mean <- predictions[, c]
      pred_var <- variances[c]

      # Convert variance to shape and scale parameters for the Gamma distribution
      shape <- (pred_mean^2) / pred_var
      scale <- pred_var / pred_mean

      # Ensure valid shape and scale values (avoid division by zero or negative values)
      shape[shape <= 0] <- 1e-6
      scale[scale <= 0] <- 1e-6

      # Compute Gamma density
      g.dense <- dgamma(observed, shape = shape, scale = scale)

      # Compute responsibility z
      z[, c] <- weights[c] * g.dense
    }

    z <- z / rowSums(z)  # Normalize z

    z[is.na(z)]<-1e-6

    # M-step: Update weights
    weights <- colMeans(z)

    alpha<-4

    weights<-weights^alpha

    weights <- weights / sum(weights)  # Normalize weights to sum to 1

    # Update variance estimate
    variances <- colSums(z * (observed - predictions)^2) / colSums(z)

    # Log-likelihood
    new_log_likelihood <- sum(log(rowSums(z)))

    # Convergence Check
    if (abs(new_log_likelihood - log_likelihood) < tol) {
      break
    }
    log_likelihood <- new_log_likelihood
  }

  return(weights)
}
