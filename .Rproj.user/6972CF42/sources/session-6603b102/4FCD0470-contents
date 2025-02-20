#' Estimate the Akaike information criterion (AIC), Bayeasian information criterion (BIC), root mean square percentage error (RMSPE) and mean squared error (MAE) for the fitted models
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x  data frame containing the daily milking records for the individual
#' @param trait The name of the column containing daily milking records
#' @importFrom quantreg AIC.rq
#' @return A data frame containing the AIC, BIC, RMSPE and MAE for each fitted model for the individual
GetLacModelsMetrics<-function(converged_models, x, trait){

  # Assuming `converged_models` is your filtered list of models
  metrics_list <- lapply(names(converged_models), function(model_name) {
    model <- converged_models[[model_name]]

    # Calculate predictions
    predictions <- predict(model, newdata = x)

    # Calculate AIC based on model type
    aic_value <- if (inherits(model, "rq")) {
      AIC.rq(model)  # For rq models
    } else if (inherits(model, c("lm", "nls"))) {
      AIC(model)  # For lm and nls models
    } else {
      NA  # Default if model type doesn't match
    }

    # Calculate AIC based on model type
    bic_value <- if (inherits(model, "rq")) {
      AIC.rq(model, k=log(length(predictions)))  # For rq models
    } else if (inherits(model, c("lm", "nls"))) {
      AIC(model, k=log(length(predictions)))  # For lm and nls models
    } else {
      NA  # Default if model type doesn't match
    }

    # Calculate RMSPE
    rmspe_value <- sqrt(mean((x[,trait] - predictions)^2))

    # Calculate MAE
    mae_value <- mean(abs(x[,trait] - predictions))

    # Return metrics as a named list, including model name
    data.frame(Model = model_name, AIC = aic_value, BIC=bic_value, RMSPE = rmspe_value, MAE = mae_value)
  })

  # Combine the list of data frames into a single data frame
  metrics_df <- do.call(rbind, metrics_list)

  return(metrics_df)

}
