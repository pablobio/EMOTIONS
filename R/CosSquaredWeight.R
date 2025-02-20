#' Estimate normalized model`s weights based on the cosine similarity for each model's predictions
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x  data frame containing the daily milking records for the individual
#' @return A vector containing the model normalized weight
CosSquaredWeight<-function(converged_models, x){

  # Step 1: Generate predictions for each model
  predictions <- lapply(converged_models, predict, newdata = x)

  # Convert predictions to a matrix for easier manipulation
  pred_matrix <- do.call(cbind, predictions)

  # Step 2: Compute the average prediction (mean across models)
  avg_prediction <- rowMeans(pred_matrix)

  # Step 3: Compute cosine similarity for each model's predictions
  cosine_similarity <- apply(pred_matrix, 2, function(pred) {
    sum(pred * avg_prediction) / (sqrt(sum(pred^2)) * sqrt(sum(avg_prediction^2)))
  })

  alpha<-2

  # Step 4: Square the cosine similarities to calculate weights
  cosine_squared_weights <- cosine_similarity^alpha

  # Step 5: Normalize the weights so they sum to 1
  normalized_weights <- cosine_squared_weights / sum(cosine_squared_weights)

  return(normalized_weights)

}
