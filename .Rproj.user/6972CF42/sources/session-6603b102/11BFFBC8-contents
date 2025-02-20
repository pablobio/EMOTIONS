#' Plot the actual daily milk daily production and the predicted values obtained by the ensemble model
#' @param data A data frame containing the daily milking records
#' @param ID The name of the column containing the unique IDs of the individuals
#' @param trait The name of the column containing daily milking records
#' @param metric The name of the strategy used obtained the predicted values through the ensemble model
#' @param dim The name of the column containing days in milk records
#' @param col The colors of the actual and predicted values
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment labs theme theme_minimal element_text
#' @return A plot with the actual and predicted daily milk production across the days in milk
#' @export
PlotWeightLac<-function(data,ID, trait,metric,dim,col=c("red","blue")){

  utils::globalVariables(c(".data"))

  data<-data$production[[ID]]

  data %>% ggplot2::ggplot() +
    # Plot the observed data points
    geom_point(aes(x = .data[[dim]], y = .data[[trait]]), color = col[1] , size = 2) +
    # Plot the ensemble regression line
    geom_line(data = data, aes(x = .data[[dim]], y = .data[[metric]]), color = col[2], size = 1) +
    labs(title = "", x = "DIM", y = "Production") +
    theme_minimal() +
    theme(axis.text.x = element_text( size=15), axis.text.y = element_text(size=15), axis.title = element_text(size=15))

}
