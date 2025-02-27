#' Create a line plot that shows the range of the ranks obtained for each model across the individuals
#' @param LacCurveFit The object obtained from the LacCurveFit function
#' @param metric The name of the metric to be use to plot the model´s ranks
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment labs theme theme_minimal element_text geom_text
#' @importFrom dplyr n_distinct
#' @return A line plot that shows the range of the ranks obtained for each model across the individuals
#' @export
ModelRankRange<-function(LacCurveFit,metric="AIC_rank"){

  data <- bind_rows(LacCurveFit$models_weight, .id = "ID")

  range_df <- data %>%
    group_by(Model) %>%
    summarise(min_value = min(!!sym(metric)), max_value = max(!!sym(metric)),
              n_individuals = n_distinct(ID))

  range_df$Model<-factor(range_df$Model,levels = range_df$Model[order(range_df$Model,decreasing = T)])


  ggplot(range_df, aes(x = min_value , y = Model, xend = max_value)) +
    geom_segment(aes(yend = Model, xend = max_value), size = 1) +
    geom_text(aes(label = paste0("n=", n_individuals), x = max_value + 0.02),size = 5, hjust = 0)+
    labs(x = "Rank", y = "Model", title = "") +
    theme(axis.text.x = element_text( size=15), axis.text.y = element_text(size=15), axis.title = element_text(size=15))


}
