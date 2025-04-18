% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotWeightLac.R
\name{PlotWeightLac}
\alias{PlotWeightLac}
\title{Plot the actual and predicted daily milk production obtained by the ensemble model}
\usage{
PlotWeightLac(
  data,
  ID,
  trait,
  metric,
  dim,
  col = c("red", "blue"),
  point_size = 2,
  line_size = 1,
  axis_text_size = 15,
  axis_title_size = 15
)
}
\arguments{
\item{data}{The object generated by the LacCurveFit function}

\item{ID}{The ID of the individual whose daily milking records will be plotted}

\item{trait}{The name of the column containing the daily milking records}

\item{metric}{The name of the strategy used to obtain the predicted values through the ensemble model}

\item{dim}{The name of the column containing the days in milk (DIM) records}

\item{col}{A vector defining the colors for the actual and predicted values}

\item{point_size}{Numeric value indicating the size of the observed data points in the plot}

\item{line_size}{Numeric value indicating the thickness of the regression line}

\item{axis_text_size}{Numeric value defining the font size of the axis tick labels}

\item{axis_title_size}{Numeric value defining the font size of the axis titles}
}
\value{
A plot showing the actual and predicted daily milk production across the days in milk
}
\description{
Plot the actual and predicted daily milk production obtained by the ensemble model
}
