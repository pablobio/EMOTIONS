# EnsembleLacs: An R package for the implementation of average modelling algorithms in lactation curves

**INDEX**

-   📝 [Introduction](#introduction)
-   🛠️ [Installation](#installation)
-   💻 [Analysis](#analysis)
    -   [Data Input](#data-input)
-   📚 [References](#references)
-   🤝 [Authors](#authors)
-   📜 [Licence](#licence)

## **Introduction**

The mathematical representation of lactation curves has important applicability in different areas of animal science. Models that simulate milk production under various conditions are valuable for physiologists, nutritionists, and geneticists to study mammary gland function and test hypotheses, providing support to management decisions related to timing and efficiency. Several models have been proposed to model lactation curves for dairy species in the last decades. These models mainly differ from each other in the regression type (linear or nonlinear), number of parameter, relationship among the parameters, and the representation of lactation patterns, such as peak yield, time at peak, and persistency. Despite the advantages to have a wide range of model options to choose and the specificity of each individual model, the selection of a single model to represent the lactation curve might raise some limitations. Usually, the selection of the best model among the available options is based in the comparison among models for a metric, such as the Akaike information criteria (AIC) and Bayesian information criteria (BIC). The metric-based approach heavily relies on the assumption that the selected metric is a good criterion for model selection. However, these metrics can sometimes fail to correctly identify the best model, especially when there is no clear "best" model. Additionally, the use of a single model, based on a metric selection, can result in overfitting and bias introduction of individual models (especially if the number of variables is large or the dataset is noisy).

Ensemble models and model averaging are powerful techniques in predictive modeling that enhance robustness, accuracy, and generalization. Instead of relying on a single model, ensemble methods combine multiple models to reduce variance, mitigate overfitting, and improve predictive performance. Model averaging, incorporates model uncertainty by weighting predictions according to their posterior probabilities, leading to more reliable estimates. These approaches are especially valuable when individual models exhibit variability in performance across different datasets or conditions. By integrating diverse model perspectives, ensemble methods improve stability and resilience, making them particularly useful in complex systems such as biological and ecological modeling, where uncertainty quantification is crucial. The use of ensemble models in the context of modelling lactation curves can be useful based on the wide variety of models available and the specificity of lactation curves for some species and/or conditions.

The EnsembleLacs package provides a series of tools for the fit of 48 different lactation curve models previously reported in the literature. Once the data if fitted for each model, averaged predictions are obtaining using bagging based on AIC, BIC, root mean square percentage error (RMSPE), mean squared error (MAE) and variance. Additionally, the package also reports the milk daily records predictions based on Bayesian Model Averaging (BMA) and the cosine similarity for each model's predictions. The visualization of the model ranking across all the individual predictions is possible through functions (RidgeModels and ModelRankRange), allowing to better understand the weights assigned to each model. In addition, the predicted and actual daily milk records can be compared in the EnsembleLacs package using the PlotWeightLac function. Finally, the package also allows the user to estimate resilience indicators based on lag1-autocorrelation, logarithm of residual variance and residual skeweness using the predicted daily milking records.
