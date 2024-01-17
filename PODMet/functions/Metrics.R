# Benchmark metrics v.0.1
# Description: an R script with a series of functions to compute a series of metrics able to estimate deconvolution algorithms
# accuracy in estimating tumor fractions.
#
# Functions:
#   ComputeMetrics(df, vector): a function that compute a series of metrics. 
# ------------------------------------------------------------------------------------------------------------------------------------------------------

#!/usr/bin/env R script

ComputeMetrics <- function(df, metric = c("aad", "medae", "pcc", "r2", "rmse")) {
  ### df: a df containing the deconvolution results you want to compute the metrics on 
  ### metric: a list of metrics
  ###   *Average Absolute Deviation (AAD): the mean of the absolute difference between observed and expected
  ###   *Median Absolute Error (MedAE): the median of the absolute difference between observed and expected
  ###   *Pearson Correlation coefficient (PCC): the pearson correlation between observed and expected
  ###   *Coefficient of determination (R2): the square value of the correlation between observed and expected
  ###   *Root Mean Square Error (RMSE): the square root of the mean of all errors square
  ### Return the selected error metric value of the input df
  
  if (metric == "aad") {
    return((1/length(df$nbl))*sum(abs((df$nbl) - (df$expected))))
  }
  else if (metric == "medae") {
    return(stats::median(abs((df$nbl) - (df$expected))))
  }
  else if (metric == "pcc") {
    return(stats::cor(df$expected, df$nbl, method = "pearson"))
  }
  else if (metric == "r2") {
      return(stats::cor(df$expected, df$nbl) ^ 2)
  }
  else if (metric == "rmse") {
    return(Metrics::rmse(df$expected, df$nbl))
  }
  else{stop("Invalid metric.")}
}
