# Plots v.0.1
# Description: an R script with a series of functions to plot your deconvolution results.
#
# Functions:
#   AddDotsAes(df, vec): adds the specified aesthetic to the input ggplot df. 
#   AddMetricsAes(df): adds the specified aesthetic to the input ggplot df. 
#   PlotBars(df): plots the df deconvolution results as bar plots. 
#   PlotErrorsAsBars(df): plots the specified error metrics of the input df as grouped histograms.
#   PlotHeatmap(df): plots the deconvolution result as heatmap, following a normalization procedure in log10 scale. 
#   PlotPCC(df, vec): plots in a scatter plot the PCC between expected and observed fractions of the input df. 
# ------------------------------------------------------------------------------------------------------------------------------------------------------

#!/usr/bin/env R script

AddDotsAes <- function(df, exp_frac) {
  ### df: a ggplot df which attributes have been already selected
  ### exp_frac: a numerical vector containing the expected tumor fraction 
  ### Returns the corresponding dot plot after transforming values in a log10 scale
  
  return(df +
           ggplot2::geom_point() + 
           scale_color_brewer(palette = "Paired") +
           scale_shape_manual(values = c(15, 16, 18)) +
           scale_x_continuous(trans = "log10", breaks = exp_frac + 0.000001, labels = exp_frac) +
           scale_y_continuous(trans = "log10", breaks = exp_frac + 0.000001, labels = exp_frac) +
           geom_abline(aes(intercept = 0 + 0.000001, slope = 1)) +
           labs(x = "Expected fraction (log 10)", y = "Observed fraction (log 10)", shape = "Depth") +
           theme_bw() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.85, hjust = 0.85), 
                 axis.ticks = element_blank(), 
                 legend.position = "bottom", 
                 strip.text = element_text(size = 10),  
                 strip.background = element_rect(colour = "black", fill = "lightgray", linetype = "solid")) +
           guides(colour = guide_legend(override.aes = list(size=3)), shape = guide_legend(override.aes = list(size=3)))
  )
}

AddMetricsAes <- function(df) {
  ### df: a ggplot df which attributes have been already selected
  ### Returns a dot plot of each metric with the following specified aesthetics
  
  return(df +
           ggplot2::geom_point() +
           scale_color_viridis_c(direction = -1) + 
           labs(x = "Depth", color = "PCC/R2", size = "1/error_metric") + 
           theme_bw() +
           theme(axis.ticks = element_blank(), 
                 panel.grid.major = element_blank(), 
                 strip.text.x = element_text(size = 10), 
                 strip.text.y = element_text(size = 10),
                 strip.background = element_rect(colour = "black", fill = "lightgray", linetype = "solid")) +
           ggh4x::facet_nested(~ "Error metrics" + variable2 ~ "Coefficients" + variable) 
  )
}

PlotBars <- function(df) {
  ### df: a df containing the deconvolution results
  ### Returns a bar plot displaying, for each sample, the inferred nbl and normal fractions
  
  df$expected <- factor(df$expected, levels = unique(df$expected))
  
  return(df %>%
           reshape2::melt() %>%
           ggplot2::ggplot(aes(x = expected, y = value, fill = variable, label = format(round(value, 4)))) +
            geom_bar(stat = "identity") +
            geom_text(size = 3, position = position_stack(vjust = 0.3)) +
            scale_fill_brewer(palette = "Paired") +        
            labs(x = "Expected fraction", y = "Observed fraction", fill = "Entity") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.85, hjust = 0.85), 
                  axis.ticks = element_blank())
         )
}

PlotErrorsAsBars <- function(df) {
  ### df: a df containing the error values AAD, MedAE, RMSE of each given sample
  ### Returns a grouped bar plot displaying the specified error metrics 
  
  return(df %>%
           reshape::melt() %>%
           ggplot2::ggplot(aes(x = variable, y = value, fill = variable)) +
           geom_bar(position = "dodge", stat = "identity") +
           labs(x = "Error metric", y = "Value") +
           scale_fill_brewer(palette = "Paired") + 
           theme_bw() +
           theme(axis.text.x = element_text(angle = 0, vjust = 0.85, hjust = 0.85), 
                 legend.position = "none")
         )
}

PlotHeatmap <- function(df) {
  ### df: a df containing the deconvolution results
  ### Returns the corresponding heatmap in a log10 scale
  
  df$expected <- factor(df$expected, levels = unique(df$expected))
  
  return(df %>%
           reshape2::melt() %>%
           ggplot2::ggplot(aes(x = expected, y = variable)) +
            geom_raster(aes(fill = log10(value + 0.000001))) +
            scale_fill_gradient(low = "gray80", high = "purple4")  +
            labs(x = "Expected tumor fractions", y = "Entity", fill = "Tumor fraction (log10)") +
            theme_bw() + 
            theme(axis.line = element_blank(), 
                  axis.text.x = element_text(angle = 90, size = 7, vjust = 0.3),
                  axis.ticks = element_blank()) +
           facet_wrap(~ modality)
         )
}

PlotPCC <- function(df, exp_frac) {
  ### df: a df containing the deconvolution results
  ### exp_frac: a numerical vector containing the expected tumor fractions of the given samples
  ### Returns the PCC scatter plot in log scale
  
  exp_frac <- sort(exp_frac)

  return(df %>%
           ggpubr::ggscatter(x = "expected", y = "nbl", shape = 19, size = 2, 
                             add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), 
                             conf.int = TRUE, 
                             xlab = "Expected fraction (log10)", ylab = "Observed fractions (log10)") +
           geom_abline(aes(intercept = 0 + 0.000001, slope = 1)) +
           scale_x_continuous(trans = "log10", breaks = exp_frac + 0.000001, labels = round(exp_frac, 5)) +
           scale_y_continuous(trans = "log10") +
           theme_bw() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.85, hjust = 0.85)) +
           ggpubr::stat_cor(aes(label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")), r.accuracy = 0.001, p.accuracy = 0.001)
  )
}
