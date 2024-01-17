#!/usr/bin Rscript

library(tidyverse)

# Defining the working directory 
setwd("C:/Users/Dell/Desktop/R implementations/10nbl22normal/9vitros/diff_parameters_data/")

# Loading complementary scripts
source("C:/Users/Dell/Desktop/R implementations/functions/DataManipulation.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Plots.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Metrics.R")

# --- Data prep ----
## 1. Loading deconvolution results
res_files <- list.files(pattern = "_res.csv")
res_list <- lapply(res_files, read.csv) #rpc, qp
res_list <- lapply(seq_along(res_list), function(x) res_list[[x]][gtools::mixedorder(res_list[[x]]$X),]) #re-ordering entries based on X

## 2. Defining important features
tools <- unique(sub("_.*", "", res_files))
feat_sel <- c("100hyper", "100top")
cpg_filter <- c("15", "30")
exp_frac <-  c(0.5, 0.25, 0.1045, 0.0399, 0.0146, 0.0053, 0.0019, 0, 1)

## 3. Creating a single every-info-df
res_df <- bind_rows(res_list) %>%
  add_column(tool = c(rep(tools[1], each = 9, 4), rep(tools[2], each = 9, 4))) %>%
  add_column(feature_select = rep(feat_sel, each = 18, 2)) %>%
  add_column(cpg_filter = rep(cpg_filter, each = 9, 4))

res_df <- res_df[, c(4,5,6,1,2,3)] #re-ordering cols
names(res_df)[names(res_df) == 'X'] <- "expected"
res_df$expected <- rep(exp_frac, 8)
res_df <- res_df[order(res_df$tool, res_df$feature_select, res_df$cpg_filter, res_df$expected), ] #sorting values based on own preference
write.csv(res_df, "res_df.csv")
#str(res_df) #for future refs

## 4. Creating an errors df
all_dfs <- split(res_df, rep(1:8, each = 9)) #list of all possible dfs

errors_df <- CreateMetricsDf(all_dfs) %>% #creating the metric df
  add_column(tool = c(rep(tools[1], 4), rep(tools[2], 4))) %>%
  add_column(feature_select = rep(feat_sel,  each = 2, 2)) %>%
  add_column(cpg_filter = rep(cpg_filter, 4))

errors_df <- errors_df[, c(6,7,8,1,2,3,4,5)] #changing the cols order
write.csv(errors_df, "errors_df.csv")
#str(errors_df) #for future refs, in case you need to change data type while plotting

## 5. Creating the benchmark df
setwd("data/")
bench_tools <- c("CIBERSORT", "EpiDISH", "MethAtlas", "PRMeth")
bench_colnames <- c("expected", "nbl", "normal")

### CIBERSORT
cbs_list <- lapply(list.files()[1:4], function(x) read.table(x, header = 1)[, 1:3])
cbs_list <- lapply(cbs_list, setNames, bench_colnames)

### MethAtlas
matlas_list <- lapply(list.files()[5:8], function(x) as.data.frame(t(read.csv(x, row.names = 1, check.names = FALSE))) %>% rownames_to_column())
matlas_list <- lapply(matlas_list, setNames, bench_colnames)

### Complete benchmark df
bench_df <- bind_rows(cbs_list, matlas_list) %>%
  add_column(tool = c(rep(bench_tools[1], each = 36), rep(bench_tools[3], each = 36))) %>%
  add_column(feature_select = rep(feat_sel, each = 9, 4)) %>%
  add_column(cpg_filter = rep(cpg_filter, each = 18, 2))

bench_df$expected <- rep(exp_frac, 8)
bench_df <- bench_df[, c(4,5,6,1,2,3)]
bench_df <- rbind(bench_df, res_df) 
bench_df <- bench_df[order(bench_df$tool, bench_df$feature_select, bench_df$cpg_filter, bench_df$expected), ] #re-ordering values

setwd("../")
write.csv(bench_df, "bench_df.csv")

## 6. Creating the benchmark error df
all_dfs <- split(bench_df, rep(1:16, each = 9)) #list of all possible dfs

bench_errors_df <- CreateMetricsDf(all_dfs) %>% #creating the error df
  add_column(tool = rep(bench_tools, each = 4)) %>%
  add_column(feature_select = rep(feat_sel, each = 2, 4)) %>%
  add_column(cpg_filter = rep(cpg_filter, 8)) 

bench_errors_df <- bench_errors_df[, c(6,7,8,1,2,3,4,5)] #changing the cols order
write.csv(bench_errors_df, "bench_errors_df.csv")
#str(bench_errors_df) #for future refs, in case you need to change data type while plotting

# --- Loading ----
res_df <- read.csv("res_df.csv", row.names = 1)
bench_df <- read.csv("bench_df.csv", row.names = 1)

metrics_df <- read.csv("metrics_df.csv", row.names = 1)
bench_metrics_df <- read.csv("bench_metrics_df.csv", row.names = 1)

# --- Benchmark ----
## Adapt codes according to your aims!

## 1. Barplots
PlotBars(bench_df %>% dplyr::filter(feature_select == "100top")) +
  ggh4x::facet_nested(~ cpg_filter ~ tool)

## 2. Scatterplot (nbl tumor fraction)
AddDotsAes(bench_df %>% dplyr::filter(feature_select == "100top") %>% 
             ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool)), exp_frac) + 
  labs(color = "Deconvolution tool") +
  ggh4x::facet_nested(~ "10nbl22normal" + cpg_filter)

## 3. Errors histograms
PlotErrorsAsBars(bench_errors_df %>% dplyr::select(-PCC, -R2)) +
  ggh4x::facet_nested(~ tool ~ cpg_filter)

## 4. Metrics dot plot
corr_df <- reshape2::melt(bench_errors_df[, c(1:5)] %>% dplyr::filter(feature_select == "100top") %>% dplyr::select(-feature_select))
err_df <- reshape2::melt(bench_errors_df[, c(1:3, 6:8)] %>% dplyr::filter(feature_select == "100top") %>% dplyr::select(-feature_select))
colnames(err_df) <- c("tool", "cpg_filter", "variable2", "value2")

AddMetricsAes(merge(corr_df, err_df) %>% ggplot2::ggplot(aes(x = variable2, y = tool, color = value, size = 1/value2))) + 
  labs(y = "Deconvolution method", x = "Error metric", size = "1/error metric") +
  ggh4x::facet_nested(~ cpg_filter ~ variable)
