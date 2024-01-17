#!/usr/bin Rscript

library(tidyverse)

# Defining the working directory 
setwd("C:/Users/Dell/Desktop/R implementations/10nbl22normal/90replicates/diff_parameters_data/")

# Loading complementary scripts
source("C:/Users/Dell/Desktop/R implementations/functions/DataManipulation.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Metrics.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Plots.R")

# --- Data prep ----
## 1. Loading deconvolution results
res_files <- list.files(pattern = "_res.csv") #rpc, qp
res_list <- lapply(res_files, read.csv)

## 2. Changing the samples names and ordered them
res_list <- lapply(seq_along(res_files), function(x) FormattingSamplesName(res_list[[x]]))

## 3. Defining important variables 
tools <- unique(sub("_.*", "", res_files))
depth_num <- unique(sub("M_.*", "", rownames(res_list[[1]])))
rep_num <- unique(sub(".+_", "", rownames(res_list[[1]])))
cpg_filter <- c("15", "30")

## 4. Creating a single tmp every-info-df
tmp_df <- bind_rows(res_list) %>%
  rownames_to_column() %>%
  add_column(tool = c(rep(tools[1], 180), rep(tools[2], 180))) %>%
  add_column(depth = rep(depth_num, each = 30, 4)) %>%
  add_column(replicate = rep(rep_num, 120)) %>%
  add_column(cpg_filter = rep(cpg_filter, each = 90, 2))

names(tmp_df)[names(tmp_df) == 'rowname'] <- "expected" #changing some col name
tmp_df$expected <- as.numeric(mgsub::mgsub(tmp_df$expected , c("\\d+M_", "_\\d+\\.+\\d+$"), c("", ""))) #changing expected col names
tmp_df$expected <- tmp_df$expected/100 #converting to fractions :))))

## 5. Creating a df containing all replicates
res_df <- tmp_df[order(tmp_df$tool, as.numeric(tmp_df$depth), tmp_df$replicate), ]
res_df <- res_df[, c(4,5,6,7,1,2,3)]
write.csv(res_df, "res_df.csv")
#str(res_df) #for future refs

### 5.1. Creating the relative metrics df
all_dfs <- split(res_df, rep(1:36, each = 10)) #list of all possible dfs

metrics_df <- CreateMetricsDf(all_dfs) %>% #creating the metric df
  add_column(tool = c(rep(tools[1], 18), rep(tools[2], 18))) %>%
  add_column(depth = rep(depth_num, each = 6, 2)) %>%
  add_column(replicate = rep(rep_num, each = 2, 6)) %>%
  add_column(cpg_filter = rep(cpg_filter, 18))

metrics_df <- metrics_df[, c(6,7,8,9,1,2,3,4,5)] #changing the cols order
metrics_df <- metrics_df[order(metrics_df$tool, as.numeric(metrics_df$depth)),] #changing the order of entries
write.csv(metrics_df, "metrics_df.csv")
#str(metrics_df) #for future refs, in case you need to change data type while plotting

## 6. Creating a df which replicates exp have been averaged
### 6.1. Averaging phenotypes expression (three replicates of each samples)
rep_list <- res_df[order(res_df$expected, res_df$cpg_filter), ] %>% #splitting the ordered dataframe by the number of replicates
  split(rep(1:120, each = 3))
  
res_avg_df <- cbind(res_df[order(res_df$expected, res_df$cpg_filter), ] %>% distinct(tool, depth, cpg_filter, expected), #creating the new "averaged" df
                    AveragingPhenos(rep_list)) 
res_avg_df <- res_avg_df[order(res_avg_df$tool, as.numeric(res_avg_df$depth), res_avg_df$cpg_filter), ][-7] #ordering the df
write.csv(res_avg_df, "res_avg_df.csv")
#str(res_avg_df) #for future refs, in case you need to change data type while plotting

### 6.2. Creating the relative metric df
all_dfs <- split(res_df, rep(1:12, each = 10)) #list of all possible dfs

metrics_avg_df <- CreateMetricsDf(all_dfs) %>% #creating the metric df
  add_column(tool = c(rep(tools[1], 6), rep(tools[2], 6))) %>%
  add_column(depth = rep(depth_num, each = 2, 2)) %>%
  add_column(cpg_filter = rep(cpg_filter, 6))

metrics_avg_df <- metrics_avg_df[, c(6,7,8,1,2,3,4,5)] #changing the cols order
metrics_avg_df <- metrics_avg_df[order(metrics_avg_df$tool, as.numeric(metrics_avg_df$depth)),] #changing the order of entries
write.csv(metrics_avg_df, "metrics_avg_df.csv")
#str(metrics_avg_df) #for future refs, in case you need to change data type while plotting

## 7. Creating the benchmark sub-dfs 
setwd("data/")
bench_tools <- c("CIBERSORT", "EpiDISH", "MethAtlas", "PRMeth")
bench_colnames <- c("X", "nbl", "normal")

### CIBERSORT
cbs_list <- lapply(list.files()[1:2], function(x) read.table(x, header = 1)[, 1:3])
cbs_list <- lapply(cbs_list, setNames, bench_colnames)
cbs_list <- lapply(seq_along(cbs_list), function(x) FormattingSamplesName(cbs_list[[x]]))

cbs_df <- bind_rows(cbs_list) %>%
  rownames_to_column("expected") %>%
  add_column(tool = rep(bench_tools[1], 180)) %>%
  add_column(depth = rep(depth_num, each = 30, 2)) %>%
  add_column(replicate = rep(rep_num, 60)) %>%
  add_column(cpg_filter = rep(cpg_filter, each = 90))

cbs_df$expected <- as.numeric(mgsub::mgsub(cbs_df$expected , c("\\d+M_", "_\\d+\\.+\\d+$"), c("", "")))
cbs_df$expected <- cbs_df$expected/100
cbs_df <- cbs_df[, c(4,5,6,7,1,2,3)]

all_dfs <- split(cbs_df, rep(1:60, each = 3))
cbs_avg_df <- cbind(cbs_df %>% dplyr::distinct(tool, depth, cpg_filter, expected), 
                    AveragingPhenos(all_dfs))[-7]

### MethAtlas
matlas_list <- lapply(list.files()[3:4], function(x) as.data.frame(t(read.csv(x, row.names = 1, check.names = FALSE))) %>% rownames_to_column())
matlas_list <- lapply(matlas_list, setNames, bench_colnames)
matlas_list <- lapply(seq_along(matlas_list), function(x) FormattingSamplesName(matlas_list[[x]]))

matlas_df <- bind_rows(matlas_list) %>%
  rownames_to_column("expected") %>%
  add_column(tool = rep(bench_tools[3], 180)) %>%
  add_column(depth = rep(depth_num, each = 30, 2)) %>%
  add_column(replicate = rep(rep_num, 60)) %>%
  add_column(cpg_filter = rep(cpg_filter, each = 90))

matlas_df$expected <- as.numeric(mgsub::mgsub(matlas_df$expected , c("\\d+M_", "_\\d+\\.+\\d+$"), c("", "")))
matlas_df$expected <- matlas_df$expected/100
matlas_df <- matlas_df[, c(4,5,6,7,1,2,3)]

all_dfs <- split(matlas_df, rep(1:60, each = 3))
matlas_avg_df <- cbind(matlas_df %>% dplyr::distinct(tool, depth, cpg_filter, expected), 
                    AveragingPhenos(all_dfs))[-7]

### Complete benchmark df
bench_df <- rbind(cbs_df, matlas_df, res_df)
bench_df <- bench_df[order(bench_df$tool, as.numeric(bench_df$depth), bench_df$cpg_filter, bench_df$replicate), ]

bench_avg_df <- rbind(cbs_avg_df, matlas_avg_df, res_avg_df)
bench_avg_df <- bench_avg_df[order(bench_avg_df$tool, as.numeric(bench_avg_df$depth), bench_avg_df$cpg_filter), ]

setwd("../")
write.csv(bench_df, "bench_df.csv")
write.csv(bench_avg_df, "bench_avg_df.csv")

## 8. Creating the benchmark's metrics df
all_dfs <- split(bench_df, rep(1:72, each = 10)) #listing all possible dfs

bench_metrics_df <- CreateMetricsDf(all_dfs) %>%
  add_column(tool = rep(bench_tools, each = 18)) %>%
  add_column(depth = rep(depth_num, each = 6, 4)) %>%
  add_column(replicate = rep(rep_num, 24)) %>%
  add_column(cpg_filter = rep(cpg_filter, each = 3, 12))

bench_metrics_df <- bench_metrics_df[, c(6,7,8,9,1,2,3,4,5)] #changing the cols order
bench_metrics_df <- bench_metrics_df[order(bench_metrics_df$tool, as.numeric(bench_metrics_df$depth)),] #changing the order of entries
write.csv(bench_metrics_df, "bench_metrics_df.csv")
#str(bench_metrics_df) #for future refs, in case you need to change data type while plotting

all_dfs <- split(bench_avg_df, rep(1:24, each = 10)) #listing all possible dfs

bench_metrics_avg_df <- CreateMetricsDf(all_dfs) %>%
  add_column(tool = rep(bench_tools, each = 6)) %>%
  add_column(depth = rep(depth_num, 8)) %>%
  add_column(cpg_filter = rep(cpg_filter, each = 3, 4))

bench_metrics_avg_df <- bench_metrics_avg_df[, c(6,7,8,1,2,3,4,5)] #changing the cols order
bench_metrics_avg_df <- bench_metrics_avg_df[order(bench_metrics_avg_df$tool, as.numeric(bench_metrics_avg_df$depth)),] #changing the order of entries
write.csv(bench_metrics_avg_df, "bench_metrics_avg_df.csv")
#str(bench_metrics_avg_df) #for future refs, in case you need to change data type while plotting

# --- Loading ----
res_df <- read.csv("res_df.csv", row.names = 1)
res_avg_df <- read.csv("res_avg_df.csv", row.names = 1)

metrics_df <- read.csv("metrics_df.csv", row.names = 1)
metrics_avg_df <- read.csv("metrics_avg_df.csv", row.names = 1)

bench_df <- read.csv("bench_df.csv", row.names = 1)
bench_avg_df <- read.csv("bench_avg_df.csv", row.names = 1)

bench_metrics_df <- read.csv("bench_metrics_df.csv", row.names = 1)
bench_metrics_avg_df <- read.csv("bench_metrics_avg_df.csv", row.names = 1)

# --- Benchmark ----
## Adapt codes according to your aims!

## 1. Barplots
bench_df$depth <- factor(bench_df$depth, levels = unique(bench_df$depth))
bench_avg_df$depth <- factor(bench_avg_df$depth, levels = unique(bench_avg_df$depth))

PlotBars(bench_df) +
  ggh4x::facet_nested(~ cpg_filter ~ tool ~ replicate ~ depth)

PlotBars(bench_avg_df) +
  ggh4x::facet_nested(~ cpg_filter ~ tool ~ depth)

## 2. Scatterplot (nbl tumor fraction)
exp_frac <- unique(bench_df$expected)

AddDotsAes(bench_df %>% 
             ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool, shape = depth)), exp_frac) + 
  labs(color = "Deconvolution tool") +
  ggh4x::facet_nested(~ depth ~ replicate ~ "10nbl22normal" + cpg_filter)

AddDotsAes(bench_avg_df %>% 
             ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool, shape = depth)), exp_frac) + 
  labs(color = "Deconvolution tool") +
  ggh4x::facet_nested(~ depth ~ "10nbl22normal" + cpg_filter)

## 3. Errors histograms
bench_metrics_df$depth <- factor(bench_metrics_df$depth, levels = unique(bench_metrics_df$depth))
bench_metrics_avg_df$depth <- factor(bench_metrics_avg_df$depth, levels = unique(bench_metrics_avg_df$depth))

PlotErrorsAsBars(bench_metrics_df %>% dplyr::select(-PCC, -R2)) +
  ggh4x::facet_nested(~ replicate ~ depth ~ tool + cpg_filter)

PlotErrorsAsBars(bench_metrics_avg_df %>% dplyr::select(-PCC, -R2)) +
  ggh4x::facet_nested(~ depth ~ tool + cpg_filter)

## 4. Metrics dot plot
### With replicates
bench_metrics_df$cpg_filter <- factor(bench_metrics_df$cpg_filter, levels = unique(bench_metrics_df$cpg_filter))

corr_df <- reshape2::melt(bench_metrics_df[, c(1:6)])
corr_df$depth <- factor(corr_df$depth, levels = unique(corr_df$depth))

err_df <- reshape2::melt(bench_metrics_df[, c(1:4, 7:9)])
err_df$depth <- factor(err_df$depth, levels = unique(err_df$depth))
colnames(err_df) <- c("tool", "depth", "replicate", "cpg_filter", "variable2", "value2")

AddMetricsAes(merge(corr_df, err_df) %>% ggplot2::ggplot(aes(x = variable2, y = tool, color = value, size = 1/value2))) + 
  labs(y = "Deconvolution method", x = "Error metric", size = "1/error metric") +
  ggh4x::facet_nested(~ replicate + depth ~ variable + cpg_filter)

### Without replicates
corr_avg_df <- reshape2::melt(bench_metrics_avg_df[, c(1:5)])
err_avg_df <- reshape2::melt(bench_metrics_avg_df[, c(1:3, 6:8)])
colnames(err_avg_df) <- c("tool", "depth", "cpg_filter", "variable2", "value2")

AddMetricsAes(merge(corr_avg_df, err_avg_df) %>% ggplot2::ggplot(aes(x = variable2, y = tool, color = value, size = 1/value2))) + 
  labs(y = "Deconvolution method", x = "Error metric", size = "1/error metric") +
  ggh4x::facet_nested(~ depth ~ variable + cpg_filter)
