#!/usr/bin Rscript

library(quadprog)
library(tidyverse)

# Defining the working directory 
setwd("C:/Users/Dell/Desktop/R implementations/10nbl22normal/90replicates/100top_30/")

# Loading complementary scripts
source("C:/Users/Dell/Desktop/R implementations/functions/DataManipulation.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Metrics.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Plots.R")
source("C:/Users/Dell/Desktop/R implementations/functions/prmeth.R")

# --- Data prep ----
## 1. Loading deconvolution results
res_files <- list.files(pattern = "_res.csv")
res_list <- lapply(res_files, read.csv)

## 2. Changing the samples names and ordered them
res_list <- lapply(seq_along(res_files), function(x) FormattingSamplesName(res_list[[x]]))

## 3. Defining important variables 
tools <- unique(sub("_.*", "", res_files))
mods <- mgsub::mgsub(res_files, c("^[A-Za-z]+_", "_res.csv$", "_"), c("", "", "."))
depth_num <- unique(sub("M_.*", "", rownames(res_list[[1]])))
rep_num <- unique(sub(".+_", "", rownames(res_list[[1]])))

## 4. Creating a single tmp every-info-df
tmp_df <- bind_rows(res_list) %>%
  rownames_to_column() %>%
  add_column(tool = c(rep(tools[1], 360), rep(tools[2], 540))) %>%
  add_column(depth = rep(depth_num, each = 30, 10)) %>%
  add_column(replicate = rep(rep_num, 300)) %>%
  add_column(modality = rep(mods, each = 90))

names(tmp_df)[names(tmp_df) == 'rowname'] <- "expected" #changing some col name
names(tmp_df)[names(tmp_df) == 'X.1'] <- "unknown"
tmp_df$expected <- as.numeric(mgsub::mgsub(tmp_df$expected , c("\\d+M_", "_\\d+\\.+\\d+$"), c("", ""))) #changing expected col names
tmp_df$expected <- tmp_df$expected/100 #converting to fractions :))))

## 5. Creating a df containing all replicates
res_df <- tmp_df[order(tmp_df$tool, as.numeric(tmp_df$depth), tmp_df$modality, tmp_df$replicate), ]
res_df <- res_df[, c(5,8,6,7,1,2,3,4)]
write.csv(res_df, "res_df.csv")
#str(res_df) #for future refs, in case you need to change data type while plotting

### 5.1. Creating the relative metrics df
all_dfs <- split(res_df, rep(1:90, each = 10)) #list of all possible dfs

metrics_df <- CreateMetricsDf(all_dfs) %>% #creating the metric df
  add_column(tool = c(rep(tools[1], 36), rep(tools[2], 54))) %>%
  add_column(modality = c(rep(mods[1:4], each = 3, 3), rep(mods[5:10], each = 3, 3))) %>%
  add_column(depth = c(rep(depth_num, each = 12), rep(depth_num, each = 18))) %>%
  add_column(replicate = rep(rep_num, 30))

metrics_df <- metrics_df[, c(6,7,8,9,1,2,3,4,5)] #changing the cols order
metrics_df <- metrics_df[order(metrics_df$tool, metrics_df$modality, as.numeric(metrics_df$depth)),] #changing the order of entries
write.csv(metrics_df, "metrics_df.csv")
#str(metrics_df) #for future refs, in case you need to change data type while plotting

## 6. Creating a df which replicates exp have been averaged
### 6.1. Averaging phenotypes expression (three replicates of each samples)
rep_list <- res_df[order(res_df$expected), ] %>% #splitting the ordered dataframe by the number of replicates
  split(rep(1:300, each = 3))

res_avg_df <- cbind(res_df[order(res_df$expected), ] %>% dplyr::distinct(tool, modality, depth, expected), #creating the new "averaged" df
                AveragingPhenos(rep_list)) 
res_avg_df <- res_avg_df[order(res_avg_df$tool, as.numeric(res_avg_df$depth), res_avg_df$modality), ] #ordering the df
write.csv(res_avg_df, "res_avg_df.csv")
#str(res_avg_df) #for future refs

### 6.2. Creating the relative metric df
all_dfs <- split(res_avg_df, rep(1:30, each = 10)) #list of all possible dfs

metrics_avg_df <- CreateMetricsDf(all_dfs) %>% #creating the metric df
  add_column(tool = c(rep(tools[1], 12), rep(tools[2], 18))) %>%
  add_column(modality = c(rep(mods[1:4], 3), rep(mods[5:10], 3))) %>%
  add_column(depth = rep(depth_num, 10)) 

metrics_avg_df <- metrics_avg_df[, c(6,7,8,1,2,3,4,5)] #changing the cols order
metrics_avg_df <- metrics_avg_df[order(metrics_avg_df$tool, metrics_avg_df$modality, as.numeric(metrics_avg_df$depth)),] #changing the order of entries
write.csv(metrics_avg_df, "metrics_avg_df.csv")
#str(metrics_df) #for future refs, in case you need to change data type while plotting

## 7. Creating the benchmark sub-dfs 
### CIBERSORT
cbs_res <- read.table("data/results_2perm_test_CIBERSORT_10NGSnbl22NGSnormal_filter30_100top_GRCh37_iter1.txt", header = 1)[, 1:3]
names(cbs_res)[names(cbs_res) == 'Mixture'] <- "X" 
cbs_df <- CreateBenchDf(cbs_res, tool = "CIBERSORT", depth_num, rep_num)
cbs_avg_df <- CreateAvgBenchDf(cbs_res, tool = "CIBERSORT", depth_num, rep_num)

### EpiDISH (RPC)
epi_df <- res_df %>% dplyr::filter(tool == "EpiDISH", modality == "rpc") %>% dplyr::select(-modality, -unknown)
epi_avg_df <- res_avg_df %>% dplyr::filter(tool == "EpiDISH", modality == "rpc") %>% dplyr::select(-modality, -unknown)

### MethAtlas
matlas_res <- as.data.frame(t(read.csv("data/test_methatlas_10NGSnbl22NGSnormal_filter30_100top_GRCh37_iter1_deconv_output.csv", row.names = 1, check.names = FALSE))) %>%
  rownames_to_column("X") 
matlas_df <- CreateBenchDf(matlas_res, tool = "MethAtlas", depth_num, rep_num)
matlas_avg_df <- CreateAvgBenchDf(matlas_res, tool = "MethAtlas", depth_num, rep_num)

### PRMeth (QP)
pr_df <- res_df %>% dplyr::filter(tool == "PRMeth", modality == "qp") %>% dplyr::select(-modality, -unknown)
pr_avg_df <- res_avg_df %>% dplyr::filter(tool == "PRMeth", modality == "qp") %>% dplyr::select(-modality, -unknown)

### 7.1. Creating the dfs
bench_df <- rbind(cbs_df, epi_df, matlas_df, pr_df)
write.csv(bench_df, "bench_df.csv")

bench_avg_df <- rbind(cbs_avg_df, epi_avg_df, matlas_avg_df, pr_avg_df)
write.csv(bench_avg_df, "bench_avg_df.csv")

## 8. Creating the benchmarks metrics dfs
bench_tools <- c("CIBERSORT", "EpiDISH", "MethAtlas", "PRMeth")

### 8.1. Not averaged df
all_dfs <- split(bench_df, rep(1:36, each = 10)) #listing all possible dfs

bench_metrics_df <- CreateMetricsDf(all_dfs) %>%
  add_column(tool = sort(rep(bench_tools, each = 9))) %>%
  add_column(depth = rep(depth_num, each = 1, 12)) %>%
  add_column(replicate = rep(rep_num, each = 3, 4))

bench_metrics_df <- bench_metrics_df[, c(6,7,8,1,2,3,4,5)] #changing the cols order
bench_metrics_df <- bench_metrics_df[order(bench_metrics_df$tool, as.numeric(bench_metrics_df$depth)),] #changing the order of entries
write.csv(bench_metrics_df, "bench_metrics_df.csv")
#str(bench_metrics_df) #for future refs, in case you need to change data type while plotting

### 8.2. Averaged df
all_dfs <- split(bench_avg_df, rep(1:12, each = 10)) #listing all possible dfs

bench_metrics_avg_df <- CreateMetricsDf(all_dfs)%>%
  add_column(tool = sort(rep(bench_tools, each = 3))) %>%
  add_column(depth = rep(depth_num, 4))

bench_metrics_avg_df <- bench_metrics_avg_df[, c(6,7,1,2,3,4,5)] #changing the cols order
bench_metrics_avg_df <- bench_metrics_avg_df[order(bench_metrics_avg_df$tool, as.numeric(bench_metrics_avg_df$depth)),] #changing the order of entries
write.csv(bench_metrics_avg_df, "bench_metrics_avg_df.csv")
#str(bench_metrics_df) #for future refs, in case you need to change data type while plotting

# --- Loading ----
res_df <- read.csv("res_df.csv", row.names = 1)
res_avg_df <- read.csv("res_avg_df.csv", row.names = 1)

metrics_df <- read.csv("metrics_df.csv", row.names = 1)
metrics_avg_df <- read.csv("metrics_avg_df.csv", row.names = 1)

bench_df <- read.csv("bench_df.csv", row.names = 1)
bench_avg_df <- read.csv("bench_avg_df.csv", row.names = 1)

bench_metrics_df <- read.csv("bench_metrics_df.csv", row.names = 1)
bench_metrics_avg_df <- read.csv("bench_metrics_avg_df.csv", row.names = 1)

# --- Intra-comparison between methods ----
## Adapt each plot code according to your aims!

## 1. Correlation vs error metric dot plot 
### 1.1. Creating errors and coefficients dfs
corr_df <- reshape2::melt(metrics_df[, c(1:6)])
corr_df$depth <- factor(corr_df$depth, levels = unique(corr_df$depth))
corr_df$replicate <- factor(corr_df$replicate, levels = unique(corr_df$replicate))

corr_avg_df <- reshape2::melt(metrics_avg_df[, c(1:5)])
corr_avg_df$depth <- factor(corr_avg_df$depth, levels = unique(corr_avg_df$depth))

err_df <- reshape2::melt(metrics_df[, c(1:4, 7:9)])
err_df$depth <- factor(err_df$depth, levels = unique(err_df$depth))
err_df$replicate <- factor(err_df$replicate, levels = unique(err_df$replicate))
colnames(err_df) <- c("tool", "modality", "depth", "replicate", "variable2", "value2")

err_avg_df <- reshape2::melt(metrics_avg_df[, c(1:3, 6:8)])
err_avg_df$depth <- factor(err_avg_df$depth, levels = unique(err_avg_df$depth))
colnames(err_avg_df) <- c("tool", "modality", "depth", "variable2", "value2")

### 1.2. Plotting
AddMetricsAes(merge(corr_df, err_df) %>% dplyr::filter(tool == "EpiDISH") %>%
                ggplot2::ggplot(aes(x = depth, y = modality, color = value, size = 1/value2))) + 
  labs(y = "Modality", size = "1/error metric")  +
  ggh4x::facet_nested(~ variable2 ~ replicate ~ "10nbl22normal"+ variable)

AddMetricsAes(merge(corr_avg_df, err_avg_df) %>% dplyr::filter(tool == "EpiDISH") %>%
                ggplot2::ggplot(aes(x = depth, y = modality, color = value, size = 1/value2))) + 
  labs(y = "Modality", size = "1/error metric")  +
  ggh4x::facet_nested(~ variable2 ~"10nbl22normal"+ variable)

## 2. Barplots
### 2.1. With replicates
res_df$expected <- factor(res_df$expected, levels = unique(res_df$expected)) 
res_df$depth <- factor(res_df$depth, levels = gtools::mixedsort(unique(res_df$depth))) 
res_df$replicate <- factor(res_df$replicate, levels = gtools::mixedsort(unique(res_df$replicate))) 

PlotBars(res_df[res_df$tool == "EpiDISH", -1][, 1:6]) + 
  theme(strip.background = element_rect(colour = "black", fill = "lightgray", linetype = "solid")) +
  ggh4x::facet_nested(~ modality ~ depth ~ replicate)

### 2.2. Without replicates
res_avg_df$expected <- factor(res_avg_df$expected, levels = unique(res_avg_df$expected)) 
res_avg_df$depth <- factor(res_avg_df$depth, levels = gtools::mixedsort(unique(res_avg_df$depth))) 

PlotBars(res_avg_df[res_avg_df$tool == "EpiDISH", -1][, 1:5]) + 
  theme(strip.background = element_rect(colour = "black", fill = "lightgray", linetype = "solid")) +
  ggh4x::facet_nested(~ modality ~ depth)

## 3. Scatterplots (nbl tumor fractions) 
### 2.1. With replicates
res_df$expected <- as.numeric(levels(res_df$expected))[res_df$expected]
exp_frac <- unique(res_df$expected)

AddDotsAes(res_df %>% dplyr::filter(tool == "EpiDISH", modality == "rpc") %>%
             ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = replicate, shape = depth)), exp_frac) +
  ggplot2::labs(color = "Modality") +
  facet_wrap(~ depth)

### 2.2. Without replicates
res_avg_df$expected <- as.numeric(levels(res_avg_df$expected))[res_df$expected]

AddDotsAes(res_avg_df %>% dplyr::filter(tool == "EpiDISH") %>%
             ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = modality, shape = depth)), exp_frac) +
  ggplot2::labs(color = "Modality") +
  ggh4x::facet_nested(~ "10nbl22normal" + depth)

# --- Benchmark ----
## Adapt each plot code according to your aims!

## 1. All coefficients vs all errors
### 1.1. Creating errors and coefficients dfs
bench_corr_df <- reshape2::melt(bench_metrics_df[, c(1:5)])
bench_corr_df$depth <- factor(bench_corr_df$depth, levels = unique(bench_corr_df$depth))
bench_corr_df$replicate <- factor(bench_corr_df$replicate, levels = unique(bench_corr_df$replicate))

bench_err_df <- reshape2::melt(bench_metrics_df[, c(1:3, 6:8)])
bench_err_df$depth <- factor(bench_err_df$depth, levels = unique(bench_err_df$depth))
bench_err_df$replicate <- factor(bench_err_df$replicate, levels = unique(bench_err_df$replicate))
colnames(bench_err_df) <- c("tool", "depth", "replicate", "variable2", "value2")

bench_corr_avg_df <- reshape2::melt(bench_metrics_avg_df[, c(1:4)])
bench_corr_avg_df$depth <- factor(bench_corr_avg_df$depth, levels = unique(bench_corr_avg_df$depth))

bench_err_avg_df <- reshape2::melt(bench_metrics_avg_df[, c(1:2, 5:7)])
bench_err_avg_df$depth <- factor(bench_err_avg_df$depth, levels = unique(bench_err_avg_df$depth))
colnames(bench_err_avg_df) <- c("tool", "depth", "variable2", "value2")

### 1.2. Plotting
AddMetricsAes(merge(bench_corr_df, bench_err_df) %>% 
                ggplot2::ggplot(aes(x = depth, y = tool, color = value, size = 1/value2))) + 
  labs(y = "Deconvolution tool", size = "1/error metric") +
  ggh4x::facet_nested(~ variable2 ~ replicate ~ "10nbl22normal" + variable) 

AddMetricsAes(merge(bench_corr_avg_df, bench_err_avg_df) %>% 
                ggplot2::ggplot(aes(x = depth, y = tool, color = value, size = 1/value2))) + 
  labs(y = "Deconvolution tool", size = "1/error metric") +
  ggh4x::facet_nested(~ variable2 ~ "10nbl22normal" + variable) 

## 2. Barplots
### 2.1. With replicates
bench_df$expected <- factor(bench_df$expected, levels = unique(bench_df$expected))
bench_df$depth <- factor(bench_df$depth, levels = gtools::mixedsort(unique(bench_df$depth)))
bench_df$replicate <- factor(bench_df$replicate, levels = gtools::mixedsort(unique(bench_df$replicate))) 

PlotBars(bench_df) +
  ggh4x::facet_nested(~ replicate ~ depth ~ "10nbl22normal" + tool)

### 2.2. Without replicates
bench_avg_df$expected <- factor(bench_avg_df$expected, levels = unique(bench_avg_df$expected))
bench_avg_df$depth <- factor(bench_avg_df$depth, levels = gtools::mixedsort(unique(bench_avg_df$depth)))

PlotBars(bench_avg_df) +
  ggh4x::facet_nested(~ tool ~ "10nbl22normal" + depth)

## 3. Scatterplot (nbl tumor fraction)
### 3.1. With replicates
bench_df$expected <- as.numeric(levels(bench_df$expected))[bench_df$expected]

AddDotsAes(bench_df %>% ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool, shape = depth)), exp_frac) +
  ggplot2::labs(color = "Deconvolution tool") +
  ggh4x::facet_nested(~replicate ~ "10nbl22normal" + depth)

### 3.2. Without replicates
bench_avg_df$expected <- as.numeric(levels(bench_avg_df$expected))[bench_avg_df$expected]

AddDotsAes(bench_avg_df %>% ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool, shape = depth)), exp_frac) +
  ggplot2::labs(color = "Deconvolution tool") +
  ggh4x::facet_nested(~ "10nbl22normal" + depth)

# --- Other (special) inter-comparative investigations ----
## Adapt each function (setting different variables, parameters, filtering) according to your aims!

## 1. CBS algorithms
### 1.1. With replicates
cbs_bench_df <- rbind(bench_df %>% dplyr::filter(tool == "CIBERSORT"), 
                      res_df[1:7] %>% dplyr::filter(tool == "EpiDISH", modality == "cbs") %>% dplyr::select(-modality))
cbs_bench_df$depth <- factor(cbs_bench_df$depth, levels = gtools::mixedsort(unique(cbs_bench_df$depth)))

PlotBars(cbs_bench_df) +
  ggh4x::facet_nested(~ depth ~ replicate ~ "10nbl22normal" + tool)

### 1.2. Without replicates
cbs_bench_avg_df <- rbind(bench_avg_df %>% dplyr::filter(tool == "CIBERSORT"), 
                          res_avg_df[1:6] %>% dplyr::filter(tool == "EpiDISH", modality == "cbs") %>% dplyr::select(-modality))
cbs_bench_avg_df$depth <- factor(cbs_bench_avg_df$depth, levels = gtools::mixedsort(unique(cbs_bench_avg_df$depth)))

PlotBars(cbs_bench_avg_df) +
  ggh4x::facet_nested(~ depth ~ "10nbl22normal" + tool)

## 2. DMR-selection algorithm
### 2.1. How many DMRs were retrieved by both tools?
mcbs_dmrs <- read.table("data/test_beta_10NGSnbl22NGSnormal_filter30_allClusters_GRCh37_0.4_1000_Signature.txt", header = 1)
pr_dmrs <- read.csv("dmrs_ref_db.csv")

same_dmrs <- intersect(pr_dmrs$X, mcbs_dmrs$NAME)
length(same_dmrs) #657, i.e. 343 exclusively find in each tool

### 2.2. Comparison between MCBS and PRMeth QP
ref_db <- read.csv("ref_db.csv", row.names = 1)
samples_m <- read.csv("samples_m.csv", row.names = 1, check.names = FALSE)

shared_dmrs_db <- as.matrix(subset(ref_db, rownames(ref_db) %in% same_dmrs))
dmrs_samples_m <- as.matrix(subset(samples_m, rownames(samples_m) %in% same_dmrs))

#### Deconvolving samples 
mcbs_res <- as.data.frame(EpiDISH::epidish(dmrs_samples_m, shared_dmrs_db, method = "CBS")$estF) %>%
  rownames_to_column("X")

mcbs_df <- FormattingSamplesName(mcbs_res) %>%
  rownames_to_column() %>%
  add_column(tool = rep("MBCS", 90)) %>%
  add_column(depth = rep(depth_num, each = 30)) %>%
  add_column(replicate = rep(rep_num, 30))

pr_res <- as.data.frame(t(qp(dmrs_samples_m, shared_dmrs_db)$H)) %>%
  rownames_to_column("X")

pr_df <- FormattingSamplesName(pr_res) %>%
  rownames_to_column() %>%
  add_column(tool = rep("PRMeth", 90)) %>%
  add_column(depth = rep(depth_num, each = 30)) %>%
  add_column(replicate = rep(rep_num, 30))

#### Creating the dmrs-df
dmrs_df <- rbind(mcbs_df, pr_df)
names(dmrs_df)[names(dmrs_df) == "rowname"] <- "expected" #editing the obtained df
dmrs_df$expected <- rep(exp_frac, each = 3, 6)
dmrs_df <- dmrs_df[, c(4,5,6,1,2,3)]
dmrs_df <- dmrs_df[order(dmrs_df$tool, as.numeric(dmrs_df$depth), dmrs_df$replicate, dmrs_df$expected), ]
write.csv(dmrs_df, "dmrs_df.csv")

all_dfs <- split(dmrs_df[order(dmrs_df$tool, as.numeric(dmrs_df$depth), dmrs_df$expected), ], rep(1:60, each = 3))
dmrs_avg_df <- cbind(dmrs_df[order(dmrs_df$tool, as.numeric(dmrs_df$depth), dmrs_df$expected), ] %>% dplyr::distinct(tool, depth, expected), 
                     AveragingPhenos(all_dfs)) %>%
  dplyr::select(-unknown)
write.csv(dmrs_avg_df, "dmrs_avg_df.csv")

#### Creating the dmrs-metrics df
dmrs_dfs <- split(dmrs_df, rep(1:18, each = 10)) 

dmrs_metrics_df <- CreateMetricsDf(dmrs_dfs) %>%
  add_column(tool = c(rep("MCBS", 9), rep("PRMeth", 9))) %>%
  add_column(depth = rep(depth_num, each = 3, 2)) %>%
  add_column(replicate = rep(rep_num, 6))

dmrs_metrics_df <- dmrs_metrics_df[, c(6,7,8,1,2,3,4,5)]
write.csv(dmrs_metrics_df, "dmrs_metrics_df.csv")

dmrs_avg_dfs <- split(dmrs_avg_df, rep(1:6, each = 10)) 

dmrs_metrics_avg_df <- CreateMetricsDf(dmrs_avg_dfs) %>%
  add_column(tool = c(rep("MCBS", 3), rep("PRMeth", 3))) %>%
  add_column(depth = rep(depth_num, 2))

dmrs_metrics_avg_df <- dmrs_metrics_avg_df[, c(6,7,1,2,3,4,5)]
write.csv(dmrs_metrics_avg_df, "dmrs_metrics_avg_df.csv")

#### Plotting
dmrs_avg_df$depth <- factor(dmrs_avg_df$depth, levels = unique(dmrs_avg_df$depth))
dmrs_df$depth <- factor(dmrs_df$depth, levels = unique(dmrs_df$depth))
dmrs_df$replicate <- factor(dmrs_df$replicate, levels = unique(dmrs_df$replicate))

PlotBars(dmrs_avg_df) +
  ggh4x::facet_nested(~ "10nbl22normal" + tool ~ depth)

PlotBars(dmrs_df) +
  ggh4x::facet_nested(~ replicate ~ depth ~ "10nbl22normal" + tool)

## 3. Inequality algorithms (1 unknown cell type)
### 3.1. With replicates
ineq_df <- rbind(res_df %>% dplyr::filter(modality == "cpIneq"), 
                 res_df %>% dplyr::filter(modality == "nmf"))
ineq_df$depth <- factor(ineq_df$depth, levels = unique(ineq_df$depth))

AddDotsAes(ineq_df %>% ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool)), exp_frac) +
  labs(color = "Deconvolution tool") +
  ggh4x::facet_nested(~ replicate ~ "10nbl22normal" + depth)

### 3.1. Without replicates
ineq_avg_df <- rbind(res_avg_df %>% dplyr::filter(modality == "cpIneq"), 
                     res_avg_df %>% dplyr::filter(modality == "nmf"))
ineq_avg_df$depth <- factor(ineq_avg_df$depth, levels = unique(ineq_avg_df$depth))

AddDotsAes(ineq_avg_df %>% ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool)), exp_frac) +
  labs(color = "Deconvolution tool") +
  ggh4x::facet_nested(~ "10nbl22normal" + depth)
