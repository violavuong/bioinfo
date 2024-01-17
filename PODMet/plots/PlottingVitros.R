#!/usr/bin Rscript

library(quadprog)
library(tidyverse)

# Defining the working directory 
setwd("C:/Users/Dell/Desktop/R implementations/10nbl22normal/9vitros/100top_30/")

# Loading complementary scripts
source("C:/Users/Dell/Desktop/R implementations/functions/DataManipulation.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Metrics.R")
source("C:/Users/Dell/Desktop/R implementations/functions/Plots.R")
source("C:/Users/Dell/Desktop/R implementations/functions/prmeth.R")

# --- Data prep ----
## 1. Loading deconvolution results
res_files <- list.files(pattern = "_res.csv")
res_list <- lapply(res_files, read.csv) #cbs, cpEq, cpIneq, rpc, dmrs_nmf, dmrs_qp, dmrs_rf, nmf, qp, rf

## 2. Defining important features
tools <- c("EpiDISH", "PRMeth")
mods <- c("cbs", "cpEq", "cpIneq", "rpc", "dmrs_nmf", "dmrs_qp", "dmrs_rf", "nmf", "qp", "rf")
exp_frac <-  c(0.0146, 0, 0.0399, 0.0053, 0.25, 0.0019, 1, 0.5, 0.1045) #10nbl22normal
#exp_frac <-  c(0.0146, 0, 0.0399, 0.0053, 0.0019, 0.5, 0.1045, 0.25, 1) #9nbl22normal

## 3. Creating a single every-info-df
res_df <- bind_rows(res_list) %>%
  add_column(tool = c(rep(tools[1], 36), rep(tools[2], 54))) %>%
  add_column(modality = rep(mods, each = 9))

res_df <- res_df[, c(5,6,1,2,3,4)] #re-ordering cols
names(res_df)[names(res_df) == 'X'] <- "expected"
names(res_df)[names(res_df) == 'X.1'] <- "unknown"
res_df$expected <- rep(exp_frac, 10)
res_df <- res_df[order(res_df$tool, res_df$modality, res_df$expected), ] #sorting values based on own preference
write.csv(res_df, "res_df.csv")
#str(res_df) #for future refs

## 4. Creating an errors df
all_dfs <- split(res_df, rep(1:10, each = 9)) #list of all possible dfs

errors_df <- CreateMetricsDf(all_dfs) %>% #creating the metric df
  add_column(tool = c(rep(tools[1], 4), rep(tools[2], 6))) %>%
  add_column(modality = mods)

errors_df <- errors_df[, c(6,7,1,2,3,4,5)] #changing the cols order
write.csv(errors_df, "errors_df.csv")
#str(errors_df) #for future refs, in case you need to change data type while plotting

## 5. Creating the benchmark df
bench_tools <- c("CIBERSORT", "EpiDISH", "MethAtlas", "PRMeth")
bench_exp_frac <- c(0.5, 0.25, 0.1045, 0.0399, 0.0146, 0.0053, 0.0019, 0, 1)

### CIBERSORT
cbs_df <- read.table("data/results_2perm_test_CIBERSORT_10NGSnbl22NGSnormal_filter30_100top_GRCh37_iter1.txt", header = 1)[, 1:3] %>%
  add_column(tool = rep(bench_tools[1], 9))

names(cbs_df)[names(cbs_df) == 'Mixture'] <- "expected"
cbs_df$expected <- bench_exp_frac

### EpiDISH (RPC)
epi_df <- res_df %>% 
  dplyr::filter(tool == "EpiDISH", modality == "rpc") %>%
  dplyr::select(-modality, -unknown)

### MethAtlas
matlas_df <- as.data.frame(t(read.csv("data/test_methatlas_10NGSnbl22NGSnormal_filter30_100top_GRCh37_iter1_deconv_output.csv", row.names = 1, check.names = FALSE))) %>%
  rownames_to_column("expected") %>%
  add_column(tool = rep(bench_tools[3], 9))

matlas_df$expected <- bench_exp_frac

### PRMeth (QP) - complete reference available
pr_df <- res_df %>%  
  dplyr::filter(tool == "PRMeth", modality == "qp") %>%
  dplyr::select(-modality, -unknown)

bench_df <- bind_rows(cbs_df, epi_df, matlas_df, pr_df) #creating the df
bench_df <- bench_df[, c(4,1,2,3)] #re-ordering columns
bench_df <- bench_df[order(bench_df$tool, bench_df$expected), ] #re-ordering values
write.csv(bench_df, "bench_df.csv")

## 6. Creating the benchmark error df
all_dfs <- split(bench_df, rep(1:4, each = 9)) #list of all possible dfs

bench_errors_df <- CreateMetricsDf(all_dfs) %>% #creating the error df
  add_column(tool = rep(bench_tools, each = 1))

bench_errors_df <- bench_errors_df[, c(6,1,2,3,4,5)] #changing the cols order
write.csv(bench_errors_df, "bench_errors_df.csv")
#str(bench_errors_df) #for future refs, in case you need to change data type while plotting

# --- Loading ----
res_df <- read.csv("res_df.csv", row.names = 1)
bench_df <- read.csv("bench_df.csv", row.names = 1)

metrics_df <- read.csv("metrics_df.csv", row.names = 1)
bench_metrics_df <- read.csv("bench_metrics_df.csv", row.names = 1)

# --- Intra-comparison between methods  ----
## Adapt each function (setting different variables, parameters, filtering) according to your aims!

## 1. Barplots
PlotBars(res_df[1:5] %>% dplyr::filter(tool == "EpiDISH")) +
  facet_wrap(~ modality)

## 2. Heatmaps
PlotHeatmap(res_df[1:5] %>% dplyr::filter(tool == "EpiDISH"))

## 3. Scatterplot (nbl tumor fractions) 
AddDotsAes(res_df[1:5] %>% dplyr::filter(tool == "EpiDISH") %>%
             ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = modality)), exp_frac) + 
                labs(color = "Modality") +
                facet_wrap(~ "10nbl22normal")

## 4. PCC plot
PlotPCC(res_df %>% dplyr::filter(tool == "EpiDISH"), exp_frac) +
  facet_wrap(~ modality) 

## 5. Histograms for error metrics
PlotErrorsAsBars(errors_df %>% dplyr::filter(tool == "EpiDISH") %>% dplyr::select(-PCC, -R2)) +
  facet_wrap(~ modality)

## 6. Dot plots for error metrics
corr_df <- reshape2::melt(errors_df[, c(1:4)])
err_df <- reshape2::melt(errors_df[, c(1:2, 5:7)])
colnames(err_df) <- c("tool", "modality","variable2", "value2")

AddMetricsAes(merge(corr_df, err_df) %>% dplyr::filter(tool == "EpiDISH") %>% 
                ggplot2::ggplot(aes(x = variable2, y = modality, color = value, size = 1/value2))) + 
                  labs(y = "Modality", x = "Error metric", size = "1/error metric") +
                  ggh4x::facet_nested(~ "EpiDISH" + variable)

# --- Benchmark ----
## Adapt each function (setting different variables, parameters, filtering) according to your aims!

## 1. Barplots
PlotBars(bench_df) +
  facet_wrap(~ tool)

## 2. Benchmark nbl dot plot
AddDotsAes(bench_df %>%
             ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool)), exp_frac) +
                labs(color = "Deconvolution tool") +
                facet_wrap(~ "10nbl22normal")

## 3. Benchmark PCC plot 
PlotPCC(bench_df, exp_frac) +
  facet_wrap(~ tool)

## 4. Benchmark histograms
PlotErrorsAsBars(bench_errors_df %>% dplyr::select(-PCC, -R2)) +
  facet_wrap(~ tool)

# --- Other (special) inter-comparative investigations ----
## Adapt each function (setting different variables, parameters, filtering) according to your aims!

## 1. CBS algorithms 
cbs_bench_df <- rbind(bench_df %>% dplyr::filter(tool == "CIBERSORT"), 
                      res_df[1:5] %>% dplyr::filter(tool == "EpiDISH", modality == "cbs") %>% dplyr::select(-modality))

PlotBars(cbs_bench_df) +
  ggh4x::facet_nested(~ "10nbl22normal" + tool)

## 2. DMR-selection algorithms
### 2.1. How many DMRs were retrieved by both tools?
mcbs_dmrs <- read.table("data/test_beta_10NGSnbl22NGSnormal_filter30_allClusters_GRCh37_0.4_1000_Signature.txt", header = 1)
pr_dmrs <- read.csv("dmrs_ref_db.csv")

same_dmrs <- intersect(pr_dmrs$X, mcbs_dmrs$NAME)
length(same_dmrs) #657, i.e. 343 exclusively find in each tool

### 2.2. Comparison between MCBS and PRMeth QP
ref_db <- read.csv("ref_db.csv", row.names = 1)
samples_m <- read.csv("samples_m.csv", row.names = 1)

shared_dmrs_db <- as.matrix(subset(ref_db, rownames(ref_db) %in% same_dmrs))
dmrs_samples_m <- as.matrix(subset(samples_m, rownames(samples_m) %in% same_dmrs))

#### Deconvolving samples 
mcbs_res <- EpiDISH::epidish(dmrs_samples_m, shared_dmrs_db, method = "CBS")$estF
pr_res <- t(qp(dmrs_samples_m, shared_dmrs_db)$H)

dmrs_exp <- c(0.0146, 0, 0.0399, 0.0053, 0.25, 0.0019, 1, 0.5, 0.1045)

#### Creating the dmrs-df
dmrs_df <- as.data.frame(rbind(mcbs_res, pr_res)) %>%
  rownames_to_column() %>%
  add_column(tool = c(rep("MCBS", 9), rep("PRMeth", 9)))

names(dmrs_df)[names(dmrs_df) == "rowname"] <- "expected" #editing the obtained df
dmrs_df$expected <- rep(dmrs_exp, 2)
dmrs_df <- dmrs_df[, c(4,1,2,3)]
dmrs_df <- dmrs_df[order(dmrs_df$tool, dmrs_df$expected), ]
write.csv(dmrs_df, "dmrs_df.csv")

#### Creating the dmrs-metrics df
dmrs_dfs <- split(dmrs_df, rep(1:2, each = 9)) 
dmrs_metrics_df <- CreateMetricsDf(dmrs_dfs) %>%
  add_column(tool = c("MCBS", "PRMeth"))
dmrs_metrics_df <- dmrs_metrics_df[, c(6,1,2,3,4,5)]
write.csv(dmrs_metrics_df, "dmrs_metrics_df.csv")

#### Plotting
PlotBars(dmrs_df) +
  ggh4x::facet_nested(~ "10nbl22normal" + tool)

## 3. Inequality algorithms (1 unknown cell type)
ineq_df <- rbind(res_df %>% dplyr::filter(modality == "cpIneq"), 
                 res_df %>% dplyr::filter(modality == "nmf"))

exp_frac <- unique(res_df$expected)

AddDotsAes(ineq_df %>% ggplot2::ggplot(aes(x = expected + 0.000001, y = nbl + 0.000001, color = tool)), exp_frac) +
  labs(color = "Deconvolution tool") +
  facet_wrap(~ "9nbl22normal")