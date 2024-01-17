# Data structure manipulation v.0.1
# Description: an R script containing a series of R functions (listed alphabetically) that allow data structure manipulation in order to 
# accommodate and simplify the downstream deconvolution analysis.
#
# Functions:
#   AveragingPhenos(list): creates a df with the averaged phenotypes for each sample. 
#   CreateBenchDf(df, vector, vector, vector): creates a ordered benchmark df.
#   CreateAvgBenchDf(df, vector, vector, vector): creates a ordered benchmark df which phenotypes have been averaged.
#   CreateRef(matrix): returns the reference database (as matrix) that the selected tool will perform the deconvolution on.
#   FormattingSamplesName(df): returns a dataframe which rownames has been modified and ordered. 
#   MatchingCGCluster(matrix, df): returns a matrix with clustered cgs.
# ------------------------------------------------------------------------------------------------------------------------------------------------------

#!/usr/bin/env R script

AveragingPhenos <- function(list) {
  ### list: a list of df which phenotypes expressions has to be averaged by replicates
  ### Returns the mentioned df, subsetted
  
  return(do.call(rbind, Map(data.frame, 
                            nbl = lapply(seq_along(list), function(x) mean(list[[x]]$nbl)), 
                            normal = lapply(seq_along(list), function(x) mean(list[[x]]$normal)), 
                            unknown = lapply(seq_along(list), function(x) mean(list[[x]]$unknown))))
         )
}

CreateBenchDf <- function(df, tool = c("CIBERSORT", "MethAtlas"), depth_num, rep_num) {
  ### df: a benchmark dataframe containing the deconvolution results
  ### bench_tools: a char vector containing the benchmark tools (ordered alphabetically)
  ### depth_num: a char vector containing all possible depths
  ### Returns the mentioned df in standard common format for downstream analysis
  
  df <- FormattingSamplesName(df) %>%
    rownames_to_column("expected") %>%
    add_column(tool = rep(tool, 90)) %>%
    add_column(depth = rep(depth_num, each=30)) %>%
    add_column(replicate = rep(rep_num, 30))
  
  df$expected <- as.numeric(mgsub::mgsub(df$expected , c("\\d+M_", "_\\d+$"), c("", ""))) #obtaining expected fractions
  df$expected <- df$expected/100 #converting to fractions
  df <- df[order(df$tool, as.numeric(df$depth), df$replicate, df$expected), ]
  return(df[, c(4,5,6,1,2,3) ]) #ordering the df
}

CreateAvgBenchDf <- function(df, tool = c("CIBERSORT", "MethAtlas"), depth_num, rep_num) {
  ### df: a benchmark dataframe containing the deconvolution results
  ### bench_tools: a char vector containing the benchmark tools (ordered alphabetically)
  ### depth_num: a char vector containing all possible depths
  ### rep_num: a char vector containing all possible replicates number
  ### Returns the mentioned df in standard common format for downstream analysis after averaging phenotype expression by replicates 
  
  df <- FormattingSamplesName(df) %>%
    rownames_to_column("expected") %>%
    add_column(tool = rep(tool, 90)) %>%
    add_column(depth = rep(depth_num, each=30)) %>%
    add_column(replicate = rep(rep_num, 30))
  
  df$expected <- as.numeric(mgsub::mgsub(df$expected , c("\\d+M_", "_\\d+$"), c("", ""))) #obtaining expected fractions
  df$expected <- df$expected/100 #converting to fractions
  
  rep_list <- df[order(df$expected, df$depth), ] %>%
    split(rep(1:30, each = 3))
  
  df <- cbind(df[order(df$expected, df$depth), ] %>% distinct(tool, depth, expected), #creating the new "averaged" df
                     AveragingPhenos(rep_list)) %>%
    dplyr::select(-unknown)
  
  return(df[order(df$tool, as.numeric(df$depth)), ]) #ordering the df
}

CreateRef <- function(ref_m) {
  ### ref_m: a matrix with cg as rows and entities/samples as cols
  ### Returns a matrix in a standard format, by collapsing entities (columns) into a single one by computing the mean of the beta values 
  ### belonging to the same cluster
  
  nbl_cols <- dplyr::select(data.frame(ref_m), tidyselect::starts_with("nbl")) 
  norm_cols <- dplyr::select(data.frame(ref_m), tidyselect::starts_with("normal"))
  ref_db <- data.frame(ref_m) %>%
    dplyr::mutate(nbl = rowMeans(nbl_cols), normal = rowMeans(norm_cols)) %>%
    dplyr::select(nbl, normal)
  ref_db <- as.matrix(ref_db)
 
  return(ref_db)
}

CreateMetricsDf <- function(list) {
  ### list: a list of all possible dfs
  ### Returns a df containing the computed metrics
  
  pcc <- lapply(seq_along(list), function(x) ComputeMetrics(list[[x]] %>% dplyr::select(expected, nbl), metric = "pcc"))
  r2 <- lapply(seq_along(list), function(x) ComputeMetrics(list[[x]] %>% dplyr::select(expected, nbl), metric = "r2")) 

  aad <- lapply(seq_along(list), function(x) ComputeMetrics(list[[x]] %>% dplyr::select(expected, nbl), metric = "aad"))
  medae <- lapply(seq_along(list), function(x) ComputeMetrics(list[[x]] %>% dplyr::select(expected, nbl), metric = "medae"))
  rmse <- lapply(seq_along(list), function(x) ComputeMetrics(list[[x]] %>% dplyr::select(expected, nbl), metric = "rmse"))
  
  return(do.call(rbind, Map(data.frame, PCC = pcc, R2 = r2, AAD = aad, MedAE = medae, RMSE = rmse)))
}

FormattingSamplesName <- function(df) {
  ### df: a dataframe containing the deconvolution result of interest
  ### Returns a dataframe which rownames have been modified and ordered to simplify downstream analysis
  
  df$X <- mgsub::mgsub(df$X, 
                       c("^([0-9M]+)_.+_([A-Za-z]+)_([0-9])[-\\._]([0-9]+)_[a-z]+(\\d+)_R1.*", "^([0-9M]+)_.+_([A-Za-z]+)_([0-9]+)_[a-z]+(\\d+)_R1.*"),
                       c("\\1_\\3.\\4_\\5", "\\1_\\3_\\4"))
  rownames(df) <- df$X
  df <- df[-1]
  return(df[gtools::mixedsort(rownames(df)), ])
}

MatchingCGCluster <- function(data_m, cluster_df) {
  ### data_m: a matrix-like data structure
  ### cluster_df: a df-like structure containing the genomic position of each cg, along side the corresponding cg cluster ID
  ### Returns a matrix with each cg associated to its corresponding cluster
  
  data_tmp <- data_m[-1,]
  colnames(data_tmp) <- paste(gsub("-", "_", data_m[1,]), c(1:ncol(data_tmp)), sep = "_")
  rownames(data_tmp) <- paste0("cg", cluster_df$clusterID)
  class(data_tmp) <- "numeric"
  
  return(data_tmp)
}