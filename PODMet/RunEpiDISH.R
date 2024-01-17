# EpiDISH R script v.2.16.0
# Description: EpiDISH is an R package to infer the proportions of a priori known cell-types present in a sample representing a mixture of
# such cell-types. EpiDISH can be used on DNAm data of whole blood, generic epithelial tissue and breast tissue. It also possible to 
# identify differentially methylated cell-types (DMC) and their directionality of change in EWAS. 
# Author: Andrew E. Teschendorff, Shijie C. Zheng (shijieczheng@gmail.com)
#
# Dependencies: 
#   R v4.1 or later
#   install.packages('MASS')
#   install.packages('e1071')
#   install.packages('quadprog')
#   install.packages('parallel')
#   install.packages('stats')
#   install.packages('matrixStats')
#   install.packages('stringr')
#   install.packages('locfdr')
#   install.packages('Matrix')
# 
# Functions:
#   epidish(beta.m, ref.m, method=c("RPC", "CBS", "CP"), maxit=50, nu.v=c(0.25, 0.5, 0.75), constraint=c("inequality", "equality"))
#     A reference-based function to infer the fractions of a priori known cell subtypes present in a sample representing a a mixture of 
#     such cell-types.
# ---------------------------------------------------------------------------------------------------------------------------------------------

#!/usr/bin/env Rscript

library(tidyverse)

# Loading complementary scripts
source(paste("C:/Users/Dell/Desktop/R implementations/functions/DataManipulation.R"))

#### Preparing the R environment for the deconvolution ####
wd <- paste("C:/Users/Dell/Desktop/R implementations/9nbl22normal/9vitros/100hyper_30/", sep = "/")
setwd(wd)

#### Reference ####
# EpiDISH takes as input only matrices; NA values are also not allowed. Refer to MethAtlas output file train_methatlas_*NGSnbl22NGSnormal_allClusters_GRCh37 
# A cluster file is needed, a .csv/.tsv file containing 4 cols, respectively chr, cg starting position, cg ending position, and cg cluster ID
cluster_path <- file.path("data/RRBS_450k_intersectClusters.tsv")
ref_path <- file.path("data/train_methatlas_9NGSnbl22NGSnormal_allClusters_GRCh37") 

cluster_df <- read.table(cluster_path, sep = "\t", header = TRUE)
ref_m <- t(read.table(ref_path, sep = "\t", header = FALSE, quote = "", dec = ".", na.strings = c("NA", ""))) #transposing to have cg as rows
ref_m[is.na(ref_m)] <- 0 #changing all NAs with 0s
ref_db <- CreateRef(MatchingCGCluster(ref_m, cluster_df))
write.csv(ref_db, file = "ref_db.csv")

#### Samples ####
# EpiDISH takes as input only matrices; NA values are also not allowed. Refer to MethAtlas output file test_beta_*NGSnbl22NGSnormal_filter30_allClusters_GRCh37
samples_path <- file.path("data/test_beta_9NGSnbl22NGSnormal_filter30_allClusters_GRCh37")
samples_tmp_m <- t(read.table(samples_path, sep = "\t", header = FALSE, quote = "", dec = ".", na.strings = c("NA", "")))
samples_tmp_m[is.na(samples_tmp_m)] <- 0
samples_tmp_m <- MatchingCGCluster(samples_tmp_m, cluster_df)

# Retrieving only the matching cgs of ref and samples 
true_clusters <- intersect(rownames(ref_db),rownames(samples_tmp_m))
samples_m <- subset(samples_tmp_m, rownames(samples_tmp_m) %in% true_clusters) 
write.csv(samples_m, file = "samples_m.csv")

#### Running EpiDISH ####
rpc_res <- EpiDISH::epidish(samples_m, ref_db, method = "RPC")$estF
cbs_res <- EpiDISH::epidish(samples_m, ref_db, method = "CBS")$estF
cpIneq_res <- EpiDISH::epidish(samples_m, ref_db, method = "CP")$estF
cpEq_res <- EpiDISH::epidish(samples_m, ref_db, method = "CP", constraint = "equality")$estF

write.csv(rpc_res, file = "EpiDISH_rpc_res.csv")
write.csv(cbs_res, file = "EpiDISH_cbs_res.csv")
write.csv(cpIneq_res, file = "EpiDISH_cpIneq_res.csv")
write.csv(cpEq_res, file = "EpiDISH_cpEq_res.csv")
