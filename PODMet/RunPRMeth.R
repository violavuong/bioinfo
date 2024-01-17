# PRMeth R script v.?
# Description: PRMeth is an R package able to deconvolve tumor mixtures using partially available DNA methylation data. Here, an iteratively 
# optimized non-negative matrix factorization framework is adopted to simultaneously infer proportions of a priori known and unknown cell-types
# Author: Dingqin He, Ming Chen, Wenjuan Wang, Chunhui Song, Yufang Qin
#
# Dependencies: 
#   R v4.1 or later
#   install.packages('matrixStats')
#   install.packages('quadprog')
# 
# Functions:
#   select_feature(sample_m, n_start=1, n_end=1000)
#     Selects and returns the top n CpG sites with the highest cv value. Only the beta value matrix is absolutely required
#
#   getCellTypeNumber(sample_m, ref_m, k)
#     Returns the total number of cell types in the tumor mixture. The max k number of cell-types must be known a priori
#
#   prmeth(sample_m, ref_m, k, iters=500, rrsDiffStop=1e-10, lambda=1)
#     Perform the proposed deconvolution algorithm. The first three parameters have no default option. Returns the inferred tumor fraction for
#     each k cell-type in a matrix-like structure 
# ----------------------------------------------------------------------------------------------------------------------------------------------

# Loading libraries
library(matrixStats)
library(quadprog)

# Loading complementary scripts
source(paste("C:/Users/Dell/Desktop/R implementations/functions/DataManipulation.R"))
source(paste("C:/Users/Dell/Desktop/R implementations/functions/prmeth.R"))
source(paste("C:/Users/Dell/Desktop/R implementations/functions/select_feature.R"))

#### Preparing the R environment for the deconvolution ####
wd <- paste("C:/Users/Dell/Desktop/R implementations/9nbl22normal/90replicates/100hyper_30", sep = "/")
setwd(file.path(wd))

#### Reference ####
# PRMeth takes as input only matrices; NA values are also not allowed. Refer to MethAtlas output file train_methatlas_*NGSnbl22NGSnormal_allClusters_GRCh37 
# A cluster file is needed, a .csv/.tsv file containing 4 cols, respectively chr, cg starting position, cg ending position, and cg cluster ID
cluster_path <- file.path("data/RRBS_450k_intersectClusters.tsv")
ref_path <- file.path("data/train_methatlas_9NGSnbl22NGSnormal_allClusters_GRCh37") 

cluster_df <- read.table(cluster_path, sep = "\t", header = TRUE)
ref_m <- t(read.table(ref_path, sep = "\t", header = FALSE, quote = "", dec = ".", na.strings = c("NA", ""))) #transposing to have cg as rows
ref_m[is.na(ref_m)] <- 0 #changing all NAs with 0s

# Building a complete reference dataset
ref_db <- CreateRef(MatchingCGCluster(ref_m, cluster_df))
write.csv(ref_db, file = "ref_db.csv")

# Building a DMR-selected reference dataset
dmrs <- select_feature(MatchingCGCluster(ref_m, cluster_df))
dmrs_ref_db <- subset(ref_db, rownames(ref_db) %in% dmrs)
write.csv(dmrs_ref_db, file = "dmrs_ref_db.csv")

#### Samples ####
# PRMeth takes as input only matrices; NA values are also not allowed. Refer to MethAtlas output file test_beta_*NGSnbl22NGSnormal_filter30_allClusters_GRCh37
samples_path <- file.path("data/test_beta_9NGSnbl22NGSnormal_filter30_allClusters_GRCh37")
samples_tmp_m <- t(read.table(samples_path, sep = "\t", header = FALSE, quote = "", dec = ".", na.strings = c("NA", "")))
samples_tmp_m[is.na(samples_tmp_m)] <- 0
samples_tmp_m <- MatchingCGCluster(samples_tmp_m, cluster_df)

# Retrieving only the matching cgs of ref and samples 
true_clusters <- intersect(rownames(ref_db),rownames(samples_tmp_m))
dmrs_clusters <- intersect(rownames(dmrs_ref_db),rownames(samples_tmp_m))

# For the complete reference
samples_m <- subset(samples_tmp_m, rownames(samples_tmp_m) %in% true_clusters) 
write.csv(samples_m, file = "samples_m.csv")

# And for the DMRs-based one
dmrs_samples_m <- subset(samples_tmp_m, rownames(samples_tmp_m) %in% dmrs_clusters)
write.csv(dmrs_samples_m, file = "dmrs_samples_m.csv")

#### Running PRMeth ####
# First starting with the newly implemented NMF (we always assume an unknown entity)
nmf_res <- t(prmeth(samples_m, ref_db, 3)$H)
write.csv(nmf_res, file = "PRMeth_nmf_res.csv")

dmrs_nmf_res <- t(prmeth(dmrs_samples_m, dmrs_ref_db, 3)$H) 
write.csv(dmrs_nmf_res, file = "PRMeth_dmrs_nmf_res.csv")

# Then with RF (no reference)
rf_res <- t(rf(samples_m, 2)$H)
colnames(rf_res) <- c("normal", "nbl")
rf_res <- rf_res[,c(2,1)]
write.csv(rf_res, file = "PRMeth_rf_res.csv")

dmrs_rf_res <- t(rf(dmrs_samples_m, 2)$H)
colnames(dmrs_rf_res) <- c("normal", "nbl")
dmrs_rf_res <- dmrs_rf_res[,c(2,1)]
write.csv(dmrs_rf_res, file = "PRMeth_dmrs_rf_res.csv")

# Lastly, with QP (complete reference)
qp_res <- t(qp(samples_m, ref_db)$H)
write.csv(qp_res, file = "PRMeth_qp_res.csv")

dmrs_qp_res <- t(qp(dmrs_samples_m, dmrs_ref_db)$H)
write.csv(dmrs_qp_res, file = "PRMeth_dmrs_qp_res.csv")
