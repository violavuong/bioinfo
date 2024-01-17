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
#   CellDMC(beta.m, pheno.v, frac.m, adjPMethod = "fdr", adjThresh = 0.05, cov.mod = NULL, sort = FALSE, mc.cores = 1)
#     A function that allows the identification of differentially methylated cell-types in EWAS, those specific cell types(s) responsible 
#     for the observed differential methylation
#
# -------------------------------------------------------------------------------------------------------------------------------------------------------

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Checking if there is at least one argument; if not, return an error message
# if (length(args) == 0) {
# stop("At least one argument must be supplied.")
# } else if (length(args) == 1) {
# args[2] = "out.txt"
# }

# Loading libraries
library(dplyr)
library(EpiDISH)
library(readr)
library(tidyverse)

# Defining and setting the working directory (wd)
wd <- paste("C:/Users/Dell/Desktop/R implementations/", sep = "/")
if (file.exists(wd)){
  setwd(file.path(wd))
} else {
  message("Creating new directory.")
  dir.create(file.path(paste("C:/Users/Dell/Desktop/R implementations/", sep = "/")), showWarnings = FALSE)
  dir.create(file.path(paste0(wd)))
  setwd(file.path(wd))
}

#### Pre-processing step: creating reference matrix, phenotype vector, and fraction matrix ####

##### 1. building the reference database: EpiDISH takes as input only matrices; NA values are also not allowed
##### the reference corresponds to train_methatlas, containing a matrix which rows correspond to samples and cols to beta values for each 
##### cg cluster present
##### 2. the cluster file is a csv/tsv file containing 4 cols, respectively chr, cg starting position, cg ending position, and cg cluster ID
##### the cluster file is not really needed, but it's useful to take track of which cg you're referring to

message("Retrieving the reference set and cluster file.")

cluster_path <- file.path("data/RRBS_450k_intersectClusters.tsv")
ref_path <- file.path("data/train_methatlas_10NGSnbl22NGSnormal_allClusters_GRCh37") 

cluster_df <- read.table(cluster_path, sep = "\t", header = TRUE)
ref_m <- t(read.table(ref_path, sep = "\t", header = FALSE, quote = "", dec = ".", na.strings = c("NA", ""))) #transposing to have cg as rows
ref_m[is.na(ref_m)] <- 0 #changing all NAs with 0s

##### creating the phenotype vector (it corresponds to your entities - nbl and normal): it corresponds to the header of the reference matrix
ref_pheno <- ref_m[1, ]

##### changing reference rows and cols names
ref_m <- ref_m[-1,]
colnames(ref_m) <- paste(gsub("-", "_", ref_pheno), c(1:ncol(ref_m)), sep = "_")
rownames(ref_m) <- paste0("cg", cluster_df$clusterID)
class(ref_m) <- "numeric"

##### creating the fraction matrix - it contains the the fractions of each entity, in this case even if it's not realistic we are going for 
##### 0 and 1 respectively for normal and nbl
frac_m <- matrix(data = NA, nrow = dim(ref_m)[2], ncol = 2)
colnames(frac_m) <- c("nbl", "normal")
rownames(frac_m) <- paste(colnames(ref_m))

for (i in 1:dim(frac_m)[1]) {
  row_name <- rownames(frac_m)[i]
  if (startsWith(row_name, "nbl")) {
    frac_m[i, "normal"] <- runif(1, 0.000001, 0.000009)
    frac_m[i, "nbl"] <- 1 - frac_m[i, "normal"] 
  } else if (startsWith(row_name, "normal")) {
    frac_m[i, "nbl"] <- runif(1, 0.000001, 0.000009)
    frac_m[i, "normal"] <- 1 - frac_m[i, "nbl"] 
  }
  else {stop("Invalid option.")}
}

#### Performing the DMR selection using CellDMC function ####
dmrs <- CellDMC(ref_m, ref_pheno, frac_m)

#### Comparing the DMRs of EpiDISH and MethylCBS ####

#### retrieving the DMRs cg names of EpiDISH
dmrs_epi <- dmrs$dmct %>%
  as.data.frame() %>%
  filter(DMC == 1) %>%
  rownames()
length(dmrs_epi)

#### retrieving the DMRs cg names of MethylCIBERSORT (corresponding to metcbs$NAME)
metcsb <- read_table("data/In silico/MethylCBS_1000_signatures.txt")

#### retrieving the common DMRs shared by both methodologies
same_dmrs <- intersect(dmrs_epi, metcsb$NAME) 
length(same_dmrs)

#### with respect to the shared cgs in both methodologies, DM cgs present exclusively in each methodology are retrieved
only_metcbs_dmrs <- subset(metcsb$NAME, ! metcsb$NAME %in% same_dmrs)
only_epi_dmrs <- subset(dmrs_epi, ! dmrs_epi %in% same_dmrs)
length(only_epi_dmrs)
length(only_metcbs_dmrs)

#### retrieving the distributions (# of hypo- and hyper-methylated)
##### shared cgs (between EpiDISH and MethylCBS)
dmrs$dmct %>% 
  as.data.frame() %>%
  filter(DMC == 1) %>%
  subset(., rownames(.) %in% same_dmrs) %>%
  count(nbl, normal)

##### only-EpiDISH cgs 
dmrs$dmct %>% 
  as.data.frame() %>%
  filter(DMC == 1) %>%
  subset(., ! rownames(.) %in% same_dmrs)  %>%
  count(nbl, normal)