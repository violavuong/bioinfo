#!/usr/bin/Rscript

# Loading libraries
suppressWarnings(suppressMessages(library(quadprog)))

source(paste("PRMeth/PRMeth/R/prmeth.R"))

args <- commandArgs(trailingOnly = TRUE)

# Five arguments must be provided; an error message (and help) will be returned otherwise
if (args[1] %in% c("-h", "--help")) {
  message("The following 5 arguments must be supplied: \n
  \t path-to-reference-file
  \t path-to-samples-directory
  \t PRMeth deconvolution modality [RF, NMF, QP]
  \t path-to-outputs-directory \n")
  q()
} else if (length(args) == 4) {
  ref_path <- file.path(args[1]) #defining the reference path
  samplesDir <- paste(args[2]) #defining the samples path
  mod <- args[3] #defining PRMeth deconvolution modality
  outDir <- paste(args[4]) #defining the outputs directory path
} else if (length(args) != 4) {
  stop("Four arguments must be supplied. For more info, run the following command: Rscript PRMeth.R -h/--help")
  q()
}

# Reference dataset
# One of the output file of the deconvolution tool MethAtlas - already collapsed, filtered, NAs-deprived
message("Formatting the reference dataset: initiated!")

ref_df <- read.table(ref_path, sep = ",", header = TRUE, quote = "", dec = ".", na.strings = c("NA", ""))
ref_df[is.na(ref_df)] <- 0
rownames(ref_df) <- ref_df$IlmnID
ref_m <- as.matrix(ref_df[,-1])

# Counting how many entities/phenos we have in the reference
n <- length(unique(colnames(ref_m)))

message("Formatting the reference dataset: complete!")

# Samples dataset
# Another output file of the deconvolution tool MethAtlas
message("Formatting the samples dataset: initiated!")

# Appending file per time to create a single file
samples <- list.files(path = samplesDir, pattern = "*.csv", full.names = TRUE, recursive = FALSE)

for (sample in samples) {
  if (sample == samples[1]) {
    samples_df <- read.csv(sample, header = TRUE, quote = "", dec = ".", na.strings = c("NA", ""))
  } else { # otherwise, append samples together
    tmp_df <- read.csv(sample, header = TRUE, quote = "", dec = ".", na.strings = c("NA", ""))
    samples_df <- cbind(samples_df, tmp_df)
    samples_df <- samples_df[,-(ncol(samples_df) - 1)]
    rm(tmp_df)
  }
}

samples_df[is.na(samples_df)] <- 0
rownames(samples_df) <- samples_df$IlmnID
samples_m <- as.matrix(samples_df[,-1])
samples_m <- subset(samples_m, rownames(samples_m) %in% rownames(ref_m))

# Creating an output file name
#sample_pref <- strsplit(list.files()[1], "_\\s*(?=[^/]+$)", perl = TRUE)
#samplesFileName <- stri_paste(c("test", tail(unlist(sample_pref), n = 2)), collapse = "_")

write.csv(samples_m, file = paste0(outDir, "samples.csv"), quote = FALSE)

message("Formatting the samples dataset: complete!")

# Running PRMeth
message("Running PRMeth deconvolution algorithm. PRMeth can be run in 3 modalities - RF, NMF, QP.")
message("PRMeth deconvolution: started!")

if (mod == "RF") {
  res <- t(rf(samples_m, n)$H)
  output <- "/res_rf.csv"
} else if (mod == "NMF") {
  res <- t(prmeth(samples_m, ref_m, n + 1)$H) # assuming one unknown entity
  output <- "/res_nmf.csv"
} else if (mod == "QP") {
  res <- t(qp(samples_m, ref_m)$H)
  output <- "/res_qp.csv"
}
write.csv(res, file = paste0(outDir, output), quote = FALSE)

message("PRMeth deconvolution: done!")
