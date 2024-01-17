#!/usr/bin/Rscript

# Loading libraries
suppressWarnings(suppressMessages(library(EpiDISH)))
suppressWarnings(suppressMessages(library(stringi)))

args <- commandArgs(trailingOnly = TRUE)

# Four arguments must be provided; an error message (and help) will be returned otherwise
if (args[1] %in% c("-h", "--help")) {
  message("The following 4 arguments must be supplied: \n
  \t path-to-reference-file
  \t path-to-samples-directory
  \t EpiDISH deconvolution modality [RPC, CBS, CP_ineq, CP_eq]
  \t path-to-outputs-directory \n")
  q()
} else if (length(args) == 4) {
  ref_path <- file.path(args[1]) #defining the reference path
  samplesDir <- file.path(args[2]) #defining the samples directory C:/Users/Dell/Documents/EpiDISH/data/samples/
  mod <- args[3] #defining EpiDISH deconvolution modality
  outDir <- paste(args[4]) #defining the outputs directory path
} else if (length(args) != 4) {
  stop("Four arguments must be supplied. For more info, run the following command: Rscript EpiDISH.R -h/--help")
  q()
}

# Reference dataset
# One of the output file of the deconvolution tool MethAtlas - already collapsed, filtered, NAs-deprived
message("Formatting the reference dataset: initiated!")

ref_df <- read.csv(ref_path, header = TRUE, quote = "", dec = ".", na.strings = c("NA", ""))
ref_df[is.na(ref_df)] <- 0
rownames(ref_df) <- ref_df$IlmnID
ref_m <- as.matrix(ref_df[,-1])

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
sample_pref <- strsplit(list.files()[1], "_\\s*(?=[^/]+$)", perl = TRUE)
samplesFileName <- stri_paste(c("test", tail(unlist(sample_pref), n = 2)), collapse = "_")

write.csv(samples_m, file = paste0(outDir, samplesFileName, ".csv"), quote = FALSE)

message("Formatting the samples dataset: complete!")

# Running EpiDISH
message("Running EpiDISH deconvolution algorithm. EpiDISH can be run in 4 modalities - RPC, CBS, CP_ineq, CP_eq.")
message("EpiDISH deconvolution: started!")

if (mod == "RPC") {
  res <- epidish(samples_m, ref_m, method = "RPC")$estF
  output <- "res_rpc.csv"
} else if (mod == "CBS") {
  res <- epidish(samples_m, ref_m, method = "CBS")$estF
  output <- "res_cbs.csv"
} else if (mod == "CP_ineq") {
  res <- epidish(samples_m, ref_m, method = "CP", constraint = "inequality")$estF
  output <- "res_cpIneq.csv"
} else if (mod == "CP_eq") {
  res <- epidish(samples_m, ref_m, method = "CP", constraint = "equality")$estF
  output <- "res_cpEq.csv"
}
write.csv(res, file = paste0(outDir, output), quote = FALSE)

message("EpiDISH deconvolution: done!")