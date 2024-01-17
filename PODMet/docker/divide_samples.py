#!/usr/bin/python3

import pandas as pd
import numpy as np
import sys

def divide_samples(f, outDir):
    """the function takes as input a single .csv file containing all samples that you want to deconvolve; 
        it returns .csv files for each samples in the specified output directory."""
    
    samples_df = pd.read_csv(f, sep = ",", header = 0, index_col = 0)

    for col in samples_df.columns:
        sample_df = samples_df[col]
        sample_df.to_csv(outDir + col + ".csv", index=True, header=True)
    
    return("Samples division into single .csv files: done!")


if __name__ == "__main__":
    f = sys.argv[1]
    outDir = sys.argv[2]
    print(divide_samples(f, outDir))