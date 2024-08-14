# subset_baypass_input-22-07-20.py

## Harrison Ostridge 20/07/2022

## This subsets BayPass input by selecting every other SNP. 
## This is done becuase running the AUX model on all sites at once takes too long for myriad (>48 hours). 
## We choose to select every other SNP to ensure each subset has a good sample across the whole genome.

# Import modules
import numpy as np
import pandas as pd
import argparse

# Argparse
parser = argparse.ArgumentParser(description="Subset BayPass input.")
parser.add_argument("-baypass_input", metavar="Input file", type=str,
                    help="Full path to _baypass.input.geno file")
parser.add_argument("-allele_counts", metavar="Input file", type=str,
                    help="Full path to _pop.allele.counts file")
parser.add_argument("-n_subsets", metavar="Number of subsets", type=int)

args = parser.parse_args()

# Rename arguments (allows easy testing of script in jupyter)
baypass_input_file=args.baypass_input
allele_counts_file=args.allele_counts
n_subsets=args.n_subsets

# Read in 
baypass_input = pd.read_csv(baypass_input_file, sep=" ", header=None)
allele_counts = pd.read_csv(allele_counts_file, sep="\t")

# Input checks
## Same number of SNPs in both files
assert baypass_input.shape[0]==allele_counts.shape[0], "Number of SNPs not equal in two files."
## Same number of populations in both files
### -2 from allele counts file becuase it has chr and pos columns also
assert baypass_input.shape[1]==(allele_counts.shape[1]-2), "Number of populations not equal in two files."
baypass_input_subsets={}
allele_counts_subsets={}

for i in range(0, n_subsets):
    baypass_input.iloc[i::n_subsets, :].to_csv(baypass_input_file+"_subset"+str(i+1)+"of"+str(n_subsets), sep=" ", header=False, index=False)
    allele_counts.iloc[i::n_subsets, :].to_csv(allele_counts_file+"_subset"+str(i+1)+"of"+str(n_subsets), sep="\t", header=True, index=False)
    
