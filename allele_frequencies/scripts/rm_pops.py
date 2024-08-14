#!/usr/bin/env python
# coding: utf-8

## rm_pops

#### Harrison Ostridge, 26/08/2022

#This removed populations from BayPass input files. This is just for follow up analyses. Note that this will mess up filteres such as requiring <30% missingness.


# Import modules
import numpy as np
import pandas as pd
import argparse

# Argparse
parser = argparse.ArgumentParser(description="Apply minimum allele count filter")
parser.add_argument("-in_prefix", metavar="Input prefix", type=str,
                    help="Full path and prefix of files which end in _pop.allele.counts_minMAC2 or _baypass.input.geno_minMAC2")
parser.add_argument("-rm_pops", metavar="populations to remove", nargs='+', 
                    help="Write out populations to remove in the standard format eg. -rm_pops e.IssaValley c.Chinko")
args = parser.parse_args()

# Rename (to make testing in jupyter notebook easier)
in_prefix=args.in_prefix
rm_pops=args.rm_pops
print(rm_pops)
# Read in 
allele_counts = pd.read_csv(in_prefix+"_pop.allele.counts_minMAC2", sep="\t")
baypass_input = pd.read_csv(in_prefix+"_baypass.input.geno_minMAC2", sep=" ", header=None)
# Get column names from population names
rm_pops_cols=list(map(lambda st: str.replace(st, ".", "_"), rm_pops))
rm_pops_cols_dac=[pop+'_chr1.dac' for pop in rm_pops_cols]
rm_pops_cols_aac=[pop+'_chr1.aac' for pop in rm_pops_cols]
rm_pops_cols=rm_pops_cols_dac+rm_pops_cols_aac
# Make masks
allele_counts_pop_mask=(allele_counts.columns.isin(rm_pops_cols)==False)
baypass_input_pop_mask=allele_counts_pop_mask[2:]
# Filter and write file
allele_counts.loc[:,allele_counts_pop_mask].to_csv(in_prefix+"_pop.allele.counts_minMAC2.rm_"+str("_".join(rm_pops)), sep="\t", header=True, index=False)
baypass_input.loc[:,baypass_input_pop_mask].to_csv(in_prefix+"_baypass.input.geno_minMAC2.rm_"+str("_".join(rm_pops)), sep=" ", header=False, index=False)
