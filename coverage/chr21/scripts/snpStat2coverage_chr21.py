#!/usr/bin/env python
# coding: utf-8

# snpStat2coverage_chr21.py

# Harrison Ostridge, 26/11/2021

# Prior to this script, ANGSD should be run chromosome by chromosome and population by population
## to estimate population allele frequencies with -snpStat 1 to output a snpStat.gz file.
## maf2dac must then be run to identify snps of interest.

print("Importing modules")

import numpy as np
import pandas as pd
import os
from functools import reduce
import gc
import matplotlib.pyplot as plt
import sys
import copy
import argparse

# Argparse
parser = argparse.ArgumentParser(description="Combine ANGSD .maf.gz files from multiple populations into a allele counts file for input into BayPass.")
parser.add_argument("-in_dir", metavar="Input directory", type=str,
                    help="Input directory containing .snpStat.gz files.", default='./')
parser.add_argument("-snp_by_pop_file", metavar="pop.allele.counts file.", type=str,
                    help="Path to file with columns 'chr' and 'pos' defining the sites o be investigated (this may be a _pop.allele.counts file outputted from maf2dac.")
parser.add_argument("-out_prefix", metavar="Output file prefix", type=str, default='snpStat2coverage.out',
                    help="The chromosome, subspecies and missing populations threshold will be added as suffixes to the file name.")
args = parser.parse_args()

# Modify/rename parameters
input_dir=args.in_dir
snp_by_pop_file=args.snp_by_pop_file
output_prefix=args.out_prefix

# Read in data

## Sites and populations
snp_by_pop=pd.read_csv(snp_by_pop_file, sep="\t")
#snp_by_pop=snp_by_pop[snp_by_pop['chr']=='chr'+chromo]
### Select derived alelel column names 
pops = pd.DataFrame(snp_by_pop.columns[snp_by_pop.columns.str.endswith("dac")]).fillna(0)
### Remove suffix
pops = list(pops[0].str.split('_chr21.dac').str[0])
### Replace '_' with '.' as this is what is used in the file prefixes
pops_file_pre=list([p.replace('_', '.') for p in pops])
pops_file_pre

## snpStats

snpstat = {}
for pop in pops_file_pre:
    # Get list of files corresponding to a population (for all chrs)
    filelist = sorted([os.path.join(root, name)
                for root, dirs, files in os.walk(input_dir)
                for name in files
                if name.startswith(pop) 
                           & name.endswith(".snpStat.gz")])
    # Read in all chr files for a population and combine into one file
    dfs = [pd.read_csv(file, sep=" |\t") for file in filelist]
    df = pd.concat(dfs, axis=0)
    df.rename(columns={'Chromo': 'chr', 'Position': 'pos'}, inplace=True)
    # Keep onlt selected sites using inner merge
    df=pd.merge(snp_by_pop[['chr','pos']],df,on=['chr', 'pos'], how='inner')
    # Replace '.' with '_' for naming 
    pop_name = pop.replace(".", "_")
    # Create a coverage column by summing the major and minor allele counts on the forward and reverse strand
    df[pop_name+'_coverage']= df.iloc[:, 2:5].sum(axis=1)
    df=df.loc[:, ['chr', 'pos', pop_name+'_coverage']]
    # Save to dictionary
    snpstat[pop_name] = df

# Merge populations
def merge_pops(snpstat):
    """Combine dataframes from multiple populations into one large dataframe"""
    print("Merging populations", flush=True)
    # Perform a left outer merge the first dataframe with the sites file and save it back in the dictionary
    ## this ensures that we are only keeping sites in the sites file
    pop1=list(snpstat.keys())[0]
    snpstat[pop1]=pd.merge(snp_by_pop[['chr','pos']],snpstat[pop1],on=['chr', 'pos'], how='left')
    # Now perform left outer joins on all the files
    pop_list = []
    for pop in snpstat:
        pop_list.append(snpstat[pop])
    merged_pops = reduce((lambda left,right: pd.merge(left,right,on=['chr', 'pos'], how='left')), pop_list)
    return merged_pops

merged_pops=merge_pops(snpstat)
merged_pops

# Write outputs
merged_pops.to_csv(output_prefix+"_coverage", sep="\t", header=True, index=False)


