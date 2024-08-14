#!/usr/bin/env python
# coding: utf-8

## minMAC

#### Harrison Ostridge, 15/06/2021

# This applied a minimum allele count (MAC) filter.
# This is to remove very low frequency alleles (singletons and doubletons mainly) to try and improve the fit of XtX distributions to chi squared. 

# Import modules
import numpy as np
import pandas as pd
import argparse

# Argparse
parser = argparse.ArgumentParser(description="Apply minimum allele count filter")
parser.add_argument("-in_prefix", metavar="Input prefix", type=str,
                    help="Full path and prefix of files output by maf2dac which end in _pop.allele.counts, _baypass.input.geno or _annovar.input.")
parser.add_argument("-minMAC", metavar="Minimum minor allele count", type=int,
                    help="e.g. 2 would require sites to have a count of 2 (removeing singletons), requiring 3 would remove doubletons etc.")
args = parser.parse_args()

# Rename (to make testing in jupyter notebook easier)
in_prefix=args.in_prefix
minMAC=args.minMAC

# Read in 
allele_counts = pd.read_csv(in_prefix+"_pop.allele.counts", sep="\t")
baypass_input = pd.read_csv(in_prefix+"_baypass.input.geno", sep=" ", header=None)
annovar_input = pd.read_csv(in_prefix+'_annovar.input', sep="\t")

# Derived allele count
dac = allele_counts.loc[:, allele_counts.columns.str.endswith("dac")].sum(axis=1)
# Ancestral allele count
aac = allele_counts.loc[:, allele_counts.columns.str.endswith("aac")].sum(axis=1)
# Reuire a site to have a DAC >= minMAC and AAC >= minMAC
minMACmask=(dac>=minMAC) & (aac>=minMAC)

# Filter and write file
allele_counts[minMACmask].to_csv(in_prefix+"_pop.allele.counts_minMAC"+str(minMAC), sep="\t", header=True, index=False)
baypass_input[minMACmask].to_csv(in_prefix+"_baypass.input.geno_minMAC"+str(minMAC), sep=" ", header=False, index=False)
annovar_input[minMACmask].to_csv(in_prefix+"_annovar.input_minMAC"+str(minMAC), sep="\t", header=True, index=False)
