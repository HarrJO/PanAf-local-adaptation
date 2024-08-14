#!/usr/bin/env python
# coding: utf-8

# maf2dac.v3

# Harrison Ostridge, 04/05/2021

# Prior to this script, ANGSD should be run chromosome by chromosome and population by population to estimate population allele frequencies.
# ANSGD must be run chromosome by chromosome as it is too computationally demanding to run this script on all chromosomes at once.
# ANGSD must also be run population by population as this is the only way we can get allele frequency estimations using genotype likelihoods for each population.
# The purpose of this script is to pull together all population allele frequency files to form a file which describe allele frequencies across all populations (for each chromosome).
# The script also filters sites to include only bialleleic sites with data for enough populations.
# This script outputs files in the format required for BayPass (to run selection tests) and ANNOVAR (to annotate SNPs).
# The script is run chromosome by chromosome and so the outputs need to be combined afterwards.

# Import modules
import numpy as np
import pandas as pd
import os
from functools import reduce
import gc
import matplotlib.pyplot as plt
import copy
import argparse

# Argparse
parser = argparse.ArgumentParser(description="Combine ANGSD .maf.gz files from multiple populations into a allele counts file for input into BayPass.")
parser.add_argument("-in_dir", metavar="Input directory", type=str,
                    help="Input directory containing .maf.gz files.", default='./')
parser.add_argument("-out_dir", metavar="Ouptput directory", type=str, default='./',
                    help="If subspecies are specified then there must be subdirectories within the output directory for each subspecies (i.e. subdirectories called c, e, n and w).")
parser.add_argument("-out_prefix", metavar="Output file prefix", type=str, default='maf2dac.out',
                    help="The chromosome, subspecies and missing populations threshold will be added as suffixes to the file name.")
parser.add_argument("-subsp", metavar="Subspecies", nargs='+', default='all',
		    help="Either 'all' or list the first letters of the subspecies you want to include.")
parser.add_argument("-chr", metavar="Chromosome", type=str, choices=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],
                    help="Only autosomes.")
parser.add_argument("-snp_p", metavar="p-value", type=float, help="SNP p-value thresold (estimated by ANGSD).", default=0.000001)
parser.add_argument("-prop_missing_pops", metavar="Proportion of populations allowed to have missing data", type=float,
                    help="Specify the proportion of populations allowed to have missing data. If you require a site to have data in 70%% of populations then you should but 0.3 here.")
args = parser.parse_args()

# Modify/rename parameters
input_dir=args.in_dir
output_dir=args.out_dir
output_prefix=args.out_prefix
subsp=tuple(args.subsp)
chromo=args.chr
SNP_p_threshold=args.snp_p
proportion_missing_pops=args.prop_missing_pops
pops_to_exclude=["c.LaBelgique", "n.Mbe", "w.OutambaKilimi"]

## If a subspecies is specified
if "".join(subsp) != 'all':
    ### Add subspecies to file name
    output_name=output_prefix+"."+"".join(subsp)+".pops_chr"+str(chromo)+"_missing.pops."+str(proportion_missing_pops)
## If no subspecies is specified
else:
    output_name=output_prefix+".pops_chr"+str(chromo)+"_missing.pops."+str(proportion_missing_pops)

# Print parameters
print("-------------Running maf2dac.v3.py-------------")
print("Subspecies: "+"".join(subsp))
print("Chromosome: "+chromo)
print("Input directory: "+input_dir)
print("Output directory: "+output_dir)
print("Output name: "+output_name)
print("SNP p-value threshold: "+str(SNP_p_threshold))
print("Proportion of missing populations: "+str(proportion_missing_pops))
print("Populations excluded: "+str(pops_to_exclude))


# Define functions
## Read in data
def read_data():
    """Read in ANGSD .maf.gz output from a predefined directory
    This function looks in all subdirectories and reads files in alphabetical order"""
    print("Reading in data", flush=True)
    # If there is a subspecies specified
    if "".join(subsp) != 'all':
        filelist = sorted([os.path.join(root, name)
                for root, dirs, files in os.walk(input_dir)
                for name in files
                if name.startswith(subsp) 
                           & name.endswith("chr" + chromo + ".mafs.gz")])
    # If there is no subspecies specified
    else:
        filelist = sorted([os.path.join(root, name)
                for root, dirs, files in os.walk(input_dir)
                for name in files
                if name.endswith("chr" + chromo + ".mafs.gz")])
    # Exclude populations
    filelist_filtered = list()
    for file in filelist:
        if any(pop in file for pop in pops_to_exclude)==False:
            filelist_filtered.append(file)

    mafs = {}
    for file in filelist_filtered:
        file_name = os.path.splitext(os.path.splitext(os.path.basename(file))[0])[0].replace(".", "_")
        mafs[file_name] = np.array(pd.read_csv(file, sep="\t"))
        print(file_name, "read in", flush=True)
    return mafs

## Plot SFS from allele frequencies
def plot_SFS(af_dict, af_column, title=None, file_prefix=None):
    """Plot SFS from dictionary of np arrays (af_dict) with allele frequencies for each population"""
    fig, ax = plt.subplots(1,2, figsize=(15,5))
    for pop in af_dict:
        ax[0].hist(list(af_dict[pop][:,af_column]), 100, histtype='step', label=pop)
        ax[1].hist(list(af_dict[pop][:,af_column]), 100, histtype='step', label=pop)
    ax[0].set_title("Full SFS")
    ax[0].legend(loc='upper center',ncol=5,fontsize=3)
    ax[1].set_title("Zoomed in SFS")
    axes = plt.gca()
    ax[1].axes.set_ylim([0,500])
    ax[1].legend(loc='upper center',ncol=5,fontsize=3)
    fig.suptitle(title, fontsize=16)
    if file_prefix is not None:
        plt.savefig(output_dir+file_prefix+output_name+'.pdf')
    plt.show()

## Polarise alleles
def polarise_alleles(mafs):
    """Use ancestral allele reported in .maf.gz output to polarise the allele frequencies resulting in DAFs"""
    print("Polariasing alleles", flush=True)
    dafs = mafs
    for pop in mafs:
        for row in range (0, len(mafs[pop])):
            # If the anc allele = the minor allele
            if mafs[pop][row, 5] == mafs[pop][row, 3]:
                # Make the der allele (previously the major) the minor
                dafs[pop][row, 3] = mafs[pop][row, 2]
                # Make the anc allele the major
                dafs[pop][row, 2] = mafs[pop][row, 5]
                # Mafe the DAF from 1 - MAF
                dafs[pop][row, 6] = 1 - mafs[pop][row, 6]
            # If the anc allele does not equal the major or minor allele
            if mafs[pop][row, 5] != mafs[pop][row, 3] and mafs[pop][row, 5] != mafs[pop][row, 2]:
                dafs[pop][row, 6] = "NA"
        # Remove sites where neither allele equals the ancestral allele
        dafs[pop] = dafs[pop][dafs[pop][:,6] != "NA"]
        # Remove uneeded columns
        dafs[pop] = np.concatenate((dafs[pop][:,0:4], dafs[pop][:,6:9]), axis=1)
    return dafs 

## Identify locally monomorphic sites
def SNP_thresholds(daf_input):
    """Report sites as fixed in a population if they fail to pass a predefined minimum MAF and SNP p value threshold"""
    print("Identifying locally monomorphic sites", flush=True)
    dafs_SNP_thresholds = {}
    for pop in daf_input:
        # Make empty array for modified derived allele frequencies (DAFs)
        dafs = np.empty(len(daf_input[pop]))
        dafs[:] = np.nan
        # Make empty array for modified derived alleles
        das = np.empty(len(daf_input[pop]), dtype=str)
        das[:] = np.nan
        # For each position (row)
        for pos in range (0, len(daf_input[pop])):
            # If the SNP p is above the threshold -> round DAF to 0 or 1
            if daf_input[pop][pos, 5] > SNP_p_threshold:
                dafs[pos] = round(daf_input[pop][pos, 4])
            # If the SNP p is belove the threshold -> just keep DAF 
            if daf_input[pop][pos, 5] < SNP_p_threshold:
                dafs[pos] = daf_input[pop][pos, 4]
            # Min MAF filter
            ## min MAF set to 1/2n (singleton)
            min_MAF = 1/(2*daf_input[pop][pos, 6])
            ## If the SNP DAF is below the MAF threshold OR DAF is above 1-minimum MAF -> round DAF to 0 or 1
            if daf_input[pop][pos, 4] < min_MAF or daf_input[pop][pos, 4] > (1-min_MAF):
                dafs[pos] = round(daf_input[pop][pos, 4])
            # If the rounded DAF is now 0, let the derived allele = Na (because there is no derrived allele)
            if dafs[pos] == 0:
                das[pos] = 'n'
            # If the rounded DAF is not 0, keep the derived allele as it is
            if dafs[pos] != 0:
                das[pos] = daf_input[pop][pos, 3]
        # Reshapse vectors so they can be merged
        dafs = dafs.reshape(-1, 1)
        das = das.reshape(-1, 1)
        # Return table with position coordinates and the ancestral and derived allele counts
        dafs_SNP_thresholds[pop]=np.concatenate((daf_input[pop][:,0:3], das, dafs, daf_input[pop][:,5:7]), axis=1)
    return dafs_SNP_thresholds

## Allele frequencies to allele counts
def daf2dac(dafs):
    """Convert derived allele requencies (DAFs) into derived allele counts (DACs) 
    by multiplying DAF by the haploid sample size"""
    print("Converting allele frequencies to counts", flush=True)
    dacs = {}
    for pop in dafs:
        # calculate ancestral allele frequency
        aafs = 1 - dafs[pop][:,4]
        # derived and ancestral counts (*2 because they're diploid)
        dac = (2*dafs[pop][:,4]*dafs[pop][:,6]).reshape(-1, 1)
        aac = (2*aafs*dafs[pop][:,6]).reshape(-1, 1)
        # return table with position coordinates and the ancestral and derived allele counts
        ## Ensure DAC comes before AAC otherwise BayPass outputs are back to front
        dacs[pop] = pd.DataFrame(np.concatenate((dafs[pop][:,0:4], dac, aac), axis=1), columns=['chr', 'pos', "aa", pop+".da", pop+".dac", pop+".aac"])
    return dacs

## Merge communities
def merge_pops(dacs):
    """Combine derived allele count arrays from multiple populations into one large dataframe"""
    print("Merging populations", flush=True)
    pop_list = []
    for pop in dacs:
        pop_list.append(dacs[pop])
    merged_pops = reduce((lambda left,right: pd.merge(left,right,on=['chr', 'pos', 'aa'], how='outer')), pop_list)
    return merged_pops

## Plot SFS from allele counts
def plot_SFS_from_counts_per_pop(merged_pops, plot_title='', file_prefix=None, ylim=None):
    """Plot the SFS from dataframe containing derived and ancestral allele counts"""
    fig, ax = plt.subplots(1,3, figsize=(20,5))
    # Calculate DAFs per population
    dacs = merged_pops.loc[:, merged_pops.columns.str.endswith("dac")].fillna(0)
    aacs = merged_pops.loc[:, merged_pops.columns.str.endswith("aac")].fillna(0)
    for pop in range(0, len(dacs.columns)):
        pop_dafs = dacs.iloc[:, pop]/(dacs.iloc[:, pop]+aacs.iloc[:, pop])
        ax[0].hist(list(pop_dafs), 100, histtype='step', label=dacs.columns[pop])
    # Calculate total DAF
    global_dac = (merged_pops.loc[:, merged_pops.columns.str.endswith("dac")].fillna(0)).sum(axis=1)
    global_total = (merged_pops.loc[:, merged_pops.columns.str.endswith("ac")].fillna(0)).sum(axis=1)
    global_daf = global_dac/global_total
    ax[1].hist(list(global_daf), 100, histtype='step', label="Total")
    ax[2].hist(list(global_daf), 100, histtype='step', label="Total")
    # Plot parameters
    #axes = fig.gca()
    if ylim is not None:
        ax[0].set_ylim([0, ylim])
        ax[1].set_ylim([0, ylim])
    ax[0].legend(ncol=4, fontsize=5)
    ax[0].set_title("SFS per population, zoomed in")
    ax[1].legend(ncol=1, fontsize=5)
    ax[1].set_title("Total SFS, zoomed in")
    ax[2].legend(ncol=1, fontsize=5)
    ax[2].set_title("Total SFS")
    fig.suptitle(plot_title, fontsize=16)
    if file_prefix is not None:
        plt.savefig(output_dir+file_prefix+output_name+'.pdf')
    plt.show()

## Remove globally monomorphic sites
def remove_monomorphic(merged_pops):
    """Remove sites which are globally monomorphic i.e. monomorphic for the same allele in every population"""
    print("Removing globally monomorphic sites", flush=True)
    # Initialise polymrphic site mask
    polymorphic_mask = np.empty(len(merged_pops), dtype='bool')
    polymorphic_mask[:] = True
    # Create extract derived and ancestral allele counts
    dacs = np.array(merged_pops.loc[:, merged_pops.columns.str.endswith("dac")])
    aacs = np.array(merged_pops.loc[:, merged_pops.columns.str.endswith("aac")])
    # For each row (i.e. genomic site/position)
    for pos in range (0, len(merged_pops)):
        ## If the sum of derived allele frequencies or the sum of ancestral allele frequencies == 0 -> site is not polymorphic
        if np.nansum(dacs[pos,:]) == 0 or np.nansum(aacs[pos,:]) == 0:
            polymorphic_mask[pos] = False
    # Apply mask to data, only keep polymorphic sites
    return merged_pops[polymorphic_mask]

## Missing populations filter
def missing_pops_filter(merged_pops, proportion_missing_pops):   
    """Remove sites with missing data in X% of all populations"""
    print("Filtering by number of missing populations", flush=True)
    # Define filter
    missing_population_filter = merged_pops.loc[:, merged_pops.columns.str.endswith("dac")].isnull().sum(axis=1) <= (sum(merged_pops.columns.str.endswith("da"))*proportion_missing_pops)
    # Apply filter
    merged_pops=merged_pops[missing_population_filter]
    return merged_pops

## Remove triallelic sites
def biallelic_filter(merged_pops):
    """Remove sites with more than 1 derived allele 
    (sites have already been filtered to have the same ancestral allele so this filters for bialleleic sites)"""
    print("Removing triallelic stites", flush=True)
    # Extract derived alleles column for each population
    merged_pops_das = np.array(merged_pops.loc[:, merged_pops.columns.str.endswith("da")])
    # Initialise biallelic filter 
    bases = ['A', 'T', 'G', 'C']
    biallelic = np.empty(len(merged_pops_das), dtype=bool)
    biallelic[:] = False
    # For each site
    for pos in range (0, len(merged_pops_das)):
        # If there is only one derrived allele, return true to bialalic mask
        if len([i for i in bases if i in merged_pops_das[pos, :]]) == 1:
            biallelic[pos] = True
    # Apply mask
    merged_pops = merged_pops[biallelic]
    return merged_pops

## Write outputs
def write_allele_counts(merged_pops, filter_table, output_dir, output_name):
    """Write all outputs"""
    print("Writing outputs", flush=True)
    # Reset indexes and make position column an integer (not float)
    merged_pops.reset_index(drop=True, inplace=True)
    merged_pops['pos'] = merged_pops['pos'].astype(int)
    # Write allele counts file: very similar to the BayPass input file but with the populations (column names) and SNP locations clearly indicated
    ## Select allele count columns
    allele_counts = merged_pops.loc[:, merged_pops.columns.str.endswith("ac")].fillna(0)
    ## Select site coordinate columns
    SNPs = merged_pops[['chr', 'pos']]
    ## Conbine site position and allele counts columns and write
    pd.concat([SNPs, allele_counts], axis=1).to_csv(os.path.join(output_dir, output_name+"_pop.allele.counts"), sep="\t", header=True, index=False)
    # BayPass
    ## Population names: BayPass output does not give the population names but the columns will be in this order
    ### Select derived alelel column names 
    pops = pd.DataFrame(merged_pops.columns[merged_pops.columns.str.endswith("da")]).fillna(0)
    ### Remove '.da' suffix
    pops[0] = pops[0].str.split('_chr').str[0]
    ### Write 
    pops.to_csv(os.path.join(output_dir, output_name+"_pop.names"), sep="\t", header=False, index=False)
    ## BayPass input file
    ### Select only allele counts columns
    merged_pops_BayPassInput_all = merged_pops.loc[:, merged_pops.columns.str.endswith("ac")].fillna(0)
    ### Write
    merged_pops_BayPassInput_all.to_csv(os.path.join(output_dir, output_name+"_baypass.input.geno"), sep=" ", header=False, index=False)
    # Filter table
    filter_table.to_csv(os.path.join(output_dir, output_name+"_filter.table.csv"), sep=",", header=True, index=True)

## Write ANNOVAR input
def write_ANNOVAR_input(merged_pops, output_dir, output_name):
    """Write output in format for ANNOVAR input with the following columns:
    chr     pos     pos     aa      da      daf"""
    # Get global derived allele count by summming all derived allele columns
    global_dac = (merged_pops.loc[:, merged_pops.columns.str.endswith("dac")].fillna(0)).sum(axis=1)
    # Get total global allele count (ancestral and derived)
    global_total = (merged_pops.loc[:, merged_pops.columns.str.endswith("ac")].fillna(0)).sum(axis=1)
    # Get DAF by dividing derived allele count by total allele count
    global_daf = global_dac/global_total
    # Make DAf into a pd colum with an apropriate header
    global_daf = pd.DataFrame(global_daf, columns=['daf'])
    # Select all the population derrrived allele columns
    com_das = merged_pops.loc[:, merged_pops.columns.str.endswith("da")]
    # Create empty np.array for derrived allele column
    das = np.empty(len(merged_pops), dtype=str)
    das[:] = np.nan
    # Select the base which is in each row (not n or NaN)
    bases = ['A', 'T', 'G', 'C']
    for pos in range (0, len(merged_pops)):
        das[pos] = [i for i in bases if i in np.array(com_das.iloc[[pos]])][0]
    # Give an appropriate header
    das = pd.DataFrame(das, columns=['da'])
    # Combine the columns together and print file
    snps = merged_pops[merged_pops.columns[0:3]]
    annovar_input = pd.concat([merged_pops[merged_pops.columns[0:2]], merged_pops[merged_pops.columns[1]],
                               merged_pops[merged_pops.columns[2]], das, global_daf], axis=1)
    annovar_input.to_csv(os.path.join(output_dir, output_name+"_annovar.input"), sep="\t", header=True, index=False)

# Main
def main():
    # Read in data
    mafs=read_data()
    # Plot SFS from allele frequencies
    plot_SFS(mafs, af_column=6, title='Input Folded SFS', file_prefix='SFSinput_')
    # Polarise alleles
    dafs=polarise_alleles(mafs)
    # Identify locally monomorphic sites
    dafs_SNP_thresholds=SNP_thresholds(dafs)
    # Allele frequencies to allele counts
    dacs=daf2dac(dafs_SNP_thresholds)
    # Merge communities
    merged_pops=merge_pops(dacs)
    ## Start making filter table
    total_sites = len(merged_pops)
    filter_table = pd.DataFrame(columns = ['Filter', 'Number of sites', 'Percentage of sites'])
    filter_table = filter_table.append({'Filter':'Total', 'Number of sites': len(merged_pops), 'Percentage of sites':((len(merged_pops)*100)/total_sites)}, ignore_index=True)
    # Remove globally monomorphic sites
    merged_pops = remove_monomorphic(merged_pops)
    ## Add to filter table
    filter_table=filter_table.append({'Filter':'Remove globally monomorphic', 'Number of sites': len(merged_pops), 'Percentage of sites':round((len(merged_pops)*100)/total_sites, 3)}, ignore_index=True)
    # Missing populations filter
    merged_pops=missing_pops_filter(merged_pops, proportion_missing_pops)
    ## Add to filter table
    filter_table=filter_table.append({'Filter':'Missing populations filter', 'Number of sites': len(merged_pops), 'Percentage of sites':round((len(merged_pops)*100)/total_sites, 3)}, ignore_index=True)
    # Remove triallelic sites
    merged_pops=biallelic_filter(merged_pops)
    ## Add to filter table
    filter_table=filter_table.append({'Filter':'Remove globally triallelic sites', 'Number of sites': len(merged_pops), 'Percentage of sites':round((len(merged_pops)*100)/total_sites, 3)}, ignore_index=True)
    # Write outputs
    write_allele_counts(merged_pops, filter_table, output_dir, output_name)
    # Write ANNOVAR input
    write_ANNOVAR_input(merged_pops, output_dir, output_name)
    # Plot SFS from allele counts
    plot_SFS_from_counts_per_pop(merged_pops, plot_title="Output unfolded SFS", ylim=300, file_prefix='SFSoutput_')
    
# Run main
if __name__ == "__main__":
    main()

