This repository contains data and code from Ostridge et al. 

Preprint: https://doi.org/10.1101/2024.07.09.601734

This code starts with the mapped on target BAM files for the exome and chr21 (Fontsere et al., 2022) from PanAf project 
(http://panafrican.eva.mpg.de/). The scripts here then perform sample filtering, population structure analyses, population allele 
frequency estimation and tests for signatures of local adaptation. This directory contains scripts that were run on an HPC and 
others that can be run locally.

Each of the subdirectories listed below have their own README file going into greater detail.

genomic_data
- This directory contains genomic data.
- All exome sequence data (BAMs and FASTQs) will be available on ENA under the accession code ENA:PRJEB76176.
- All chr21 sequence data (BAMs and FASTQs) is available on ENA under the accession code ENA:PRJEB46115.

angsd
- This directory contains all the job scripts used to run ANGSD (Korneliussen, Albrechtsen & Nielsen 2014) on an HPC. ANGSD is 
a method that estimates genotype likelihoods (GLs) and allele frequencies. Estimated GLs can be used by other programmes in the ANGSD 
family (e.g., PCAngsd, NGSadmix and NGSrelate) to preform various population genomic analyses.
- ANGSD is run at each filtering stage to generate files required for the next stage.
- After filtering, ANGSD is run again for population structure analyses and to estimate population allele frequencies for the 
BayPass analysis.

pcangsd
- PCAngsd (Meisner and Albrechtsen 2018) is used to perform PCA using GLs estimated in ANGSD. 
- This script contains bash job scripts used to run this on a HPC.
- This is used in sample filtering to identify sample swaps and outliers indicative of contamination (./sample_filtering) and 
for population structure analyses (./population_structure).
- This also contains experimental scripts running PCAngsd's selection method, however, our discrete sampling strategy (which is not ideal for this method) and poor correspondence with BayPass core (genetics-only) results means we have little faith in these results. 

ngsadmix
- NGSadmix (Skotte Korneliussen and Albrechtsen 2013) is used to estimate individual admixture proportions using GLs estimated 
in ANGSD. 
- This script contains bash job scripts used to run this on a HPC.
- This is used in sample filtering (to identify sample swaps and outliers indicative of contamination) and for population 
structure analyses.

ngselate
- NGSrelate (Hanghøj et al. 2019) is used to estimate kinship coefficients using GLs estimated in ANGSD. 
- This script contains bash job scripts used to run this on a HPC.
- This is used in sample filtering to remove related samples.

sample_filtering
- This directory contains Rmarkdowns that plot and interrogate the results from PCAngsd/NGSadmix/NGSrelate and filter out 
samples with exceptionally low coverage, evidence of contamination or related samples.
- Exome data filtering is performed first, chr21 samples are then filtered for only those that passed the exome filtering and 
do not show any evidence of contamination in the chr21 data.

population_structure
- This directory contains Rmarkdowns that plot and interrogate the results from PCAngsd/NGSadmix for the filtered samples to 
investigate population structure and to define 'populations' by combining closely related sample sites.

allele_frequencies
- ANGSD is run population-by-population (and chromosomes-by-chromosome) to estimate population minor allele frequencies (job 
scripts in ./angsd); this directory contains custom python scripts and Rmarkdowns that collate all these files, reformat them and 
perform SNP filtering.
- The python scripts were run on a HPC and run chromosomes-by-chromosome (for the exome) so the process could be run in 
parallel. The large ANGSD output files means this process is very slow.

coverage
- This directory contains the scripts used to calculate the total coverage at each SNP from ANGSD output files, this is used 
later (in ./baypass) to perform post hoc tests for coverage differences across SNPs.

baypass
- BayPass (Gautier 2015) can be used to for evidence of local adaptation in the form of exceptionally high allele frequency 
differentiation (under the core model) or allele frequencies correlating with environmental covariables (under the AUX model).
- ./baypass/running_baypass contains the job scripts used to run BayPass on a HPC.
- ./baypass/analysing_baypass_output contains Rmarkdowns used to analyses the results of BayPass runs.

gowinda
- Gowinda is used to perform gene set enrichment analysis on the candidate SNPs identified in BayPass.
- ./gowinda/baypass_core contains the scripts required to generate input files, run gowinda and analyses the output. 
./gowinda/baypass_aux contains the scripts required to run gowinda and analyses the output - this directory used input files generates 
in ./gowinda/baypass_core so some of these scripts must be run first

pgls
- Phylogenetic Generalized Least Squares (PGLS) was run to confirm the BayPass GEA candidates
- ./pgls/scripts contains scripts required to run and plot the output of PGLS


