- This directory contains the scripts used to calculate the total coverage at each SNP from ANGSD output files, this is used 
later (in ./baypass) to perform post hoc tests for coverage differences across SNPs.
	
./$DATA/scripts/snpStat2coverage*.py
- This script takes the snpStat.gz files from running ANGSD for each population to estimate allele frequencies (in 
../angsd/$DATA/output) which reports the number of reads corresponding to each allele at each site in the population and sums it 
across all alleles and populations to get the total coverage depth at that site in the subspecies-dataset.
- The different scripts for exome and chr21 data are just due to inconsistent naming across the two, they are the same script 
in principle.
- ../allele_frequencies/scripts/maf2dac.v3.2.py must be run first to selct the SNPs of interest
	
./$DATA/scripts/ *.sh
Job scripts for running ./$DATA/scripts/snpStat2coverage*.py on the HPC (this isn't very computationally expensive).

./exome/scripts/samtools_depth/run_samtools.depth_number.of.bases.covered.sh
- Bash script that uses samtools depth to calculate the total number of sites covered by at least one read
- This is used in the f3 sample filtering script (../sample_filtering/exome/scripts/f3_exome.sample.filtering.Rmd).
- Outputs found in /exome/output