#!/bin/bash -l

# Job name
#$ -N snpStat2coverage_chr21.py_chr21.f7.ce.pops_minInd.6or50pct_missing.pops.0.3_minMAC2

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_chr21_output/allele.frequencies/chr21.f7.ce.pops_minInd.6or50pct_missing.pops.0.3

# Time requested
#$ -l h_rt=48:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=200G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 1

#################################### JOB ####################################
date
# Load python
module load python3/3.8

# Run snpStat2coverage_chr21.py

/usr/bin/time --verbose python3 -u \
/home/ucfajos/analysis/phase1and2_chr21_analysis/coverage/scripts/snpStat2coverage_chr21.py \
-in_dir /home/ucfajos/Scratch/output/phase1and2_chr21_output/angsd_output/chr21.f7.pops.chrs_minInd.6or50pct_doMajorMinor.1_HWE.p.1e-3_beagle \
-out_prefix chr21.f7.ce.pops_chr21_missing.pops.0.3_pop_minMAC2 \
-snp_by_pop_file chr21.f7.ce.pops_chr21_missing.pops.0.3_pop.allele.counts_minMAC2

date
#################################### DONE ####################################
