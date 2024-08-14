#!/bin/bash -l

# Job name
#$ -N snpStat2coverage.py_f5.ce.pops_minInd.6or50pct_missing.pops.0.3_minMAC2

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.ce.pops_minInd.6or50pct_missing.pops.0.3

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

# Run snpStat2coverage.py

/usr/bin/time --verbose python3 -u \
/home/ucfajos/analysis/phase1and2_exome_analysis/coverage/scripts/snpStat2coverage.py \
-in_dir /home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f5.0.5x.pops.chrs_minInd.6or50pct_doMajorMinor.1_HWE.p.1e-3_beagle \
-out_prefix f5.ce.pops.all.chrs_missing.pops.0.3_minMAC2 \
-snp_by_pop_file f5.ce.pops.all.chrs_missing.pops.0.3_pop.allele.counts_minMAC2

date
#################################### DONE ####################################
