#!/bin/bash -l

# Job name
#$ -N maf2dac.v3.2.py_chr21.f7.w.pops_minInd.6or50pct_missing.pops.0.3

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_chr21_output/allele.frequencies/chr21.f7.w.pops_minInd.6or50pct_missing.pops.0.3

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

# Run maf2dac.v3.2.py

/usr/bin/time --verbose python3 -u \
/home/ucfajos/analysis/phase1and2_exome_analysis/allele.frequencies/scripts/mapped.on.target/maf2dac.v3.2.py \
-in_dir /home/ucfajos/Scratch/output/phase1and2_chr21_output/angsd_output/chr21.f7.pops.chrs_minInd.6or50pct_doMajorMinor.1_HWE.p.1e-3_beagle \
-out_dir ./ \
-out_prefix chr21.f7 \
-subsp w \
-chr 21 \
-snp_p 0.000001 \
-prop_missing_pops 0.3


date
#################################### DONE ####################################
