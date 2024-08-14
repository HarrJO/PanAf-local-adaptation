#!/bin/bash -l

# Job name
#$ -N baypass_f5.w.pops_minInd.6or50pct_missing.pops.0_minMAC2_core_seed.10k

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/baypass_output/mapped.on.target/ac/f5.w.pops_minInd.6or50pct_missing.pops.0_minMAC2/core

# Time requested
#$ -l h_rt=12:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=10G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 8


#################################### JOB ####################################

INDIR="/home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.w.pops_minInd.6or50pct_missing.pops.0"

# Run BayPass
/usr/bin/time --verbose \
~/bin/baypass_2.2/sources/i_baypass \
-nthreads 8 \
-gfile ${INDIR}/f5.w.pops.all.chrs_missing.pops.0_baypass.input.geno_minMAC2 \
-outprefix f5.w.pops_missing.pops.0_minMAC2_seed.10k \
-seed 10000

#################################### DONE ####################################
