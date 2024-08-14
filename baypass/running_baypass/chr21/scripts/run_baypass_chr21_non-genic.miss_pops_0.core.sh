#!/bin/bash -l

# Job name
#$ -N baypass.{POPS}_minInd.6or50pct_missing.pops.0_minMAC2.non-genic_{FLANKS}.flanks_core{SEEDNAME}

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_chr21_output/baypass_output/ac/{POPS}_minInd.6or50pct_missing.pops.0_minMAC2.non-genic_{FLANKS}.flanks/core

# Time requested
#$ -l h_rt=12:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=10G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 8


#################################### JOB ####################################

INDIR="/home/ucfajos/Scratch/output/phase1and2_chr21_output/allele.frequencies/{POPS}_minInd.6or50pct_missing.pops.0"

# Run BayPass
/usr/bin/time --verbose \
~/bin/baypass_2.2/sources/i_baypass \
-nthreads 8 \
-gfile ${INDIR}/{POPS}_chr21_missing.pops.0.0_baypass.input.geno_minMAC2.non-genic_{FLANKS}.flanks \
-outprefix {POPS}_missing.pops.0_minMAC2.non-genic_{FLANKS}.flanks{SEEDNAME} \
-seed {SEED}

#################################### DONE ####################################

