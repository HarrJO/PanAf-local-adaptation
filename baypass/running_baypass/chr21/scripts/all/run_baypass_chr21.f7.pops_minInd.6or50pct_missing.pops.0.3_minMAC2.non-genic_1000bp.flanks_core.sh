#!/bin/bash -l

# Job name
#$ -N baypass.chr21.f7.pops_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_1000bp.flanks_chr21.omega_core_seed.10k

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_chr21_output/baypass_output/ac/chr21.f7.pops_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_1000bp.flanks/core

# Time requested
#$ -l h_rt=12:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=10G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 8


#################################### JOB ####################################

INDIR="/home/ucfajos/Scratch/output/phase1and2_chr21_output/allele.frequencies/chr21.f7.pops_minInd.6or50pct_missing.pops.0.3"
OMEGA="/home/ucfajos/Scratch/output/phase1and2_chr21_output/baypass_output/ac/chr21.f7.pops_minInd.6or50pct_missing.pops.0_minMAC2/core/chr21.f7.pops_missing.pops.0_minMAC2_mat_omega.out"

# Run BayPass
/usr/bin/time --verbose \
~/bin/baypass_2.2/sources/i_baypass \
-nthreads 8 \
-gfile ${INDIR}/chr21.f7.pops_chr21_missing.pops.0.3_baypass.input.geno_minMAC2.non-genic_1000bp.flanks \
-omegafile ${OMEGA} \
-outprefix chr21.f7.pops_missing.pops.0.3_minMAC2.non-genic_1000bp.flanks_chr21.omega_seed.10k \
-seed 10000

#################################### DONE ####################################

