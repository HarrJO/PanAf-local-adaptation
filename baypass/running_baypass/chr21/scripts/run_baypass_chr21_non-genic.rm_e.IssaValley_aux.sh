#!/bin/bash -l

# Job name
#$ -N baypass_{POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_{FLANKS}bp.flanks.{COV}{SEEDNAME}

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_chr21_output/baypass_output/ac/{POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_{FLANKS}bp.flanks.rm_e.IssaValley/aux/{COV}

# Time requested
#$ -l h_rt=24:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=15G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 8


#################################### JOB ####################################

INDIR="/home/ucfajos/Scratch/output/phase1and2_chr21_output/allele.frequencies/{POPS}_minInd.6or50pct_missing.pops.0.3"
OMEGA="../../../{POPS}_minInd.6or50pct_missing.pops.0_minMAC2.non-genic_{FLANKS}bp.flanks.rm_e.IssaValley/core/{POPS}_missing.pops.0_minMAC2.non-genic_{FLANKS}bp.flanks.rm_e.IssaValley{SEEDNAME}_mat_omega.out"
EFILE="/home/ucfajos/analysis/phase1and2_exome_analysis/baypass/baypass_env_input/{COV}.{SUBSP}.baypass_env_input.txt"

# Run BayPass
/usr/bin/time --verbose \
~/bin/baypass_2.2/sources/i_baypass \
-nthreads 8 \
-gfile ${INDIR}/{POPS}_chr21_missing.pops.0.3_baypass.input.geno_minMAC2.non-genic_{FLANKS}bp.flanks.rm_e.IssaValley \
-omegafile ${OMEGA} \
-efile ${EFILE} \
-scalecov \
-auxmodel \
-outprefix {POPS}_missing.pops.0.3_minMAC2.non-genic_{FLANKS}bp.flanks.{COV}{SEEDNAME} \
-seed {SEED}

#################################### DONE ####################################
