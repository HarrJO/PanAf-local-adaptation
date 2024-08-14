#!/bin/bash -l

# Job name
#$ -N baypass_{POPS}_minInd.6or50pct_missing.pops.{MISS_POPS1}_minMAC2_aux.{COV}{SUBSET}{SEEDNAME}

# Set working directory
#$ -wd /home/ucfajos/Scratch/output_2024/phase1and2_exome_output/baypass_output/mapped.on.target/ac/{POPS}_minInd.6or50pct_missing.pops.{MISS_POPS1}_minMAC2/auxi/{COV}

# Time requested
#$ -l h_rt={HOURS}:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=15G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 8


#################################### JOB ####################################

INDIR="/home/ucfajos/Scratch/output_2024/phase1and2_exome_output/allele.frequencies/mapped.on.target/{POPS}_minInd.6or50pct_missing.pops.{MISS_POPS1}"
OMEGA="/home/ucfajos/Scratch/output_2024/phase1and2_exome_output/baypass_output/mapped.on.target/ac/{POPS}_minInd.6or50pct_missing.pops.0_minMAC2/core/{POPS}_missing.pops.0_minMAC2{SEEDNAME}_mat_omega.out"
EFILE="/home/ucfajos/analysis/phase1and2_exome_analysis/baypass/baypass_env_input/{COV}.{SUBSP}.baypass_env_input.txt"

# Run BayPass
/usr/bin/time --verbose \
~/bin/baypass_2.2/sources/i_baypass \
-nthreads 8 \
-gfile ${INDIR}/{POPS}.all.chrs_missing.pops.{MISS_POPS1}_baypass.input.geno_minMAC2{SUBSET} \
-omegafile ${OMEGA} \
-efile ${EFILE} \
-scalecov \
-auxmodel \
-outprefix {POPS}_missing.pops.{MISS_POPS1}_minMAC2_{COV}{SUBSET}{SEEDNAME} \
-seed {SEED}

cp /home/ucfajos/analysis/phase1and2_exome_analysis/baypass/baypass_env_input/{COV}.{SUBSP}.txt ./

#################################### DONE ####################################
