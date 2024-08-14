#!/bin/bash -l

# Job name
#$ -N baypass_f5.ce.pops_minInd.6or50pct_missing.pops.0.3_minMAC2_aux.f_over_sum_known_trees.rm_e.IssaValley_seed.10k

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/baypass_output/mapped.on.target/ac/f5.ce.pops_minInd.6or50pct_missing.pops.0.3_minMAC2.rm_e.IssaValley/aux/f_over_sum_known_trees.rm_e.IssaValley

# Time requested
#$ -l h_rt=48:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=15G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 8


#################################### JOB ####################################


############ NB: COV will contain the suffix '.rm_e.IssaValley' ############


INDIR="/home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.ce.pops_minInd.6or50pct_missing.pops.0.3"
OMEGA="/home/ucfajos/Scratch/output/phase1and2_exome_output/baypass_output/mapped.on.target/ac/f5.ce.pops_minInd.6or50pct_missing.pops.0_minMAC2.rm_e.IssaValley/core/f5.ce.pops_missing.pops.0_minMAC2.rm_e.IssaValley_seed.10k_mat_omega.out"
EFILE="/home/ucfajos/analysis/phase1and2_exome_analysis/baypass/baypass_env_input/f_over_sum_known_trees.rm_e.IssaValley.ce.baypass_env_input.txt"

# Run BayPass
/usr/bin/time --verbose \
~/bin/baypass_2.2/sources/i_baypass \
-nthreads 8 \
-gfile ${INDIR}/f5.ce.pops.all.chrs_missing.pops.0.3_baypass.input.geno_minMAC2.rm_e.IssaValley \
-omegafile ${OMEGA} \
-efile ${EFILE} \
-scalecov \
-auxmodel \
-outprefix f5.ce.pops_missing.pops.0.3_minMAC2_f_over_sum_known_trees.rm_e.IssaValley_seed.10k \
-seed 10000

#################################### DONE ####################################
