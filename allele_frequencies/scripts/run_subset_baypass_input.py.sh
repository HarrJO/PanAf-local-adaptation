#!/bin/bash -l

# run_subset_baypass_input.py.sh
# 01/02/2024

# Just copy and paste the below into the commandline

# Load python
module load python3/3.8

for POPS in f5.pops
	do
	# Run script
	python3 -u /home/ucfajos/analysis/phase1and2_exome_analysis/allele.frequencies/scripts/mapped.on.target/subset_baypass_input.py \
	-baypass_input /home/ucfajos/Scratch/output_2024/phase1and2_exome_output/allele.frequencies/mapped.on.target/${POPS}_minInd.6or50pct_missing.pops.0.3/${POPS}.all.chrs_missing.pops.0.3_baypass.input.geno_minMAC2 \
	-allele_counts /home/ucfajos/Scratch/output_2024/phase1and2_exome_output/allele.frequencies/mapped.on.target/${POPS}_minInd.6or50pct_missing.pops.0.3/${POPS}.all.chrs_missing.pops.0.3_pop.allele.counts_minMAC2 \
	-n_subsets 3
	done
