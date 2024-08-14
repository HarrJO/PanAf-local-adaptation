#!/bin/bash -l

# run_minMAC.py.sh
# 15/06/2021

# Load python
module load python3/3.8

for minMAC in 2
	do
	# All
	## Small dataset
	python3 -u /home/ucfajos/analysis/phase1and2_exome_analysis/allele.frequencies/scripts/mapped.on.target/minMAC.py \
	-in_prefix /home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.pops_minInd.6or50pct_missing.pops.0/f5.pops.all.chrs_missing.pops.0 \
	-minMAC $minMAC

	## Large dataset
	python3 -u /home/ucfajos/analysis/phase1and2_exome_analysis/allele.frequencies/scripts/mapped.on.target/minMAC.py \
	-in_prefix /home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.pops_minInd.6or50pct_missing.pops.0.3/f5.pops.all.chrs_missing.pops.0.3 \
	-minMAC $minMAC
	done
