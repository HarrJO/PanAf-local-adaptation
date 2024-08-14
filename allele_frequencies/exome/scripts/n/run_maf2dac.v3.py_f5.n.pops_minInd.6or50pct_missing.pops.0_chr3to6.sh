#!/bin/bash -l

# Job name
#$ -N maf2dac.v3.py_f5.n.pops_minInd.6or50pct_missing.pops.0_chr3to6

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.n.pops_minInd.6or50pct_missing.pops.0

# Time requested
#$ -l h_rt=12:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=100G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 1

#################################### JOB ####################################
date
# Load python
module load python3/3.8

# Run maf2dac.v3.py

for CHR in {3..6}
	do
	/usr/bin/time --verbose python3 -u \
	/home/ucfajos/analysis/phase1and2_exome_analysis/allele.frequencies/scripts/mapped.on.target/maf2dac.v3.py \
	-in_dir /home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f5.0.5x.pops.chrs_minInd.6or50pct_doMajorMinor.1_HWE.p.1e-3_beagle \
	-out_dir ./ \
	-out_prefix f5 \
	-subsp n \
	-chr $CHR \
	-snp_p 0.000001 \
	-prop_missing_pops 0

	echo "Finished chromosome" ${CHR}

	done

date
#################################### DONE ####################################
