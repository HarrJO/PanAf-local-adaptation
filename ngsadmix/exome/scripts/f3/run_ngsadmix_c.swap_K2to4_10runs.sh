#!/bin/bash -l

# Job name
#$ -N ngsadmix_f3.0.5x.c.swap_K2to4_repeat.x10

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/ngsadmix_output/mapped.on.target/f3.0.5x.c.swap_minInd.5_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.05

# Time requested
#$ -l h_rt=12:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=5G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10


#################################### JOB ####################################
date

INDIR="/home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f3.0.5x.c.swap_minInd.5_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.05_beagle"

# Load modules
module load htslib

# Run NGSadmix
for K in {2..4}
	do
	for RUN in {1..10}
		do
		/usr/bin/time --verbose \
		/home/ucfajos/bin/angsd/misc/NGSadmix \
		-likes ${INDIR}/f3.0.5x.c.swap.beagle.gz \
		-K $K \
		-P 10 \
		-o f3.0.5x.c.swap_K${K}_run${RUN} 
		done
	done

date

#################################### DONE ####################################

