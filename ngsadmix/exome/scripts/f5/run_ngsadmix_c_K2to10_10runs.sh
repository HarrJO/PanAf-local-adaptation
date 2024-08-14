#!/bin/bash -l

# Job name
#$ -N ngsadmix_f5.0.5x.c_K1to10_repeat.x10

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/ngsadmix_output/mapped.on.target/f5.0.5x.subsp_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01/c

# Time requested
#$ -l h_rt=48:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=5G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10


#################################### JOB ####################################
date

INDIR="/home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f5.0.5x.subsp_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01_beagle/c"

# Load modules
module load htslib

# Run NGSadmix
for K in {2..10}
	do
	for RUN in {1..10}
		do
		/usr/bin/time --verbose \
		/home/ucfajos/bin/angsd/misc/NGSadmix \
		-likes ${INDIR}/f5.0.5x.c.beagle.gz \
		-K $K \
		-P 10 \
		-o f5.0.5x.c_K${K}_run${RUN} 
		done
	done

date

#################################### DONE ####################################

