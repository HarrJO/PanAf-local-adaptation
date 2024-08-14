#!/bin/bash -l

# Job name
#$ -N pcangsd_f5.0.5x.subsp.from.all_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/pcangsd_output/mapped.on.target/f5.0.5x.subsp.from.all_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01

# Time requested
#$ -l h_rt=02:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=1G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10

################################### INFO ####################################
# Here I try an alternative method to produce PCAs per subspecies by using the beagle file estimated with all samples and then providing PCAngsd with a individual filter for each subspecies
# I am trying this because some samples do not cluster well on the between subspecies scale but look ok within subspecies 
#################################### JOB ####################################

INDIR='/home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f5.0.5x.all_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01_beagle'
FILT_DIR='/home/ucfajos/analysis/phase1and2_exome_analysis/pcangsd/ind.filters/mapped.on.target/f5'

# Load modules
module load python3/3.7

for SUBSP in c e n w
	do
	/usr/bin/time --verbose \
	python3 ~/bin/pcangsd/pcangsd.py \
	-beagle ${INDIR}/f5.0.5x.all.beagle.gz \
	-o ${SUBSP}/f5.0.5x.${SUBSP} \
	-threads 10 \
	-filter ${FILT_DIR}/${SUBSP}.inds
	done

#################################### DONE ####################################
