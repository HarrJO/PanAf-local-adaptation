#!/bin/bash -l

# Job name
#$ -N pcangsd_f0.all_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/pcangsd_output/mapped.on.target/f0.all_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01

# Time requested
#$ -l h_rt=03:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=1G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10


#################################### JOB ####################################

INDIR='/home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f0.all_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01_beagle'

# Load modules
module load python3/3.7

/usr/bin/time --verbose \
python3 ~/bin/pcangsd/pcangsd.py \
-beagle ${INDIR}/f0.all.beagle.gz \
-o f0.all \
-threads 10

#################################### DONE ####################################
