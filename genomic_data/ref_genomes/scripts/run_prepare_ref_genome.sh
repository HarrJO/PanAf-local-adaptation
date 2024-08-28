#!/bin/bash -l

# Job name
#$ -N prep.ref

# Time requested
#$ -l h_rt=5:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=1G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Set working directory
#$ -wd /home/ucfajos/Scratch/ref_genomes/

# Number of cores
#$ -pe smp 20


#################################### JOB ####################################

# INFO: method taken from https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk

# Load modules
module load bwa
module load samtools

# bwa index
bwa index -a bwtsw hg19.fa

# samtools fadix
samtools faidx hg19.fa
