#!/bin/bash -l

# Job name
#$ -N prep.anc

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

cat ~/Scratch/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa >\
~/Scratch/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa

# bwa index
bwa index -a bwtsw ~/Scratch/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa

# samtools fadix
samtools faidx ~/Scratch/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa

#################################### DONE ####################################
