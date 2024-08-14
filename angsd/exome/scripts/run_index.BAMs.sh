#!/bin/bash -l

# Job name
#$ -N index_bams_mapped.on.target

# Set working directory
#$ -wd /home/ucfajos/Scratch/data/phase1and2_exomes/BAM.mapped.on.target

# Time requested
#$ -l h_rt=02:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=5G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 1

#################################### INFO ####################################
#################################### JOB ####################################

module load samtools
for BAM in /home/ucfajos/Scratch/data/phase1and2_exomes/BAM.mapped.on.target/*.bam; do samtools index ${BAM} 
done

#################################### DONE ####################################
