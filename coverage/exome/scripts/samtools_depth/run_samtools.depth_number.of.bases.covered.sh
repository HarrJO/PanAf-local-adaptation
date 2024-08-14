#!/bin/bash -l

# Job name
#$ -N samtools.depth_number.of.bases.covered

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/coverage/mapped.on.target/number.of.bases.covered

# Time requested
#$ -l h_rt=12:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=1G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10

#################################### INFO ####################################
# This script is for calculating the number of sites with at least 1 read.
# I wanted to see if troublesome samples could be explained because they may have a high average coverage but lots of missing sites
# This is done by running samtools depth which calculates read depth at all sites with any coverage.
## By counting the number of rows in the file you have the number of bases covered.
#################################### JOB ####################################

# Load modules
module load samtools/1.9/gnu-4.9.2

# Parameters
OUTPUT=/home/ucfajos/Scratch/output/phase1and2_exome_output/coverage/mapped.on.target/number.of.bases.covered/number.of.bases.covered.per.sample

# Add flanks to exome coordinate BED
## Remove header
#sed '1,2d' /home/ucfajos/Scratch/data/exome_coordinates/S07604514/S07604514_Regions.bed \
#>/home/ucfajos/Scratch/output/phase1and2_exome_output/target.regions_output/bed.files/S07604514_Regions_no.head.bed

# cd into input directory so when looping over file names the directory structure isn't included
cd /home/ucfajos/Scratch/data/phase1and2_exomes/BAM.mapped.on.target/

# Make filelist so I know the order of the columns
FILELIST=`ls *.bam`

# Headers for output
echo sample number_of_bases > $OUTPUT

# For each BAM
for FILE in `ls *.bam`;
        do
        # Prepare output file name
        SAMPLE="${FILE/_exome.addRG.bam/}"
        
	# Samtoolsdepth
	samtools depth \
	/home/ucfajos/Scratch/data/phase1and2_exomes/BAM.mapped.on.target/${FILE} \
	>/home/ucfajos/Scratch/output/phase1and2_exome_output/coverage/mapped.on.target/number.of.bases.covered/${SAMPLE}.depth

	# Calculate the number of rows (number of bases with reads
	NROW=`wc -l < /home/ucfajos/Scratch/output/phase1and2_exome_output/coverage/mapped.on.target/number.of.bases.covered/${SAMPLE}.depth`

	# Add number of bases to output file
	echo $SAMPLE $NROW >> $OUTPUT
	
	# rm depth file as they are quite large and no longer needed
	rm /home/ucfajos/Scratch/output/phase1and2_exome_output/coverage/mapped.on.target/number.of.bases.covered/${SAMPLE}.depth
	done

#################################### DONE ####################################
