#!/bin/bash -l

# Job name
#$ -N angsd_f5.0.5x.all_doMajorMinor.1_sites.and.pops.to.match.baypass.input.minMAC2_beagle_2024

# Set working directory
#$ -wd /home/ucfajos/Scratch/output_2024/phase1and2_exome_output/angsd_output/mapped.on.target/f5.0.5x.all_doMajorMinor.1_sites.and.pops.to.match.baypass.input.minMAC2_beagle

# Time requested
#$ -l h_rt=48:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=5G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10

#################################### INFO ####################################
# This script is run with doMajorMinor 1 (to ensure no excess of high freq DAFs) and a HWE filter (to ensure no bump at 0.5 due to paralogs).
# GLF is in beagle format for PCAngsd
# All samples have been fully filtered
# This output if for running the PCAngsd selection detection method to compare to BayPass
## I therefore filter to only include sites in the baypass input/output
## NB: I still need read/mapping filters to remove bad reads but I don't need the site filters (minMAF ect)
#################################### JOB ####################################

# Load modules
module load htslib/1.7 

# Make sites file
## Remove header and select only the first two columns (chr and pos)
awk 'NR > 1 { print }' \
/home/ucfajos/Scratch/output_2024/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.pops_minInd.6or50pct_missing.pops.0.3/f5.pops.all.chrs_missing.pops.0.3_pop.allele.counts_minMAC2 | \
awk -v OFS='\t' '{print $1,$2}' > sites.file

## Index sites file
/home/ucfajos/bin/angsd/angsd sites index sites.file

# Run ANGSD
/usr/bin/time --verbose \
~/bin/angsd/angsd \
-b /home/ucfajos/analysis/phase1and2_exome_analysis/angsd/bam.filelists/mapped.on.target/f5_pops.inc.in.baypass/bam.filelist.all_2024 \
-out f5.0.5x.all \
-ref ~/Scratch/data_2024/ref_genomes/hg19.fa \
-anc ~/Scratch/data_2024/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa \
-uniqueOnly 1 \
-remove_bads 1 \
-only_proper_pairs 1 \
-trim 0 \
-C 50 \
-baq 1 \
-GL 2 \
-doGlf 2 \
-minMapQ 30 \
-nThreads 10 \
-doMajorMinor 1 \
-doMaf 2 \
-doHWE 1 \
-doSnpStat 1 \
-sites sites.file

#################################### DONE ####################################
