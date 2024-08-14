#!/bin/bash -l

# Job name
#$ -N angsd_chr21.f7.subsp_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01_beagle

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_chr21_output/angsd_output/chr21.f7.subsp_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01_beagle

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
# The output is used to plot a PCA and perform the first stage of filtering (removing PC1 outliers)
#################################### JOB ####################################

# Load modules
module load htslib

# Subspecies; c=central, e=eastern, n=Nigeria-Cameroon, w=western
for POP in c e n w
	do
	/usr/bin/time --verbose \
	~/bin/angsd/angsd \
	-b /home/ucfajos/analysis/phase1and2_chr21_analysis/angsd/bam.filelists/f7/bam.filelist.${POP} \
	-out ${POP}/chr21.f7.${POP} \
	-ref ~/Scratch/data/ref_genomes/hg19.fa \
	-anc ~/Scratch/data/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa \
	-uniqueOnly 1 \
	-remove_bads 1 \
	-only_proper_pairs 1 \
	-trim 0 \
	-C 50 \
	-baq 1 \
	-minInd 15 \
	-skipTriallelic 1 \
	-GL 2 \
	-doGlf 2 \
	-minMapQ 30 \
	-nThreads 10 \
	-doMajorMinor 1 \
	-doMaf 2 \
	-SNP_pval 0.000001 \
	-minMaf 0.01 \
	-doHWE 1 \
	-minHWEpval 0.001
	done


#################################### DONE ####################################
