#!/bin/bash -l

# Job name
#$ -N angsd_f4.0.5x.c.Conkouati.w.MtSangbe_doMajorMinor.1_snp.p.1e-6_minMAF.5pct_glf

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f4.0.5x.c.Conkouati.w.MtSangbe_doMajorMinor.1_snp.p.1e-6_minMAF.5pct_glf

# Time requested
#$ -l h_rt=03:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=5G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10

#################################### INFO ####################################
# The purpose of this script to to make .glf.gz files for ngsRelate
# I have used exactly the same parmeters as Claudia to make my results comparable (emailed to me 10/03/2021)
# I am running ngsrelate again because swapped samples were reassigned to these populations and I want to check they are not relatives 
#################################### JOB ####################################

# Load modules
module load htslib

# Run ANGSD
for POP in c.Conkouati w.MtSangbe 
	do
	/usr/bin/time --verbose \
	~/bin/angsd/angsd \
	-b /home/ucfajos/analysis/phase1and2_exome_analysis/angsd/bam.filelists/mapped.on.target/f4_pc1_hc.1pct_coverage.0.5x_rm.rel_rm.outliers_sort.swaps/bam.filelist.${POP} \
	-out ${POP}/${POP} \
	-ref ~/Scratch/data/ref_genomes/hg19.fa \
	-anc ~/Scratch/data/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa \
	-uniqueOnly 1 \
	-remove_bads 1 \
	-only_proper_pairs 1 \
	-trim 0 \
	-C 50 \
	-baq 1 \
	-skipTriallelic 1 \
	-gl 2 \
	-minMapQ 30 \
	-nThreads 10 \
	-doGlf 3 \
	-doMajorMinor 1 \
	-doMaf 1 \
	-minMaf 0.05 \
	-SNP_pval 1e-6
	done

#################################### DONE ####################################
