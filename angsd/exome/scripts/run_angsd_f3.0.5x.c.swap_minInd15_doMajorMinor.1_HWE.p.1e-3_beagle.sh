#!/bin/bash -l

# Job name
#$ -N angsd_f3.0.5x.c.swap_minInd.5_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.05_beagle

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f3.0.5x.c.swap_minInd.5_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.05_beagle

# Time requested
#$ -l h_rt=06:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=5G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 10

#################################### INFO ####################################
# This script is run with doMajorMinor 1 (to ensure no excess of high freq DAFs) and a HWE filter (to ensure no bump at 0.5 due to paralogs).
# GLF is in beagle format for PCAngsd
# I am running ANGSD for only the central sample which has been missassigned as Campo Maan and Conkouati and Loango as PCA and NGSadmix shows it is from one of these.
## I can then run PCAngsd and NGSadmix to work out which of the two communities the sample is from 
# I decrease minInd to 5 and increase minMaf to 0.05 to reflect the smaller sample size
#################################### JOB ####################################

# Load modules
module load htslib

# Run ANGSD
/usr/bin/time --verbose \
~/bin/angsd/angsd \
-b /home/ucfajos/analysis/phase1and2_exome_analysis/angsd/bam.filelists/mapped.on.target/f3_pc1_hc.1pct_coverage.0.5x_rm.rel_rm.outliers/bam.filelist.c.Conkouati.c.Loango.CMNP1-24 \
-out f3.0.5x.c.swap \
-ref ~/Scratch/data/ref_genomes/hg19.fa \
-anc ~/Scratch/data/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa \
-uniqueOnly 1 \
-remove_bads 1 \
-only_proper_pairs 1 \
-trim 0 \
-C 50 \
-baq 1 \
-minInd 5 \
-skipTriallelic 1 \
-GL 2 \
-doGlf 2 \
-minMapQ 30 \
-nThreads 20 \
-doMajorMinor 1 \
-doMaf 2 \
-SNP_pval 0.000001 \
-minMaf 0.05 \
-doHWE 1 \
-minHWEpval 0.001 

#################################### DONE ####################################
