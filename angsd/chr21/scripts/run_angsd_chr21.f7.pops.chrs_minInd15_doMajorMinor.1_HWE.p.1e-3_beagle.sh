#!/bin/bash -l

# Job name
#$ -N angsd_chr21.f7.pops_minInd.6or50pct_doMajorMinor.1_HWE.p.1e-3_beagle

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_chr21_output/angsd_output/chr21.f7.pops.chrs_minInd.6or50pct_doMajorMinor.1_HWE.p.1e-3_beagle

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
# I am ruuning this on populations rather than communities as some communities were merged after population structure analysis 
# This is run chr by chr so the next stage i.e. estimateing allele counts and filtering sites (maf2dac) can be run chr by chr for speed
#################################### JOB ####################################

# Load modules
module load htslib

# Subspecies; c=central, e=eastern, n=Nigeria-Cameroon, w=western
for POP in c.Bateke c.Conkouati c.Goualougo c.LaBelgique c.Loango c.Lope c.MtsdeCristal e.Bili e.Budongo e.Bwindi e.Chinko e.Gishwati e.IssaValley e.Kabogo e.Ngogo e.Nyungwe e.Regomuki e.RubiTele n.Gashaka n.KorupMtCameroon n.Mbe w.Bafing w.BakounSobory w.Bia w.Boe w.Comoe w.Dindefelo w.Djouroutou w.EastNimba w.Grebo w.Kayan w.MtSangbe w.OutambaKilimi w.Sangaredi w.Sapo w.Sobeya w.Tai
	do
	# Number of individuals from number of bam files in bamfile list
	NIND=`wc -l < /home/ucfajos/analysis/phase1and2_chr21_analysis/angsd/bam.filelists/f7/bam.filelist.${POP}`
	# minInd is 50% of the nInd or 6, whichever is higher
	MININD=`awk "BEGIN {if ((${NIND}/2)>6) print ${NIND}/2; else print 6}"`
		/usr/bin/time --verbose \
		~/bin/angsd/angsd \
		-b /home/ucfajos/analysis/phase1and2_chr21_analysis/angsd/bam.filelists/f7/bam.filelist.${POP} \
		-out ${POP}/${POP}.chr21 \
		-ref ~/Scratch/data/ref_genomes/hg19.fa \
		-anc ~/Scratch/data/ref_genomes/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_fullgenome.fa \
		-uniqueOnly 1 \
		-remove_bads 1 \
		-only_proper_pairs 1 \
		-trim 0 \
		-C 50 \
		-baq 1 \
		-minInd $MININD \
		-skipTriallelic 1 \
		-GL 2 \
		-minMapQ 30 \
		-nThreads 10 \
		-doMajorMinor 1 \
		-doMaf 2 \
		-SNP_pval 1 \
		-doSnpStat 1 \
		-doHWE 1 \
		-minHWEpval 0.001 
	done






#################################### DONE ####################################
