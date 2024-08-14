#!/bin/bash -l

# Job name
#$ -N ngsrelate_f4.c.Conkouati.w.MtSangbe_doMajorMinor.1_snp.p.1e-6_minMAF.5pct

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/ngsrelate_output/mapped.on.target/f4.0.5x.c.Conkouati.w.MtSangbe_doMajorMinor.1_snp.p.1e-6_minMAF.5pct

# Time requested
#$ -l h_rt=00:30:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=1G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Number of cores
#$ -pe smp 1

#################################### INFO ###################################
# Prior to this, have we generated a file with allele frequencies (*.mafs.gz) and a file with genotype likelihoods (*.glf.gz)
# We are running ngsRelate community by community as substructure messes up the results (according to Claudia)
# The parameters used are exactly the same as Claudia sent me on 10/03/2021 (note the slight difference in that I need to do `-cut f7` rather than f6 as my .mafs.gz file 
## has a column for the ancestral allele too)
# This is the second time I have run this, I am only running it on populations where swapped individuals have been reassigned to check they are not relatives
#################################### JOB ####################################

INDIR='/home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f4.0.5x.c.Conkouati.w.MtSangbe_doMajorMinor.1_snp.p.1e-6_minMAF.5pct_glf'
BAMLIST_DIR='/home/ucfajos/analysis/phase1and2_exome_analysis/angsd/bam.filelists/mapped.on.target/f4_pc1_hc.1pct_coverage.0.5x_rm.rel_rm.outliers_sort.swaps'

for POP in c.Conkouati w.MtSangbe
	do
	# We extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
	zcat ${INDIR}/${POP}/${POP}.mafs.gz | cut -f7 | sed 1d >${POP}/${POP}.freq
	# Get the number of samples in the community
	NIND=`wc -l < ${BAMLIST_DIR}/bam.filelist.${POP}`
	# Run ngsRelate
	/usr/bin/time --verbose \
	/home/ucfajos/bin/ngsRelate/ngsRelate \
	-g ${INDIR}/${POP}/${POP}.glf.gz \
	-z ${BAMLIST_DIR}/bam.filelist.${POP} \
	-n ${NIND} \
	-f ${POP}/${POP}.freq \
	-O ${POP}/${POP}.ngsrelate.out
	done

#################################### DONE ####################################
