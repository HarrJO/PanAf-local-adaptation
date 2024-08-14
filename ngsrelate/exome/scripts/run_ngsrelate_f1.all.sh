#!/bin/bash -l

# Job name
#$ -N ngsrelate_f1.coms_doMajorMinor.1_snp.p.1e-6_minMAF.5pct

# Set working directory
#$ -wd /home/ucfajos/Scratch/output/phase1and2_exome_output/ngsrelate_output/mapped.on.target/f1.coms_doMajorMinor.1_snp.p.1e-6_minMAF.5pct

# Time requested
#$ -l h_rt=01:00:00

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
#################################### JOB ####################################

INDIR='/home/ucfajos/Scratch/output/phase1and2_exome_output/angsd_output/mapped.on.target/f1.coms_doMajorMinor.1_snp.p.1e-6_minMAF.5pct_glf'
BAMLIST_DIR='/home/ucfajos/analysis/phase1and2_exome_analysis/angsd/bam.filelists/mapped.on.target/f1_pc1_hc.1pct_coverage.0.5x/'

for POP in c.Bateke c.CampoMaan c.Conkouati c.Goualougo c.Invindo c.LaBelgique c.Loango c.Lope c.MtsdeCristal e.Bili e.Budongo e.Bwindi e.Chinko e.Gishwati e.IssaValley e.Ituri e.Kabogo e.Ngiri e.Ngogo e.Nyungwe e.Regomuki e.RubiTele n.Gashaka n.Korup n.Mbe n.MtCameroon w.Azagny w.Bafing w.Bakoun w.Bia w.Boe w.BoundialeOdienne w.Comoe2 w.ComoeCNPN w.ComoeEAST w.ComoeGEPRENAF w.ComoeWEST w.Dindefelo w.Djouroutou w.EastNimba w.Grebo w.Kayan w.Loma w.MtSangbe w.OutambaKilimi w.Sangaredi w.Sapo w.Sobeya w.Sobory w.TaiEco w.TaiR w.ZooAnkasa
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
