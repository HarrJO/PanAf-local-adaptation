This directory contains the genomic data used.

./bams
- Contains the exome and chr21 mapped BAM files mapped to hg19 containing only on target reads.
- It may not actually contain them due to a lack of storage space, if so, a README in that directory will say where to find them.
- Each file corresponds to a single sample.
- Metadata describing each sample can be found in ../sample_filtering

./ref_genomes
- Contains the reference genomes and ancestral state files provided to ANGSD.
- ./ref_genomes/scripts contains the scripts used to format and index these files.
- Details for downloading, formatting and indexing these files are found below

# Reference file
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

## Indexing genome
### Submitted 'run_prepare_ref_genome.sh'

### Load modules
module load bwa
module load samtools

### bwa index
bwa index -a bwtsw hg19.fa

### samtools fadix
samtools faidx hg19.fa 

# Ancestral file
cd ~/Scratch/ref_genomes
wget ftp://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
wget ftp://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.README
wget ftp://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.MD5SUM
tar -xvjf homo_sapiens_ancestor_GRCh37_e71.tar.bz2
cd homo_sapiens_ancestor_GRCh37_e71
mkdir GL000_files
mv *GL000* GL000_files/

### I then used nano to edit all the fasta headers into the format "chr1" etc

### Submitted 'run_prepare_ref_genome.sh'

cat homo_sapiens_ancestor_*.fa > homo_sapiens_ancestor_fullgenome.fa
bwa index -a bwtsw homo_sapiens_ancestor_fullgenome.fa
samtools faidx homo_sapiens_ancestor_fullgenome.fa


./wANNOVAR _output
- Contains the output from running wANNOVAR (https://wannovar.wglab.org/) using the candidate SNPs identified in GYPA/GYPB and HBB/HBD identified in the BayPass AUX analysis.

./schmidt_et.al.2019
- Contains data downloaded from this publication which is used to select 1-1 homologs for the gowinda analysis.
- (https://datadryad.org/stash/dataset/doi:10.5061/dryad.zcrjdfn6m)

./pawar.et.al.2022
- Contains data generated for this publication. This data is used to look at PBSnj and 3P-CLR values at key malaria candidates in ./baypass/analysing_baypass_output/malaria_candidates/scripts/.

./exome_coordinates
- Contains bed files defining the exome target space.

./gene_cards
- Used to search for the functions of the few genes with the strongest evidence of selection.
- I started a conversation with the people at GeneCards through the chat window at https://www.genecards.org/ (you can also email support@genecards.org). It took a long time but I finally got them to send me a csv with gene names and descriptions of each gene. This means that rather than manually searching each gene name in the web interface I can just search this file.

./gene_expression_data
- Data with gene expresion results sued to investigate expression levels of key malaria related genes in ./baypass/analysing_baypass_output/malaria_candidates/scripts/.

./gene_expression_data/Fair_et_al._2020
- Gene expression data from primary heart samples.
- "39 post-mortem heart tissue biopsies were collected from captive born chimpanzees, 18 of which have been previously described (Pavlovic et al., 2018)" - there are 39 sample columns 

./high_coverage_genomes/hg19_mapped_VCFs
- Contains data from de Manuel et al. (2016) and Prado-Martinez et al. (2013) aligned to hg19.
- Details of how the HBB/HBD and GPYA/GYPB specific files were generate are in ./baypass/analysing_baypass_output/malaria_candidates/scripts/structure_and_function_HBB-HBD_GYPA-GYPB.Rmd.

./high_coverage_genomes/scripts/running_GATK.Rmd
- This script shows how I generated the VCF from BAM files aligned to hg19.

