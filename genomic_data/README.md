This directory contains the genomic data used.

./bams
- This directory is where the PanAf BAMs for chr21 and exome data containing only mapped on target reads (mapped to hg19) can be stored to run ANGSD. 
- Exome BAMs will be publicly available on publication on ENA under the accession code ENA:PRJEB76176
- Chr21 BAMs are publicly available on ENA under the accession code ENA:PRJEB46115

./ref_genomes
- This directory is where the reference genome and ancestral ancestral state files can be stored to run ANGSD.
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
