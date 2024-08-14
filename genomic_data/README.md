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
