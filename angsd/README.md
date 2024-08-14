- This directory contains all the job scripts used to run ANGSD (Korneliussen, Albrechtsen & Nielsen 2014) on an HPC. ANGSD is a 
method that estimates genotype likelihoods (GLs) and allele frequencies. Estimated GLs can be used by other programmes in the ANGSD 
family (e.g., PCAngsd, NGSadmix and NGSrelate) to preform various population genomic analyses.
- ANGSD is run at each filtering stage to generate files required for the next stage.
- After filtering, ANGSD is run again for population structure analyses and to estimate population allele frequencies for the BayPass 
analysis.

./exome/scripts
- This contains all the job scripts for running ANGSD on the exome data.

./exome/scripts/run_index.BAMs.sh
- This scripts indexes the scripts to prepare them for angsd.

./exome/scripts/run_angsd_f*
- These scripts run ANGSD for the exome data.
- File names:
	• As is made clear in ../sample_filtering/exome/scripts, I label each sample filtering stage as fx, i.e., unfiltered = f0, one 
filter stage = f1 etc.. The final filter stage for the exome is f5.
	• 0.5x refers to the fact that the minimum average coverage filter is 0.5-fold.
	• "all" means all samples were given together, "subsp" means ANGSD was run for each subspecies separately, "coms" means ANGSD 
was run for each sample site (aka. community) separately and "pops" means ANGSD was run for each population (defined in 
../population_structure/exome/scripts/population_structure.Rmd) separately.
	• Other aspects of file names correspond to ANGSD parameters and should be self explanatory.

./exome/scripts/run_realSFS.subsp.fst_2024.sh
- Calculates Fst between subspecies working on PanAf genotype likelihoods estimated with ANGSD

./chr21/scripts
- This contains all the job scripts for running ANGSD on the chr21data.

./chr21/scripts/run_angsd_chr21.f*
- These scripts run ANGSD for the exome data.
- The file name convention is the same as for the exomes, however, the final filtered stage for chr21 data is chr21.f7 (this is all 
explained in ../sample_filtering/chr21/scripts)

./$DATA/bam.filelists
- This is where the lists of bamfiles to be used in each analysis are stored. The subdirectories correspond to different filter 
levels.

./$DATA/output
- This is where ANGSD outputs are stored.

