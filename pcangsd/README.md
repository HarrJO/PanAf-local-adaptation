- PCAngsd (Meisner and Albrechtsen 2018) is used to perform PCA using GLs estimated in ANGSD. 
- This script contains bash job scripts used to run this on a HPC.
- This is used in sample filtering to identify sample swaps and outliers indicative of contamination (./sample_filtering) and for 
population structure analyses (./population_structure).

./exome/scripts/run_pcangsd_f*
- These scripts run PCAngsd for the exome data.
- File names:
	- As is made clear in ../sample_filtering/exome/scripts, I label each sample filtering stage as fx, i.e., unfiltered = f0, one 
filter stage = f1 etc.. The final filter stage for the exome is f5.
	- "all" means all samples were given together and "subsp" means ANGSD was run for each subspecies separately.

./chr21/scripts/run_pcangsd_f*
- These scripts run PCAngsd for the chr21data.
The file name convention is the same as for the exomes, however, the final filtered stage for chr21 data is f7 (this is all explained 
in ../sample_filtering/chr21/scripts)
