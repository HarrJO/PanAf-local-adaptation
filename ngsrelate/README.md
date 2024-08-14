- NGSrelate (Hanghøj et al. 2019) is used to estimate kinship coefficients using GLs estimated in ANGSD. 
- This script contains bash job scripts used to run this on a HPC.
- This is used in sample filtering to remove related samples.

- This was run only for the exome (chr21 data for samples were filtered based on the kinship coefficients estimated from the exome 
data).

./exome/scripts/run_ngsrealte_f*
- File names
	• As is made clear in ../sample_filtering/exome/scripts, I label each sample filtering stage as fx, i.e., unfiltered = f0, one 
filter stage = f1 etc.. The final filter stage for the exome is f5.
"all" means ngsrelate is run for all sample sites (aka communities). If sample sites are names then the analysis is only run for these