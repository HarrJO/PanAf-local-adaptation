- NGSadmix (Skotte Korneliussen and Albrechtsen 2013) is used to estimate individual admixture proportions using GLs estimated in 
ANGSD. 
- This script contains bash job scripts used to run this on a HPC.
- This is used in sample filtering (to identify sample swaps and outliers indicative of contamination) and for population structure 
analyses.

- Unlike for PCAngsd, this was run only for the exome data as it was decided that NGSadmix does not provide much more information than PCAngsd.

./exome/scripts/f*/run_ngsadmix*
- These scripts run NGSadmix for the exome data.
- Directory names:
	• As is made clear in ../sample_filtering/exome/scripts, I label each sample filtering stage as fx, i.e., unfiltered = f0, one 
filter stage = f1 etc.. The final filter stage for the exome is f5.
