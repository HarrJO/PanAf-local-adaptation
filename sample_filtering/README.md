- This directory contains Rmarkdowns that plot and interrogate the results from PCAngsd/NGSadmix/NGSrelate and filter out 
samples with exceptionally low coverage, evidence of contamination or related samples.
- Exome data filtering is performed first, chr21 samples are then filtered for only those that passed the exome filtering and 
do not show any evidence of contamination in the chr21 data.

- File names:
		○ I label each sample filtering stage as fx, i.e., unfiltered = f0, one filter stage = f1 etc.. 
		○ The final filter stage for the exome is f5 and for the chr21 it is f7.

./$DATA/scripts
- Rmarkdowns which analyse and plot the results from PCAngsd/NGSadmix/NGSrelate.
- These are thoroughly annotated.

meta.data/
- This meta data will be included as a supplementary file in the Ostridge et al. publication