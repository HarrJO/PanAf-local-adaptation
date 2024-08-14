- ANGSD is run population-by-population (and chromosomes-by-chromosome) to estimate population minor allele frequencies (job 
scripts in ./angsd); this directory contains custom python scripts and Rmarkdowns that collate all these files, reformat them and 
perform SNP filtering.
- The python scripts were run on a HPC and run chromosomes-by-chromosome (for the exome) so the process could be run in 
parallel. The large ANGSD output files means this process is very slow.
	
- From here on, analyses are run on sperate "subspecies-datasets" all subspecies (All), (2) central and eastern together 
(Central-Eastern) (3) Nigeria-Cameroon (Nigeria-Cameroon) and (4) western (Western). The central and eastern subspecies were combined 
because of their very low levels of genetic differentiation.
- The following prefixes are used to denote these datasets (exome/chr21):
		○ All: f5.pops/chr21.f7.pops
		○ Central-Eastern: f5.ce.pops/chr21.f7.ce.pops
		○ Nigeria-Cameroon: f5.n.pops/chr21.f7.n.pops
		○ Western: f5.w.pops/chr21.f7.w.pops

./scripts/
- This directory contains the scripts relevant for both the exome and chr21 data.
- Shell scripts are job scripts which apply the custom python scripts in their name (e.g. run_minMAC.py.sh just runs minMAC.py 
on a HPC).
	
./scripts/maf2dac.v3.2.py
- This script takes the allele frequencies estimated by running ANGSD for each population and chromosomes (.maf.gz files), 
filters sites and produces a population allele count dataset for input into BayPass.
- This is a very slow process as the .maf.gz files are very big. This is because these files are not filtered to only include 
polymorphic sites because sites which are locally monomorphic (i.e. monomorphic within a population) may prove to be globally 
polymorphic (i.e. polymorphic when all populations are considered).
		○ A better way of doing this in the future might be to run ANGSD on all samples togter to identify sites with are 
globally polymorphic then apply ANGSD only to these regions using the -sites flag.
- The script is heavily annotated with help messages for each flag, but in summary:
		○ Populations from specific subspecies can be selected.
		○ Specific populations can be excluded (e.g.those with a sample size <8).
		○ SNP p-value thresholds can be selected to define SNPs as polymorphic within a population.
		○ The number of populations allowed to have missing data at a SNP can be specified, I run it once allowing no missing 
data (for estimating the covariance matric under the BayPass core model) and one allowing up of 30% of populations to have missing 
data (for testing for signatures of selection as BayPass is robust to some missing data at this stage).
			§ A quicker way of doing this would be to run it allowing 30% missing data and then just filter SNPs 
afterwards.

./scripts/minMAC.py
- This script filters files outputted by./scripts/maf2dac.v3.2.py using a minimum minor allele count (minMAC) threshold. I use 
minMAC 2 which removes all SNPs with a MAC <2 across all populations in a subspecies-dataset.
	
./scripts/rm_pops.py
- This script simply removes columns corresponding to particular populations, this is used to remove Issa Valley and run tests 
to assess the influence of this one population on the BayPass AUX results.

./scripts/subset_baypass_input.py
- This script can be used to subset the BayPass input file if BayPass needs to be run multiple times in parallel due to 
computational constraints.
- Splitting into two subsets works by selecting every other SNP for one and the rest for another, three subsets involves 
selecting everything third SNP and so on. This method means the subsets are spread across the genome.
- For me, this was only required when running BayPass with multiple environmental variables at once.

./scripts/allele_frequency_tools.R
- Defines custom functions used in Rmarkdown scripts 

./exome/scripts/
- Contains subdirectories with job scripts used to run ./scripts/maf2dac.v3.2.py on the exome data.

./chr21/scripts/
- Contains subdirectories with job scripts used to run ./scripts/maf2dac.v3.2.py on the chr21data.

./chr21/scripts/filter_non-geneic_SNPs.Rmd
Rmarkdown that filters population allele count files to only include non-genic regions of chr21 (non-genic-chr21)
