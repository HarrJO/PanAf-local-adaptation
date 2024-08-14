- BayPass (Gautier 2015) can be used to for evidence of local adaptation in the form of exceptionally high allele frequency 
differentiation (under the core model) or allele frequencies correlating with environmental covariables (under the AUX model).

./baypass/running_baypass
- Contains the job scripts used to run BayPass on a HPC.

./baypass/running_baypass/exome/scripts
- Directories correspond to subspecies-datasets. Each contains job scripts used to run BayPass.
- BayPass is first run using SNPs with no missing data under the core model to estimate the population allele frequency covariance 
matrix.
- BayPass is then run using SNPs with up to 30% missing data and using the estimated omega under the core model (to identify 
signatures of selection in the form of especially high allele frequency differentiation) or the AUX model (to perform a 
genotype-environment association analysis).
	• Because Nigeria-Cameroon only has two populations, BayPass was only run on this subspecies-dataset under the core model with 
no missing data.

./baypass/running_baypass/exome/scripts/run_baypass_aux.sh
- A generic job script for running the AUX model

./baypass/running_baypass/exome/scripts/qsub.run_baypass_aux.sh
- Shell script used to submit multiple jobs to the HPC under the AUX model by editing the generic AUX script 
(./baypass/running_baypass/exome/scripts/run_baypass_aux.sh). I wrote this after finishing the analysis under the core model so an 
equivalent script doesn't exist for the core model. Core model runs were done by manually editing the shell scripts and submitting 
them to the cluster.

./baypass/running_baypass/exome/output
- Contains all BayPass output separated into directories corresponding to particular subspecies-datasets and amount of missing data 
allowed (none or 30%).
	• Subdirectories separate into results under the core or AUX (auxi) model.
		○ Subdirectories in auxi correspond to the particular environmental covariable(s) used.
- Alongside output files are files which report the script and console output from running BayPass.

./baypass/running_baypass/chr21
- This directory mirrors ./baypass/running_baypass/exome only applied to non-genic-chr21 rather than the exome.

./baypass/running_baypass/bin/baypass_2.2
- Contains the binary files for BayPass.

./baypass/running_baypass/baypass_env_input
- Contains the environmental covariate files required for the AUX model (generated in 
../environmental_data/scripts/2.format_env_file.Rmd).



./baypass/analysing_baypass_output
- Contains Rmarkdowns and R scripts used to analyse the BayPass outputs. These are heavily annotated.

./baypass/analysing_baypass_output/baypass_core
- All related to the BayPass core model

./baypass/analysing_baypass_output/baypass_core/scripts/format_baypass_core_output.Rmd
- This script should be run first, it formats the BayPass core outputs (./baypass/running_baypass/$DATA/output) for the exome and for 
non-genic-chr21. This is only done for the runs which allowed up to 30% missing data at a SNP (apart from Nigeria-Cameroon which has 
only two populations anyway).
	• Adds the genomic position of each SNP (i.e. chromosome and position).
	• XtX files.
		○ Combines results from the three independent BayPass runs and calculates a median XtX value for each SNP.
		○ Adds coverage depth to the data frame.
		○ Assigns SNPs to genes using GTF file from ../gowinda/baypass_core.
	• Allele frequency files (exome only).
		○ Adds population names.
		○ Only a single file from a 'focal run' is used, i.e., the average allele frequencies across multiple independent 
runs is not used.
- This script is slow (I later wrote the equivalent for the AUX output in python to speed things up) but is necessary for all further 
analysis.
- Outputs to ./baypass/analysing_baypass_output/baypass_core/outputs/formatted_baypass_core_output and Outputs to 
./baypass/analysing_baypass_output/baypass_core/outputs/annotation.

./baypass/analysing_baypass_output/baypass_core/scripts/baypass_core_candidate.Rmd
- This script selects candidates using the formatted BayPass output from the exome and chr21 generated in 
./baypass/analysing_baypass_output/baypass_core/scripts/format_baypass_core_output.Rmd.
- This is done by using non-genic-chr21 as a null distribution to define estimated FPR thresholds for XtX values within coverage bins 
to account for this potentially confounding effect.
- Outputs to ./baypass/analysing_baypass_output/baypass_core/outputs/baypass_core_output_with_fprs.

./baypass/analysing_baypass_output/baypass_core/scripts/baypass_core_candidate_genome_distribution.Rmd
- This script plots the distribution of candidate SNPs across the genome.

./baypass/analysing_baypass_output/baypass_core/scripts/baypass_core_candidate_allele_frequency_patterns.Rmd
- This script investigates allele frequency patterns at candidate SNPs, mainly looking for evidence of clusters of populations/SNPs 
with similar patterns.

./baypass/analysing_baypass_output/baypass_core/scripts/baypass_core.top_n_genes.Rmd
- This script focuses on the putative functions of the genes with the very strongest evidence of selection.
- Outputs to ./baypass/analysing_baypass_output/baypass_core/outputs/top_snps.

./baypass/analysing_baypass_output/baypass_core/scripts/*.R
- This script define many custom functions used in these Rmarkdowns.


./baypass/analysing_baypass_output/baypass_aux
- All related to the BayPass AUX model

./baypass/analysing_baypass_output/baypass_aux/scripts/format_baypass_aux_output.ipynb
- Jupyter Notebook which formats BayPass AUX output (./baypass/running_baypass/$DATA/output) in the same way as the core output is 
formatted
	• Adds the genomic position of each SNP (i.e. chromosome and position).
	• Combines results from the three independent BayPass runs and calculates a median Bayes Factor (BF) and corelation 
coefficients (beta) value for each SNP.
	• Adds coverage depth to the data frame.
- This notebook uses the formatted core output to add gene annotations (so that must be run first).
- The script is written in a way to accommodates a range of different ways of conducting the BayPass analysis which were subsequently 
abandoned in favour of the current method.
- Outputs to ./baypass/analysing_baypass_output/baypass_aux/output/formatted_baypass_aux_output

./baypass/analysing_baypass_output/baypass_aux/scripts/baypass_aux_tools.R
- This contains custom R functions used in Rmarkdown scripts.


./baypass/analysing_baypass_output/baypass_aux/scripts/f_over_sum_known_trees/f_over_sum_known_trees.baypass_aux_distributions.Rmd
- This scrip plots the distributions of BayPass AUX output statistics.

./baypass/analysing_baypass_output/baypass_aux/scripts/f_over_sum_known_trees/f_over_sum_known_trees.baypass_aux_candidate.Rmd
- This script selects candidates using the formatted BayPass output from the exome and chr21 generated in 
./baypass/analysing_baypass_output/baypass_aux/scripts/format_baypass_aux_output.ipynb.
- This is done as for the core model, by using non-genic-chr21 as a null distribution to define estimated FPR thresholds for BF values 
within coverage bins to account for this potentially confounding effect.
- Outputs to ./baypass/analysing_baypass_output/baypass_aux/output/baypass_aux_output_with_fprs.

./baypass/analysing_baypass_output/baypass_aux/scripts/f_over_sum_known_trees/f_over_sum_known_trees.baypass_aux_distribution_of_candidates_in_the_genome.Rmd
- This script plots the distribution of candidate SNPs across the genome.

./baypass/analysing_baypass_output/baypass_aux/scripts/f_over_sum_known_trees/f_over_sum_known_trees.baypass_aux.candidate_allele_frequencies.Rmd
- This script investigates allele frequency patterns at candidate SNPs, mainly looking for evidence of clusters of populations/SNPs 
with similar patterns.

./baypass/analysing_baypass_output/baypass_aux/scripts/f_over_sum_known_trees/f_over_sum_known_trees.baypass_aux.top_n_genes.Rmd
- This script focuses on the putative functions of the genes with the very strongest evidence of selection.

./baypass/analysing_baypass_output/baypass_aux/scripts/f_over_sum_known_trees/f_over_sum_known_trees.baypass_aux_ce_rm_e.IssaValley.Rmd
- This script analyses the outputs for the test run with the Issa Valley population removed to determine the effect this population 
has on the candidates.


./baypass/analysing_baypass_output/scripts
- Contains Rmarkdowns which produce various plots related to the core and AUX models

./baypass/analysing_baypass_output/scripts/baypass_output_summary.Rmd
- This script plots various figure giving an overview of the BayPass analysis.
- Many of these plots are similar to those in there scripts only formatted for publication

./baypass/analysing_baypass_output/scripts/baypass_tools.R
- This contains custom R functions used in Rmarkdown scripts for both the core and AUX models.

./baypass/analysing_baypass_output/scripts/baypass_output_exome_vs_chr21.Rmd
- This script compared the BayPass core and AUX results between the exome and non-genic-chr21
- This is similar in spirit to the methodology used to define FDR thresholds but looks across multiple bins in an attempt to do something similar to the DAF difference work in Schmidt et al. (2019) and ../allele_frequencies/scripts/plot.exome.and.chr21.SFS.Rmd

./baypass/analysing_baypass_output/scripts/number_of_SNPs_and_genes.Rmd
- This script looks at the number of SNPs and genes covered by the PanAf exomes and high coverage genomes 

./baypass/analysing_baypass_output/scripts/baypass_output_core_vs_aux.Rmd
- This script compares the values per SNP from the core and aux analysis 

./baypass/analysing_baypass_output/malaria_candidates/scripts
- This contains scripts for analysing key malaria related candidate genes identified under the AUX model (HBB/HBD and GYPA/GYPB).

./baypass/analysing_baypass_output/malaria_candidates/scripts/baypass_output_HBB-HBD_GYPA-GYPB.Rmd
- This Rmarkdown looks at the evidence of selection at these loci (HBB/HBD and GYPA/GYPB) in the BayPass core and AUX model tests 

./baypass/analysing_baypass_output/malaria_candidates/scripts/structure_and_function_HBB-HBD_GYPA-GYPB.Rmd 
- This Rmarkdown investigates the function of these loci (HBB/HBD and GYPA/GYPB) and investigates evidence of structural variation 
using high coverage genomes aligned to hg19  (de Manuel et al. 2016; Prado-Martinez et al. 2013) 

