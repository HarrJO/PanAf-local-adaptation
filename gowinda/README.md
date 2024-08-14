- Gowinda is used to perform gene set enrichment analysis on the candidate SNPs identified in BayPass.
- This directory contains all scripts required to run gowinda and analyse outputs.


./gowinda/scripts/run_gowinda.sh
- This script is the general gowinda script used.

./gowinda/scripts/bin
- Contains all the gowinda binary files.

./gowinda/baypass_core
- Contains the scripts required to generate input files, run gowinda and analyses the output. 

./gowinda/baypass_core/scripts/run_gowinda.Rmd
- Rmarkdown setting up and running myriad.
- This must be run first to generate all the appropriate files (for both core and AUX candidates).
- ../baypass/analysing_baypass_output/baypass_core/scripts/baypass_core_candidate.Rmd must be run to generate candidate SNP files.

./gowinda/baypass_core/scripts/gowinda_outputs.Rmd
- Rmarkdown plotting gowinda outputs

./gowinda/baypass_core/scripts/plot_gowinda.R
- R script defining custom functions for plotting gowinda outputs.


./gowinda/baypass_aux
- Contains the scripts required to run gowinda and analyses the output - this directory used input files generates in 
./gowinda/baypass_core so some of these scripts must be run first

./gowinda/baypass_aux/scripts/run_gowinda.aux.general.Rmd
- Rmarkdown that runs gowinda for the candidates identified under the AUX model.
- ./gowinda/baypass_core/scripts/run_gowinda.Rmd must be run first to generate the appropriate annotation and association files.
- ../baypass/analysing_baypass_output/baypass_aux/scripts/f_over_sum_known_trees/f_over_sum_known_trees.baypass_aux_candidate.Rmd must 
be run to generate candidate SNP files.

./gowinda/baypass_aux/scripts/gowinda_outputs.aux.f_over_sum_known_trees.paper_figures.Rmd
- Plots the gowinda output for the AUX candidates.

