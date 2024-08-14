#!/bin/bash -l

# run_rm_pops.py.sh
# 26/08/2022

# Load python
module load python3/3.8


python3 -u /home/ucfajos/analysis/phase1and2_exome_analysis/allele.frequencies/scripts/mapped.on.target/rm_pops.py \
-in_prefix /home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.ce.pops_minInd.6or50pct_missing.pops.0/f5.ce.pops.all.chrs_missing.pops.0 \
-rm_pops e.IssaValley

python3 -u /home/ucfajos/analysis/phase1and2_exome_analysis/allele.frequencies/scripts/mapped.on.target/rm_pops.py \
-in_prefix /home/ucfajos/Scratch/output/phase1and2_exome_output/allele.frequencies/mapped.on.target/f5.ce.pops_minInd.6or50pct_missing.pops.0.3/f5.ce.pops.all.chrs_missing.pops.0.3 \
-rm_pops e.IssaValley

