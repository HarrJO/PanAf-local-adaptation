#!/bin/bash -l

# qsub.run_baypass_chr21_non-genic.rm_e.IssaValley_aux.sh
# 29/03/2022

COV='f_over_sum_known_trees.rm_e.IssaValley'
FLANKS=1000
INFILE=run_baypass_chr21_non-genic.rm_e.IssaValley_aux.sh


for SUBSP in ce
do
        if [[ $SUBSP == all ]]
        then
                POPS=chr21.f7.pops
        else
                POPS=chr21.f7.${SUBSP}.pops
        fi

	# Make output directory (if it already exists this will do nothing)
        mkdir -p /home/ucfajos/Scratch/output/phase1and2_chr21_output/baypass_output/ac/${POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_${FLANKS}bp.flanks.rm_e.IssaValley/aux/${COV}

        for SEED in 5001 100 10000
        do

                if [[ $SEED == 5001 ]]
                then
                        SEEDNAME=''
                fi

                if [[ $SEED == 10000 ]]
                then
                        SEEDNAME=_seed.10k
                fi

                if [[ $SEED == 100 ]]
                then
                        SEEDNAME=_seed.100
                fi
	
		# Edit script
		sed "s/{COV}/$COV/g" ${INFILE} | \
		sed "s/{FLANKS}/$FLANKS/g" | \
		sed "s/{SUBSP}/$SUBSP/g" | \
	        sed "s/{POPS}/$POPS/g" | \
		sed "s/{SEED}/$SEED/g" | \
		sed "s/{SEEDNAME}/$SEEDNAME/g" \
		> ${SUBSP}/run_baypass_${POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_${FLANKS}bp.flanks.rm_e.IssaValley_aux.${COV}.sh

                # !!! QSUB !!!
                qsub -m as ${SUBSP}/run_baypass_${POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_${FLANKS}bp.flanks.rm_e.IssaValley_aux.${COV}.sh
        done

done

