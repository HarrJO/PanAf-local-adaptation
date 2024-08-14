#!/bin/bash -l

# qsub.run_baypass_chr21_non-genic.miss_pops_0.3.core.sh
# 23/09/2022

FLANKS=1000bp

for SUBSP in ce w
do
        sed "s/{SUBSP}/$SUBSP/g" run_baypass_chr21_non-genic.miss_pops_0.3.core.sh | sed "s/{FLANKS}/$FLANKS/g" > run_baypass_TEMP1.sh

        if [[ $SUBSP == all ]]
        then
                POPS=chr21.f7.pops
        else
                POPS=chr21.f7.${SUBSP}.pops
        fi

        sed -i "s/{POPS}/$POPS/g" run_baypass_TEMP1.sh

	# Make output directory (if it already exists this will do nothing)
        mkdir /home/ucfajos/Scratch/output/phase1and2_chr21_output/baypass_output/ac/${POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_${FLANKS}.flanks
	mkdir /home/ucfajos/Scratch/output/phase1and2_chr21_output/baypass_output/ac/${POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_${FLANKS}.flanks/core

        for SEED in 5001 100 10000
        do

                sed "s/{SEED}/$SEED/g" run_baypass_TEMP1.sh > run_baypass_TEMP2.sh

                if [[ $SEED == 5001 ]]
                then
                        sed -i "s/{SEEDNAME}//g" run_baypass_TEMP2.sh
                fi

                if [[ $SEED == 10000 ]]
                then
                        sed -i "s/{SEEDNAME}/_seed.10k/g" run_baypass_TEMP2.sh
                fi

                if [[ $SEED == 100 ]]
                then
                        sed -i "s/{SEEDNAME}/_seed.100/g" run_baypass_TEMP2.sh
                fi

                cp run_baypass_TEMP2.sh ${SUBSP}/run_baypass_${POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_${FLANKS}.flanks_core.sh

                # !!! QSUB !!!
                qsub -m as ${SUBSP}/run_baypass_${POPS}_minInd.6or50pct_missing.pops.0.3_minMAC2.non-genic_${FLANKS}.flanks_core.sh
        done

done

rm run_baypass_TEMP?.sh
