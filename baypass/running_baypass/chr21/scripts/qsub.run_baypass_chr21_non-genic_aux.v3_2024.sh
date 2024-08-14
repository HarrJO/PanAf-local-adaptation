#!/bin/bash -l

# qsub.run_baypass_chr21_non-genic_aux.sh
# 29/03/2022

COV='paper_2_beh'
FLANKS=1000
INFILE=run_baypass_chr21_non-genic_aux.v3_2024.sh
HOURS=48 # Values <10 should have leading zeros, e.g. 4 should be 04


for SUBSP in all ce w
do
        if [[ $SUBSP == all ]]
        then
                POPS=chr21.f7.pops
        else
                POPS=chr21.f7.${SUBSP}.pops
        fi
	
	if [[ $SUBSP == n ]]
        then
                MISS_POPS1=0
                MISS_POPS2=0.0
        else
                MISS_POPS1=0.3
                MISS_POPS2=0.3
        fi

	# Make output directory (if it already exists this will do nothing)
        mkdir -p /home/ucfajos/Scratch/output_2024/phase1and2_chr21_output/baypass_output/ac/${POPS}_minInd.6or50pct_missing.pops.${MISS_POPS1}_minMAC2.non-genic_${FLANKS}bp.flanks/auxi/${COV}

        #for SUBSET in _subset1of4 _subset2of4 _subset3of4 _subset4of4
        #for SUBSET in _subset1of2 _subset2of2
        for SUBSET in ''
        do

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
			sed "s/{SEEDNAME}/$SEEDNAME/g" | \
			sed "s/{MISS_POPS1}/$MISS_POPS1/g" | \
	                sed "s/{MISS_POPS2}/$MISS_POPS2/g" | \
			sed "s/{SUBSET}/$SUBSET/g" | \
	                sed "s/{HOURS}/$HOURS/g" \
			> ${SUBSP}/run_baypass_${POPS}_minInd.6or50pct_missing.pops.${MISS_POPS1}_minMAC2.non-genic_${FLANKS}bp.flanks_aux.${COV}${SUBSET}.sh

	                # !!! QSUB !!!
	                #qsub -m as ${SUBSP}/run_baypass_${POPS}_minInd.6or50pct_missing.pops.${MISS_POPS1}_minMAC2.non-genic_${FLANKS}bp.flanks_aux.${COV}${SUBSET}.sh
		done
        done

done

