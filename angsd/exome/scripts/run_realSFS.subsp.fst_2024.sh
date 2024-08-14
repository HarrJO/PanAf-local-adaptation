#!/bin/bash -l

# Job name
#$ -N subsp.fst

# Time requested
#$ -l h_rt=06:00:00

# RAM for each core/thread (integer followed by M, G, or T)
#$ -l mem=1G

# TMPDIR space needed (default is 10G)
#$ -l tmpfs=10G

# Set working directory
#$ -wd /home/ucfajos/Scratch/output_2024/phase1and2_exome_output/realsfs/mapped.on.target/f5.0.5x.subsp.fst

# Number of cores
#$ -pe smp 10

#################################### JOB ####################################

module load htslib/1.7

# ANGSD directory
ANGSD_DIR=/home/ucfajos/Scratch/output_2024/phase1and2_exome_output/angsd_output/mapped.on.target/f5.0.5x.subsp_minInd.15_doMajorMinor.1_HWE.p.1e-3_snp.p.1e-6_minMAF.0.01_beagle_dosaf_2024/

POP1=c
POP2=e
POP3=n
POP4=w

listA=(${POP1} ${POP1} ${POP1} ${POP2} ${POP2} ${POP3})
listB=(${POP2} ${POP3} ${POP4} ${POP3} ${POP4} ${POP4})

for i in {0..5}
	do
	echo ${listA[$i]}' vs. '${listB[$i]}
	
	# Calculate 2D SFS prior
	/usr/bin/time --verbose \
	~/bin/angsd/misc/realSFS ${ANGSD_DIR}/${listA[$i]}/f5.0.5x.${listA[$i]}.saf.idx ${ANGSD_DIR}/${listB[$i]}/f5.0.5x.${listB[$i]}.saf.idx >${listA[$i]}.${listB[$i]}.ml
	
	# Index
	/usr/bin/time --verbose \
	~/bin/angsd/misc/realSFS fst index ${ANGSD_DIR}/${listA[$i]}/f5.0.5x.${listA[$i]}.saf.idx ${ANGSD_DIR}/${listB[$i]}/f5.0.5x.${listB[$i]}.saf.idx \
	-sfs ${listA[$i]}.${listB[$i]}.ml \
	-fstout ${listA[$i]}.${listB[$i]}
	
	# Get global Fst estimate
	/usr/bin/time --verbose \
	~/bin/angsd/misc/realSFS fst stats ${listA[$i]}.${listB[$i]}.fst.idx > ${listA[$i]}.${listB[$i]}.global.fst
done

#################################### DONE ####################################
