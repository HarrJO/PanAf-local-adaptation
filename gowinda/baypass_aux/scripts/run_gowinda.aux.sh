#!/bin/sh

# Read in command line arguments
while getopts ":b:c:g:a:o:" opt; do
        case $opt in
        b) BG_SNP="$OPTARG"
        ;;
        c) CAND_SNP="$OPTARG"
        ;;
        g) GENE_SET="$OPTARG"
        ;;
        a) ANNOTATION="$OPTARG"
        ;;
        o) OUTPUT="$OPTARG"
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
        ;;
        esac
done

# Run gowinda
## Remove log file if it already exists
rm ${OUTPUT}_log.txt
echo '### Background SNP file' ${BG_SNP} '\n' >${OUTPUT}_log.txt
echo '### Candidate SNP file' ${CAND_SNP} '\n' >>${OUTPUT}_log.txt
echo '### Gene set file ' ${GENESET} '\n' >>${OUTPUT}_log.txt
echo '### Annotation file ' ${ANNOTATION} '\n' >>${OUTPUT}_log.txt
echo '### Output file' ${OUTPUT} '\n' >>${OUTPUT}_log.txt
## Run gowinda
java -Xmx4g -jar bin/Gowinda-1.12.jar \
--snp-file ${BG_SNP} \
--candidate-snp-file ${CAND_SNP} \
--gene-set-file ${GENE_SET} \
--annotation-file ${ANNOTATION} \
--simulations 100000 \
--min-significance 1 \
--gene-definition updownstream5000 \
--threads 8 \
--mode gene \
--min-genes 3 \
--output-file ${OUTPUT} 2>&1 | tee -a ${OUTPUT}_log.txt
