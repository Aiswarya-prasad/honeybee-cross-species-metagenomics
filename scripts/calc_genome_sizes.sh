#!/bin/bash
bins=$1
outfile=$2
cd /scratch/aprasad/211018_Medgenome_india_samples
echo "ID, length" > Figures/Genome_sizes.csv
for file in ${bins}/*.fa; do echo "$(echo ${file##*mes/} | sed -e 's/\.fa//'),$(cat $file | grep -v ">" | wc -c)" >> ${outfile}; done