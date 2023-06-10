#!/bin/bash

cd /scratch/aprasad/211018_Medgenome_india_samples
echo "ID, length" > Figures/Genome_sizes.csv
for file in 07_AnnotationAndPhylogenies/00_genomes/*.fa; do echo "$(echo ${file##*mes/} | sed -e 's/\.fa//'),$(cat $file | grep -v ">" | wc -c)" >> Figures/Genome_sizes.csv; done