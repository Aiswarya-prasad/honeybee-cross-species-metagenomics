#!/bin/bash


outdir=$1
marker=$2
mkdir -p ${outdir}
for dir in results/07_MAG_binng_QC/03_checkm_results/*/bins/*;
do
    bin=$(basename $dir)
    bin_new=${bin/./_}
    sample=${bin_new%%_*}
    cp results/07_MAG_binng_QC/03_checkm_results/${sample}/bins/${bin}/genes.faa \
        ${outdir}/${bin_new}.faa
    cp results/07_MAG_binng_QC/03_checkm_results/${sample}/bins/${bin}/genes.gff \
        ${outdir}/${bin_new}.gff
done
touch ${marker}