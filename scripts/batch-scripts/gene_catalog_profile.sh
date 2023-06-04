#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name profile_genes_M1-1
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 30G
#SBATCH --time 04:00:00 
#SBATCH --error /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/01_profiling/M1-1_profile_genes.log.cluster.e
#SBATCH --output /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/01_profiling/M1-1_profile_genes.log.cluster.o

source ~/.bashrc
conda activate 20230313_mapping_env

cd /scratch/aprasad/20230313_apis_species_comparison

sample="M1-1"
reads1="results/01_cleanreads/${sample}_R1_repaired.fastq.gz"
reads2="results/01_cleanreads/${sample}_R2_repaired.fastq.gz"
gene_catalog="results/08_gene_content/gene_catalog_cdhit9590.fasta"

if [ ! -f ${gene_catalog}.pac ]
then
    bwa index ${gene_catalog}
else
    echo "index exists"
fi

bam="results/08_gene_content/01_profiling/${sample}_mapped.bam"
flagstat="results/08_gene_content/01_profiling/${sample}_mapped.flagstat"
depth="results/08_gene_content/01_profiling/${sample}_mapped.depth"
filter_script="scripts/filter_bam.py"

bwa mem -a -t 4 ${gene_catalog} ${reads1} ${reads2} \
        | samtools view -F 4 -h - |  python3 ${filter_script} -e 5 -m 50 | \
        samtools sort -O bam -@ 4 > ${bam}

samtools flagstat -@ 4 ${bam} > ${flagstat}

samtools depth -s -a ${bam} > ${depth}