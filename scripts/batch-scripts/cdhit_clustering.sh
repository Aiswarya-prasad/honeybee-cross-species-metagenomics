#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name cdhit-clustering
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 512G
#SBATCH --time 70:00:00 
#SBATCH --error /...<project_dir_path>.../20230313_apis_species_comparison/results/08_gene_content/cdhit_clustering.err
#SBATCH --output /...<project_dir_path>.../20230313_apis_species_comparison/results/08_gene_content/cdhit_clustering.out

source ~/.bashrc
conda activate 20230313_genes_env

cd /...<project_dir_path>.../20230313_apis_species_comparison

gene_catalog_ffn="results/08_gene_content/gene_catalog_all.ffn"
gene_catalog_faa="results/08_gene_content/gene_catalog_all.faa"
cdhit_genes="results/08_gene_content/gene_catalog_cdhit9590.fasta"
cdhit_genes_faa="results/08_gene_content/gene_catalog_cdhit9590.faa"
headers="results/08_gene_content/cdhit9590.headers"

cat results/06_metagenomicORFs/*/*.ffn > ${gene_catalog_ffn}
cat results/06_metagenomicORFs/*/*.faa > ${gene_catalog_faa}
echo "starting cd-hit-est"
cd-hit-est -i ${gene_catalog_ffn} -o ${cdhit_genes} \
    -c 0.95 -T 64 -M 0 -G 0 -aS 0.9 -g 1 -r 1 -d 0

grep "^>" ${cdhit_genes} | \
cut -f 2 -d ">" | \
cut -f 1 -d " " > ${headers}

seqtk subseq ${gene_catalog_faa} ${headers} > ${cdhit_genes_faa}