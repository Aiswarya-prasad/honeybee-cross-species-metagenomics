#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name parse-cdhit-clustering
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 512G
#SBATCH --time 10:00:00 
#SBATCH --error /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/parse-cdhit-clustering.err
#SBATCH --output /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/parse-cdhit-clustering.out

source ~/.bashrc
conda activate 20230313_scripts_env

cd /scratch/aprasad/20230313_apis_species_comparison

cluster_file="results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta.clstr"
cluster_out="results/08_gene_content/00_cdhit_clustering/clusters_information.tsv"

python3 scripts/parse_cluster_file.py --cluster_file {input.cdhit_clustering} --cluster_out {output.cdhit_clusters}
