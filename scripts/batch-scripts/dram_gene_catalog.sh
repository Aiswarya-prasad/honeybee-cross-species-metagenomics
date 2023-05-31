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
#SBATCH --error /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/cdhit_clustering.err
#SBATCH --output /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/cdhit_clustering.out

source ~/.bashrc
# conda activate 20230313_dram_env
source /etc/profile.d/lmodstacks.sh
dcsrsoft use old
export PATH=/dcsrsoft/spack/external/dram/v1.2.4/bin:$PATH
module load gcc/9.3.0 python
module load hmmer mmseqs2 prodigal infernal trnascan-se barrnap
which DRAM.py

cd /scratch/aprasad/20230313_apis_species_comparison

gene_catalog_ffn="results/08_gene_content/gene_catalog_all.ffn"
gene_catalog_faa="results/08_gene_content/gene_catalog_all.faa"
cdhit_genes="results/08_gene_content/gene_catalog_cdhit9590.fasta"

DRAM.py annotate-genes \
        -i <input path to the genes faa, NO REGEX> \
        -o <output path> \
        --threads 20 \
        --custom_fasta_loc ./fasta_1.fa \
        --custom_db_name A_seqs \
        --custom_fasta_loc ./fasta_2.fa \
        --custom_db_name  B_seqs \

DRAM.py distill -i annotation/annotations.tsv -o genome_summaries