#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name dram_annotate
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 512G
#SBATCH --time 70:00:00 
#SBATCH --error /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/dram_annotation.err
#SBATCH --output /scratch/aprasad/20230313_apis_species_comparison/results/08_gene_content/dram_annotation.out

source ~/.bashrc
# conda activate 20230313_dram_env
source /etc/profile.d/lmodstacks.sh
dcsrsoft use old
export PATH=/dcsrsoft/spack/external/dram/v1.2.4/bin:$PATH
module load gcc/9.3.0 python
module load hmmer mmseqs2 prodigal infernal trnascan-se barrnap
which DRAM.py

cd /scratch/aprasad/20230313_apis_species_comparison

gene_catalog_faa="results/08_gene_content/gene_catalog_all.faa"
cdhit_genes="results/08_gene_content/gene_catalog_cdhit9590.fasta"
outdir="results/08_gene_content/DRAM_annotation"
outdir_distill="results/08_gene_content/DRAM_distill"

DRAM.py annotate_genes \
        -i ${gene_catalog_faa} \
        -o ${outdir} \
        --threads 20 \
        --verbose

DRAM.py distill -i ${outdir}/annotations.tsv -o ${outdir_distill}