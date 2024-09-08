#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name run_gff_quant
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 100G
#SBATCH --time 10:00:00 
#SBATCH --error /...<project_dir_path>.../20230313_apis_species_comparison/results/08_gene_content/08_gff_quant/logs/08_gff_quant.err
#SBATCH --output /...<project_dir_path>.../20230313_apis_species_comparison/results/08_gene_content/08_gff_quant/logs/08_gff_quant.out

# Doesn't work. Not using.

# Load modules
source ~/.bashrc
conda activate 20230313_scripts_env

# results/08_gene_content/01_profiling_bowtie2/A1-1_filt_genes.bed
# results/08_gene_content/01_profiling_bowtie2/A1-1_filt_genes.gff
# results/08_gene_content/01_profiling_bowtie2/A1-1_mapped_Q20.bam

for file in results/08_gene_content/01_profiling_bowtie2/*_filt_genes.bed
do
    sample=$(basename $file)
    sample=${sample%_filt_genes.bed}
    # sample=${sample%_filt_genes.gff}
    echo $sample
    # gffindex $file
    # gffquant $file results/08_gene_content/01_profiling_bowtie2/${sample}_mapped_Q20.bam -o results/08_gene_content/08_gff_quant/${sample}_gff_quant --mode genome --ambig_mode 1overN
    gffquant --db $file --db_format bed results/08_gene_content/01_profiling_bowtie2/${sample}_mapped_Q20.bam -o results/08_gene_content/08_gff_quant/${sample}_gff_quant --mode small_genome --ambig_mode 1overN
done