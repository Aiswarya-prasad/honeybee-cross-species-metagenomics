#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name metabat2-try
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 50G
#SBATCH --time 02:00:00 
#SBATCH --error /scratch/aprasad/20230313_apis_species_comparison_backmapping/results/07_MAG_binng_QC/metabat2_try.err
#SBATCH --output /scratch/aprasad/20230313_apis_species_comparison_backmapping/results/07_MAG_binng_QC/metabat2_try.out

source ~/.bashrc
conda activate /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs/4e478f9b8402bd147ebb1ddae7028160_

metabat2 -i results/07_MAG_binng_QC/00_assembled_scaffolds/C1-1/C1-1_scaffolds.fasta \
    -a results/07_MAG_binng_QC/01_backmapping/C1-1_depths/compare_conditions/SampleSize/sub_40/C1-1_depthfile.txt \
    -o results/07_MAG_binng_QC/02_bins/SampleSize/40/C1-1 \
    --minContig 2000 --maxEdges 500 -x 1 --minClsSize 200000 \
    --onlyLabel --saveCls --noBinOut -v