#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name popcogent
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 300G
#SBATCH --time 70:00:00 
#SBATCH --error /scratch/aprasad/20230313_apis_species_comparison/results/12_popcogent/logs/popcogent_%j.err
#SBATCH --output /scratch/aprasad/20230313_apis_species_comparison/results/12_popcogent/logs/popcogent_%j.out

source ~/.bashrc
conda activate 20230313_popcogent_env

genus=$1
cd /scratch/aprasad/20230313_apis_species_comparison/scripts/PopCOGenT/src/PopCOGenT
bash PopCOGenT.sh -c /scratch/aprasad/20230313_apis_species_comparison/results/12_popcogent/${genus}/config.sh | tee /scratch/aprasad/20230313_apis_species_comparison/results/12_popcogent/${genus}/popcogent.log

