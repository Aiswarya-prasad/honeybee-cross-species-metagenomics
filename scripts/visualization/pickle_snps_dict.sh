#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name snps_dict
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 200G
#SBATCH --time 20:00:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/batch_log/snps_dict.e
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/batch_log/snps_dict.e

source ~/.bashrc
conda activate 20230313_scripts_env

python3 scripts/visualization/pickle_snps_dict.py