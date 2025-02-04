#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name minpath
#SBATCH --cpus-per-task 5
#SBATCH --mem 50G
#SBATCH --time 20:00:00 
#SBATCH --error additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/minpath.err
#SBATCH --output additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/minpath.out

# Load modules
source ~/.bashrc
conda activate 20230313_scripts_env

python3 additional_analyses/scripts/run_minpath.py
