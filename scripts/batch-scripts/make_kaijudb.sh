#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name kaiju_db_index
#SBATCH --cpus-per-task 5
#SBATCH --mem 400G
#SBATCH --time 70:00:00 
#SBATCH --error /...<kaiju_dir_path>.../kaiju_db/kaiju_db_nr_makedb_index.err
#SBATCH --output /...<kaiju_dir_path>.../kaiju_db/kaiju_db_nr_makedb_index.out

source ~/.bashrc
conda activate kaiju

# # cd /...<kaiju_dir_path>.../kaiju_db/
# # mkdir nr
# # cd /...<kaiju_dir_path>.../kaiju_db/nr
# # kaiju-makedb -s nr

# cd /...<kaiju_dir_path>.../kaiju_db/
# mkdir nr_euk
# cd /...<kaiju_dir_path>.../kaiju_db/nr_euk
# kaiju-makedb -s nr_euk

cd /...<kaiju_dir_path>.../kaiju_db/nr
kaiju-makedb --index-only -s nr

# kaiju-makedb --index-only -s nr_euk