#!/bin/bash

######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name kaiju_db
#SBATCH --cpus-per-task 5
#SBATCH --mem 200G
#SBATCH --time 70:00:00 
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/kaiju_db/kaiju_db_nr.err
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/kaiju_db/kaiju_db_nr.out

source ~/.bashrc
conda activate kaiju

cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/kaiju_db/
mkdir nr
cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/kaiju_db/nr
kaiju-makedb -s nr