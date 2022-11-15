###
#  For now testing this script interactively with:
#  (base) aprasad@curnagl /scratch/aprasad/211018_Medgenome_india_samples-master$
#  Sinteractive

#  Sinteractive is running with the following options:

#  -c 1 --mem 8G -J interactive -p interactive -t 1:00:00

#  salloc: Pending job allocation 21926451
#  salloc: job 21926451 queued and waiting for resources
#  salloc: job 21926451 has been allocated resources
#  salloc: Granted job allocation 21926451
#  salloc: Waiting for resource configuration
#  salloc: Nodes dna064 are ready for job
###

#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name species_validation_test
#SBATCH --account pengel_spirit
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 50G
#SBATCH --time 10:00:00
#SBATCH --output /scratch/aprasad/211018_Medgenome_india_samples/logs/species_validation_test/%x_%j.out
#SBATCH --error /scratch/aprasad/211018_Medgenome_india_samples/logs/species_validation_test/%x_%j.err

####--------------------------------------
##preparation
##set you bash variables in order to quickly call them in the script
####--------------------------------------

mkdir -p logs/species_validation_test/

group=g__Gilliamella
project_path=/scratch/aprasad/211018_Medgenome_india_samples
conda_env_scripts=/scratch/aprasad/built-envs/93bf18c4f86aa15bcc121aa5905f1e20

input_seq_dir=${project_path}/database/MAGs_database_Orthofinder/${group}/single_ortholog_sequences/

cd ${project_path}
conda activate ${conda_env_scripts}
module load gcc blast-plus muscle
