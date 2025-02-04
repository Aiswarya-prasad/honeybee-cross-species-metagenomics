#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name bifidobacterium_orthofinder
#SBATCH --account pengel_spirit
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 100G
#SBATCH --time 20:00:00
#SBATCH --output additional_analyses/results/add_isolates/g__Bifidobacterium/orthofinder_add_isolates.out
#SBATCH --error additional_analyses/results/add_isolates/g__Bifidobacterium/orthofinder_add_isolates.err


cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison
source ~/.bashrc
conda activate /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs/67f21e12bc08b775581c5bda1e986b1a_
orthofinder -t 16 -b /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/additional_analyses/results/add_isolates/g__Bifidobacterium/Results_g__Bifidobacterium/ -f additional_analyses/results/add_isolates/g__Bifidobacterium/isolates_to_add/ | tee additional_analyses/results/add_isolates/g__Bifidobacterium/orthofinder_add_isolates.log