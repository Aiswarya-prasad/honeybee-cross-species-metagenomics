###
#  For now testing this script interactively with:
#  (base) aprasad@curnagl /scratch/aprasad/211018_Medgenome_india_samples-master$
#  Sinteractive -A pengel_spirit -m 24G -t 2:00:00

#  Sinteractive is running with the following options:

# -A pengel_spirit -c 1 --mem 16G -J interactive -p interactive -t 2:00:00
#
# salloc: Pending job allocation 21926775
# salloc: job 21926775 queued and waiting for resources
# salloc: job 21926775 has been allocated resources
# salloc: Granted job allocation 21926775
# salloc: Waiting for resource configuration
# salloc: Nodes dna064 are ready for job
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

group=g__Gilliamella
project_path=/scratch/aprasad/211018_Medgenome_india_samples
conda_env_midas=/scratch/aprasad/built-envs/36fa9795c7bb520ef68880086b2f369c
conda_env_mapping=/scratch/aprasad/built-envs/d9a6c36f05b352a49e2ad7fed4c52259
# conda_env_mapping=/scratch/aprasad/built-envs/1c5cbf5f292d5200a343c355530053d1 delete that


ref_info=06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
input_og_seq_dir=${project_path}/database/MAGs_database_Orthofinder/${group}/single_ortholog_sequences
mkdir -p ${project_path}/12_species_validation
magOTU_seqs_dir_path=${project_path}/12_species_validation/${group}
mkdir -p ${magOTU_seqs_dir_path}

orf_db=${project_path}/12_species_validation/orfs_db.ffn
if [[ ! -f ${orf_db} ]];
then
  echo 'Creating ORF db which is a concatenated file of all identified genes from each sample'
  cat ${project_path}/12_species_validation/metagenomic_orfs/*/*_orfs.ffn > ${project_path}/12_species_validation/orfs_db.ffn
fi



cd ${project_path}
# conda does not work on the interactive compute node shell

# The following have been reloaded with a version change:
#   1) python/3.8.13 => python/3.9.13

conda activate ${conda_env_midas}

echo 'Performing step1: subset fasta-files in the input directory'

python scripts/subset_magOTU_OG_ffns.py --group ${group} --ref_info ${ref_info} --input_og_seq_dir ${input_og_seq_dir} --magOTU_seqs_dir_path ${magOTU_seqs_dir_path} --log_path ${project_path}/12_species_validation/${group}_messages.log

echo "total OGs for ${group} was $(ls ${input_og_seq_dir}/*.ffn | wc -l)"
for DIR in ${magOTU_seqs_dir_path}/*;
do
  echo "$(ls ${DIR}/*ffn | wc -l) OGs survived for $(basename ${DIR}) of ${group}"
done

echo 'Performing step2: blasting species core sequences against ORF database'

conda activate ${conda_env_mapping}

orf_db_INDEXFILE=$orf_db".nhr"
if [ ! -f "$orf_db_INDEXFILE" ]; then
    echo "Indexing MAG-db file for blasting"
    makeblastdb -in "$orf_db" -dbtype nucl
fi


# for DIR in ${magOTU_seqs_dir_path}/*;
cd ${project_path}
for DIR in  ${magOTU_seqs_dir_path}/163_2  ${magOTU_seqs_dir_path}/165_1  ${magOTU_seqs_dir_path}/166_1  ${magOTU_seqs_dir_path}/167_1;
do
  for magOTU_dir_iter in 12_species_validation/${group}/*; do echo ${magOTU_dir_iter}; ls ${magOTU_dir_iter}/*.ffn | wc -l; ls ${magOTU_dir_iter}/*.blastn | wc -l; done
  echo "Blasting core sequences in directory: "${DIR}
  cd ${DIR}
  COUNTER=0
  for i in $(ls *ffn);
  do
    (( COUNTER++ ))
    name=$(basename $i)
    BLAST_OUTFILE=${DIR}/${name%%.ffn}".blastn"
    blastn -db ${orf_db} -query $i -outfmt 5 -evalue 1e-5 -perc_identity 70 > ${BLAST_OUTFILE}
    if (( "${COUNTER}" % 10 == 0 ));
    then
      echo "Finished blasting ${COUNTER} files.. out of $(ls ${DIR}/*.ffn | wc -l)"
    fi
  done
  echo "Progress on magOTU directories:"
  cd ${project_path}
done


conda activate ${conda_env_midas}

echo 'Performing step3: recruiting ORF sequences from ORF database file, based on blast-files'

python scripts/recruit_orfs.py --group ${group} --ref_info ${ref_info} --orf_db ${orf_db} --magOTU_seqs_dir_path ${magOTU_seqs_dir_path} --log_path ${project_path}/12_species_validation/${group}_messages.log

conda activate ${conda_env_mapping}

echo 'Performing step4: adding recruited ORFs to core sequence gene alignments, and calculating max percentage identity'

python scripts/orf_aln_perc_id.py --group ${group} --ref_info ${ref_info} --orf_db ${orf_db} --magOTU_seqs_dir_path ${magOTU_seqs_dir_path} --log_path ${project_path}/12_species_validation/${group}_messages.log --perc_id ${project_path}/12_species_validation --input_og_seq_dir ${input_og_seq_dir}
