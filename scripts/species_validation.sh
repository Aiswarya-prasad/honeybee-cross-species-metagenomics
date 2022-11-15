###
#  For now testing this script interactively with:
#  (base) aprasad@curnagl /scratch/aprasad/211018_Medgenome_india_samples-master$
#  Sinteractive -A pengel_spirit -m 16G

#  Sinteractive is running with the following options:

# -A pengel_spirit -c 1 --mem 16G -J interactive -p interactive -t 1:00:00
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

mkdir -p logs/species_validation_test/

group=g__Gilliamella
project_path=/scratch/aprasad/211018_Medgenome_india_samples
conda_env_midas=/scratch/aprasad/built-envs/36fa9795c7bb520ef68880086b2f369c
conda_env_mapping=/scratch/aprasad/built-envs/17ad0633c9d4b928ee617c30e2bae4b5


ref_info=06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
input_seq_dir=${project_path}/database/MAGs_database_Orthofinder/${group}/single_ortholog_sequences
mkdir -p ${project_path}/12_species_validation
magOTU_seqs_dir_path=${project_path}/12_species_validation/${group}
mkdir -p ${magOTU_seqs_dir_path}

echo 'Creating ORF db which is a concatenated file of all identified genes from each sample'
cat ${project_path}/12_species_validation/metagenomic_orfs/*_orfs.ffn > ${project_path}/12_species_validation/orfs_db.ffn
orf_db=${project_path}/12_species_validation/orfs_db.ffn


cd ${project_path}
# conda does not work on the interactive compute node shell
conda activate ${conda_env_mapping}
# The following have been reloaded with a version change:
#   1) python/3.8.13 => python/3.9.13
orf_db_INDEXFILE=$orf_db".nhr"
if [ ! -f "$orf_db_INDEXFILE" ]; then
    echo "Indexing MAG-db file for blasting"
    makeblastdb -in "$orf_db" -dbtype nucl
fi

conda activate ${conda_env_midas}

echo 'Performing step1: subset fasta-files in the input directory'

python scripts/subset_magOTU_OG_ffns.py --group ${group} --ref_info ${ref_info} --input_seq_dir ${input_seq_dir} --magOTU_seqs_dir_path ${magOTU_seqs_dir_path} --log_path ${project_path}/12_species_validation/${group}_messages.log

echo "total OGs for ${group} was $(ls $input_seq_dir/*.ffn | wc -l)"
for DIR in ${magOTU_seqs_dir_path}/*;
do
  echo "$(ls ${DIR}/*ffn | wc -l) OGs survived for $(basename ${DIR}) of ${group}"
done

echo 'Performing step2: blasting species core sequences against ORF database'

for DIR in ${magOTU_seqs_dir_path}/*;
do
  echo "Blasting core sequences in directory: "${DIR}
  cd "${DIR}"
  COUNTER=0
  for i in $(ls *ffn);
  do
    (( COUNTER++ ))
    name=$(basename $i)
    BLAST_OUTFILE=${DIR}/${name%%.ffn}".blastn"
    blastn -db ${orf_db} -query $i -outfmt 5 -evalue 1e-5 -perc_identity 70 > ${BLAST_OUTFILE}
    if (( "${COUNTER}" % 10 == 0 ));
    then
      echo "Finished blasting ${COUNTER} files.."
    fi
  done
done

cd ${project_path}

echo 'Performing step3: recruiting ORF sequences from ORF database file, based on blast-files'
