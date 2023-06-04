#!/bin/bash
# script to back up to NAS (via KE's workstation using rsync)

: '
This script can be used to regularly backup a given source directory
to a destination directory either on the same system or via the network
(using rsync) but it needs to be configured such that password is not
expected each time.

the first argument is the source directory (everything inside it will be
transferred - remove or add patterns to --exculde if needed)

Make sure that the destination path exists!

usage: ./backup.sh /scratch/<> /nas/<> <path to log file>

Example: 
bash scripts/backup.sh /scratch/aprasad/20230313_apis_species_comparison/ ~/nas_recherche/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/20230313_apis_species_comparison/ /scratch/aprasad/backup_20230313_apis_species_comparison.log
'

date_str=$(date +'%y%m%d')
time_str=$(date +'%R')

source=$1
destination=$2
# Write log file to current directory
logfile=$3
touch ${logfile}

echo "#######################"
echo "backup being done on" | tee -a ${logfile}
echo ${date_str} | tee -a ${logfile}
echo "at" | tee -a ${logfile}
echo ${time_str} | tee -a ${logfile}
echo "#######################" | tee -a ${logfile}
echo "#######################" | tee -a ${logfile}
echo "#######################" | tee -a ${logfile}
echo "#############################################################" | tee -a ${logfile}
echo "starting backup of" | tee -a ${logfile}
echo $source | tee -a ${logfile}
echo "#############################################################" | tee -a ${logfile}

# do not add the --update option (files that are touched in the remote server will be considered newer and not be updated
#                                 also, files that are copied fresh will be considered "newer" rather than "uptodate")
# rsync -i --delete --backup --backup-dir="BACKUP/" --exclude "*.git*" --exclude "*_ortho_sequences*" --exclude "BACKUP/" --exclude "Report_cache/" --exclude "00_rawdata/" --exclude "00_Rawdata/" --exclude "00_RawData/" --exclude ".snakemake/" -avP ${source}/ ${destination}/ | tee -a ${logfile}

# rsync -avPni --delete --backup --backup-dir="BACKUP/" --exclude "BACKUP/" --exclude ".snakemake/" ${source}/ ${destination}/ | tee -a ${logfile}
# rsync -avPi --delete --backup --backup-dir="BACKUP/" --exclude "BACKUP/" --exclude ".snakemake/" ${source}/ ${destination}/ | tee -a ${logfile}
rsync -avPi --delete --backup --backup-dir="BACKUP/" --exclude "BACKUP/" --exclude ".snakemake/" --exclude "*.bam" --exclude 00_*/ --exclude 01_c*/ --exclude 01_trimmedconcatreads --exclude 01_trimmedconcatreads ${source}/ ${destination}/ | tee -a ${logfile}

echo "#############################################################" | tee -a ${logfile}
echo "backup done" | tee -a ${logfile}
echo "#############################################################" | tee -a ${logfile}
