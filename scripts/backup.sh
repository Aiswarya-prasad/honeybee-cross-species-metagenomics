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

usage: ./backup.sh /scratch/<> /nas/<>
'

date_str=$(date +'%y%m%d')
time_str=$(date +'%R')

source=$1
destination=$2
# Write log file to current directory
touch ${destination}/nas_backup_${date_str}.log

echo "#######################"
echo "backup being done on" | tee -a ${destination}/nas_backup_${date_str}.log
echo ${date_str} | tee -a ${destination}/nas_backup_${date_str}.log
echo "at" | tee -a ${destination}/nas_backup_${date_str}.log
echo ${time_str} | tee -a ${destination}/nas_backup_${date_str}.log
echo "#######################" | tee -a ${destination}/nas_backup_${date_str}.log
echo "#######################" | tee -a ${destination}/nas_backup_${date_str}.log
echo "#######################" | tee -a ${destination}/nas_backup_${date_str}.log
echo "#############################################################" | tee -a ${destination}/nas_backup_${date_str}.log
echo "starting backup of" | tee -a ${destination}/nas_backup_${date_str}.log
echo $source | tee -a ${destination}/nas_backup_${date_str}.log
echo "#############################################################" | tee -a ${destination}/nas_backup_${date_str}.log

# do not add the --update option (files that are touched in the remote server will be considered newer and not be updated
#                                 also, files that are copied fresh will be considered "newer" rather than "uptodate")
rsync -i --delete --backup --backup-dir="BACKUP/" --exclude "*_ortho_sequences*" --exclude "BACKUP/" --exclude "Report_cache/" --exclude "00_rawdata/" --exclude "00_Rawdata/" --exclude "00_RawData/" --exclude ".snakemake/" -avP ${source}/ ${destination}/ | tee -a ${destination}/nas_backup_${date_str}.log

echo "#############################################################" | tee -a ${destination}/nas_backup_${date_str}.log
echo "backup done" | tee -a ${destination}/nas_backup_${date_str}.log
echo "#############################################################" | tee -a ${destination}/nas_backup_${date_str}.log
