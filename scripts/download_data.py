#!/usr/bin/env python

# ============================================================================ #
# download_db.py: download example data for species-validation after cloning from github
#
# Author: Kirsten Ellegaard (kirsten.ellegaard@unil.ch)
#
# ============================================================================ #

import os
import sys
import tempfile
import shutil
import subprocess
import hashlib
import time
import urllib.request
import argparse

# edited by AISWARYA to rename metafile - later check that this works

#Links, file-names, sizes and md5 check-sums
zenodo_links = {
    "genome_db": "https://sandbox.zenodo.org/record/769286/files/genome_db_210402.tar.gz",
    "species_validation": "https://sandbox.zenodo.org/record/769286/files/species_validation.tar.gz",
    "metagenomic_reads": "https://sandbox.zenodo.org/record/769286/files/metagenomic_reads.fastq.tar.gz",
}

file_names = {
    "genome_db": "genome_db_210402.tar.gz",
    "species_validation": "species_validation.tar.gz",
    "metagenomic_reads": "metagenomic_reads.fastq.tar.gz",
}
file_sizes = {
    "genome_db": 370,
    "species_validation": 220,
    "metagenomic_reads": 1100,
}
check_sums = {
    "genome_db": "9fa43b1bfa981115409bada6f52d58e1",
    "species_validation": "b7df6f8ada50125f3ba38e1d13875cfb",
    "metagenomic_reads": "ff640807b464f73c6e9a2a8a6e574f48",
}

#function to print progress bar for python 3
def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r %d%%, %d MB, %d KB/s, %d seconds passed" % (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

#function to download data
def save_f(url, filename):
    if "--no-download-progress" in sys.argv:
        urllib.request.urlretrieve(url, filename)
    else:
        urllib.request.urlretrieve(url, filename, reporthook)

#function to calculate md5
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

#Parse input arguments, check for existing files/directories

parser = argparse.ArgumentParser(description="This script downloads data from zenodo, for running various metagenomic pipelines. Specify at least one data-set from the following options:")
parser.add_argument("--genome_db", help="Download the gut microbiota genomic database", action="store_true")
parser.add_argument("--species_validation", help="Download species validation example data", action="store_true")
parser.add_argument("--metagenomic_reads", help="Download metagenomic reads example data", action="store_true")
args = parser.parse_args()

genome_db_dirs = ['faa_files', 'ffn_files', 'bed_files', 'gff_files', 'Orthofinder']
current_dir = os.getcwd()
if args.genome_db:
    for dir in genome_db_dirs:
        if os.path.isdir(dir):
            print("The following directory already exist in the run-dir:")
            print(dir)
            print("Exiting script!")
            sys.exit(1)
if not any(vars(args).values()):
    parser.print_help(sys.stderr)
    sys.exit(1)

#Download compressed files, check md5 and unpack

all_args = vars(args)
for arg in all_args.keys():
    if (all_args[arg] == False): continue
    #Download compressed file
    print('Downloading gzipped file:', file_names[arg], '(~' + str(file_sizes[arg]) +'Mb)')
    save_f(zenodo_links[arg], file_names[arg])
    #Checking md5
    print('\nChecking MD5..')
    download_md5 = md5(file_names[arg])
    if (download_md5 == check_sums[arg]):
        print('OK!')
    else:
        print('MD5 verification failed for', file_names[arg])
        print('Removing corrupt file, and exiting script')
        os.remove(file_names[arg])
    #Extract files from archive
    print('Extracting files from:', file_names[arg])
    extract_cmd = "tar -zxvf " + file_names[arg] + " -C "+ current_dir
    try:
        FNULL = open(os.devnull, 'w')
        process = subprocess.Popen(extract_cmd.split(),stderr=FNULL,stdout=FNULL)
        output, error = process.communicate()
    except:
        print("Error: failed to extract files\n")
        exit()
    if process.returncode:
        print("Error: failed to extract files\n")
        exit()
    os.remove(file_names[arg])
    os.rename('genome_db_metafile_210402.txt', 'genome_db_210402_metafile.txt')

print('All done!')
