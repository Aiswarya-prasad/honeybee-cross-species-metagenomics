#!/usr/bin/env python3
import sys
import os
import argparse
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML


def get_magOTUs_of_group(group_provided, ref_info):
    magOTU_list = []
    with open(ref_info, "r") as ref_info_fh:
        for line in ref_info_fh:
            line = line.strip()
            if line.startswith("ID"):
                header = line+"\n"
                continue
            genome_id = line.split("\t")[0]
            group_read = line.split("\t")[18]
            magOTU = line.split("\t")[11]
            if group_read == group_provided:
                if not magOTU in magOTU_list:
                    magOTU_list.append(magOTU)
            else:
                continue
    return(magOTU_list)

def get_hit_ids(blast_file):
    blastfile_handle = open(blast_file)
    blast_records = NCBIXML.parse(blastfile_handle)
    hit_ids = list()
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            hit_name_full = alignment.title
            hit_name_split = hit_name_full.split()
            hit_id=hit_name_split[-1]
            hit_ids.append(hit_id)
    return(hit_ids)

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--group', metavar="group", required=True, help="Name of genus or group that the script is being run for (historically, phylotype)", action="store")
requiredNamed.add_argument('--ref_info',metavar="ref_info",required=True, help="Info about MAG status with MAG name in first column and other information accroding to checkpoint output", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
requiredNamed.add_argument('--magOTU_seqs_dir_path', metavar="magOTU_seqs_dir_path", required=True, help="Path to store magOTU-wise seperated sequences (inside which magOTU_seqs_dir(s) exist)", action="store")
requiredNamed.add_argument('--orf_db', metavar="orf_db", required=True, help="File containing concatenated filtered metagenomic ORFs", action="store")
requiredNamed.add_argument('--log_path', metavar="log_path", required=True, help="Path to store magOTU-wise seperated sequences (inside which magOTU_seqs_dir(s) exist)", action="store")


args = parser.parse_args()
group = args.group
ref_info = args.ref_info
magOTU_seqs_dir_path = args.magOTU_seqs_dir_path
orf_db = args.orf_db
log_path = args.log_path


if (os.stat(log_path).st_size != 0):
    log_fh = open(log_path, "a")
else:
    log_fh = open(log_path, "w")

magOTU_list = get_magOTUs_of_group(group, ref_info)

magOTU_seqs_dir_dict = {}

for magOTU in magOTU_list:
    print(f"Working on magOTU: {magOTU}")
    magOTU_seqs_dir = os.path.join(magOTU_seqs_dir_path, magOTU)
    magOTU_seqs_dir_dict[magOTU] = magOTU_seqs_dir
    blastn_suffix = "*blastn"
    blastn_files = glob.glob(os.path.join(magOTU_seqs_dir, blastn_suffix))
    print(f"For {group} and magOTU {magOTU}, {len(blastn_files)} OG blastn files were found")
    hit_id_OG = dict()
    OG_hit_id = dict()
    count_progress = 0
    for file in blastn_files:
        count_progress += 1
        if (os.stat(file).st_size != 0):
            OG = os.path.basename(file).split(".blastn")[0]
            hit_ids = get_hit_ids(file)
            OG_hit_id[OG]=dict()
            for id in hit_ids:
                hit_id_OG[id]=OG
                OG_hit_id[OG][id]=1
    if (count_progress % 100 == 0):
        print("Finished parsing",count_progress,"blastn files")
    hit_seq_objects = dict()
    print("Getting seq-records from ORF db")
    for seq_record in SeqIO.parse(orf_db, "fasta"):
        if(seq_record.id in hit_id_OG):
            seq_length = len(seq_record.seq)
            if (seq_length > 200):
                hit_seq_objects[seq_record.id] = seq_record
    print("Printing recruited ORFs to files")
    nb_orfs_recruited=0
    for OG in OG_hit_id.keys():
        ORF_seq_objects = []
        for hit_id in OG_hit_id[OG].keys():
            if (hit_id in hit_seq_objects):
                ORF_seq_objects.append(hit_seq_objects[hit_id])
        nb_objects = len(ORF_seq_objects)
        if (nb_objects != 0):
            orf_outfile = os.path.join(magOTU_seqs_dir, OG + "_orfs.ffn")
            if (os.path.isfile(orf_outfile)):
                 log_fh.write(f"NOTE: The following orf-file already exists: {orf_outfile}\n")
                 print(f"NOTE: The following orf-file already exists: {orf_outfile}")
            else:
                nb_orfs_recruited += len(ORF_seq_objects)
                SeqIO.write(ORF_seq_objects, orf_outfile, "fasta")
    print(f"Recruited {nb_orfs_recruited} ORFs to magOTU {magOTU}")
    log_fh.write(f"Recruited {nb_orfs_recruited} ORFs to magOTU {magOTU}\n")
log_fh.close()
