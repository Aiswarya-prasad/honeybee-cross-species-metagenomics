#!/usr/bin/env python3
import sys
import os.path
from Bio import SeqIO
import argparse

#Usage: python3 filter_orthologs.py 1_fper_gilli_single_ortho.txt

def get_min_seq_length(ffn_file):
    small_seq_count = 0
    for seq_record in SeqIO.parse(ffn_file, "fasta"):
        seq_length = len(seq_record.seq)
        if seq_length < 300:
            small_seq_count += 1
    return(small_seq_count)

parser = argparse.ArgumentParser()
parser.add_argument('--single_ortho', action='store', help='*_ single_ortho.txt as input')
parser.add_argument('--perc_id', action='store', help='name of file with percentage similarities')
parser.add_argument('--extracted_ffndir', action='store', help='name of directory where annotation output directories are')
parser.add_argument('--ortho_filt', action='store', help='name of file to which outputs are writted')
args = parser.parse_args()

single_ortho = args.single_ortho
perc_id = args.perc_id
extracted_ffndir = args.extracted_ffndir
ortho_filt = args.ortho_filt
ortho_filt_perc_id = ortho_filt.split(".txt")[0]+"_perc_id"+".txt"

#Read the orthofinder-file, save genes as lists in dict
OG_fams = dict()
OG_fam_list = list()
try:
    fh_orthofile = open(single_ortho)
    # fh_orthofile = open(snakemake.input.single_ortho)
except:
    print("Please provide an orthofinder-file with single-copy core gene orthologs as input file for the script")
    print("Exiting script!")
    exit()
for line in fh_orthofile:
    line = line.strip()
    split_line = line.split(':')
    OG_id = split_line.pop(0)
    OG_fams[OG_id] = split_line
    OG_fam_list.append(OG_id)
fh_orthofile.close()

#If the file "perc_id.txt" is present in the run-dir, read and flag OG-ids for which the max perc-id is > 95%
filt_perc_id = dict()
with open(perc_id) as fh_perc_id:
    for line in fh_perc_id:
        line = line.strip()
        split_line = line.split("\t")
        OG_id = split_line[0]
        max_perc_id = split_line[3]
        if (float(max_perc_id) > 0.95):
            filt_perc_id[OG_id] = max_perc_id

#Filter the ortholog file by seq-length (min 300bp) and max inter-SDP perc-id is applicable  (95%)
filt_outfile = ortho_filt
filt_outfile_perc_id = ortho_filt_perc_id
with open(filt_outfile, "w") as out_fh:
    with open(filt_outfile_perc_id, "w") as out_fh_perc_id:
        for OG in OG_fam_list:
            ffn_file = os.path.join(extracted_ffndir, OG + ".ffn")
            nb_small_seq = get_min_seq_length(ffn_file)
            if nb_small_seq == 0:
                OG_genes = " ".join(OG_fams[OG])
                out_fh.write(f"{OG}:{OG_genes}\n")
                if OG in filt_perc_id:
                    continue
                else:
                    out_fh_perc_id.write(f"{OG}:{OG_genes}\n")
