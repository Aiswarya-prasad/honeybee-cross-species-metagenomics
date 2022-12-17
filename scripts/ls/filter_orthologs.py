#!/usr/bin/env python3
import sys
import os.path
from Bio import SeqIO

#Usage: python3 filter_orthologs.py 1_fper_gilli_single_ortho.txt

def get_min_seq_length(ffn_file):
    small_seq_count = 0
    for seq_record in SeqIO.parse(ffn_file, "fasta"):
        seq_length = len(seq_record.seq)
        if seq_length < 300:
            small_seq_count += 1
    return(small_seq_count)

#Read the orthofinder-file, save genes as lists in dict
OG_fams = dict()
OG_fam_list = list()
try:
    fh_orthofile = open(snakemake.input.ortho_single)
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
with open(snakemake.input.perc_id) as fh_perc_id:
    for line in fh_perc_id:
        line = line.strip()
        split_line = line.split("\t")
        OG_id = split_line[0]
        max_perc_id = split_line[3]
        if (float(max_perc_id) > 0.95):
            filt_perc_id[OG_id] = max_perc_id

#Filter the ortholog file by seq-length (min 300bp) and max inter-SDP perc-id is applicable  (95%)
filt_outfile = snakemake.output.ortho_filt
fh_outfile = open(filt_outfile, 'w')
for OG in OG_fam_list:
    if OG in filt_perc_id:  continue
    ffn_file = snakemake.params.ffn_dir + OG + ".ffn"
    nb_small_seq = get_min_seq_length(ffn_file)
    if nb_small_seq == 0:
        OG_genes = " ".join(OG_fams[OG])
        fh_outfile.write(f"{OG}:{OG_genes}\n")
fh_outfile.close()
