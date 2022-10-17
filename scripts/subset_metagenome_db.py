#!/usr/bin/env python3

import os
from Bio import SeqIO
import argparse

def parse_prokka_MAG_header(header):
    """
    read header with prokka prefix and contig number
    and get MAG name eg. MAG_C1.1_12 from gnl|Prokka|MAG_C1.1_12_1
    """
    "_".join(header.split("|")[-1].split("_")[:-1])

def make_db_red(input_fasta, glist):
    output_fasta = input_fasta+"_reduced"
    records = (r for r in SeqIO.parse(input_fasta, "fasta") if parse_prokka_MAG_header(r.id) in glist)
    number_red = SeqIO.write(records, output_fasta, "fasta")
    print(f"adding {number_red} records to reduced database")

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--input', metavar="mag_database", required=True, help="Database of concatenated mags", action="store")
requiredNamed.add_argument('--ref_info',metavar="ref_info",required=True, help="Info about MAG status with MAG name in first column and 1 0r 0 denating if it is ref id in the third (tsv format)", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
args = parser.parse_args()

input_database = args.input
ref_info = args.ref_info

rep_genomes = set()
out_ref = ref_info.split("_ref_info.txt")[0]+"_reduced_ref_info.txt"
with open(ref_info, "r") as ref_info_fh:
    with open(out_ref, "w") as out_ref_fh:
        header = ref_info_fh.readline()
        for line in ref_info_fh:
            line = line.strip()
            genome_id = line.split("\t")[0]
            cluster = line.split("\t")[1]
            rep_genome_status = int(line.split("\t")[2])
            group = line.split("\t")[3]
            if rep_genome_status == 1:
                rep_genomes.add(genome_id)

print(f"making db_red {db_red}")
make_db_red(input_database, rep_genomes)
