#!/usr/bin/env python3

import os
from Bio import SeqIO
import argparse

def parse_prokka_get_MAG_name(header):
    """
    read header with prokka prefix and contig number
    and get MAG name eg. MAG_C1.1_12 from gnl|Prokka|MAG_C1.1_12_1
    it removes the final number after the underscore which goes from 1 to n
    where n is the number of contigs in that MAG numbered by prokka
    """
    if "|" in header:
        parsed_header = "_".join(header.split("|")[-1].split("_")[:-1])
    else:
        parsed_header = header
    return(parsed_header)

def make_db_red(input_fasta, glist, output_fasta):
    for seq in SeqIO.parse(input_fasta, "fasta"):
        print(seq)
        print(parse_prokka_get_MAG_name(seq.id))
        break
    records = (r for r in SeqIO.parse(input_fasta, "fasta") if parse_prokka_get_MAG_name(r.id) in glist)
    number_red = SeqIO.write(records, output_fasta, "fasta")
    print(f"adding {number_red} records to reduced database")

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--input', metavar="mag_database", required=True, help="Database of concatenated mags", action="store")
requiredNamed.add_argument('--ref_info',metavar="ref_info",required=True, help="Info about MAG status with MAG name in first column and 1 0r 0 denating if it is ref id in the third (tsv format)", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
requiredNamed.add_argument('--output',metavar="ref_info",required=True, help="Info about MAG status with MAG name in first column and 1 0r 0 denating if it is ref id in the third (tsv format)", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
args = parser.parse_args()

input_database = args.input
ref_info = args.ref_info
output_database = args.output

rep_genomes = set()
# ref_info = "06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv"
with open(ref_info, "r") as ref_info_fh:
    for line in ref_info_fh:
        line = line.strip()
        if line.startswith("ID"):
            header = line+"\n"
            continue
        genome_id = line.split("\t")[0]
        cluster = line.split("\t")[11]
        group = line.split("\t")[18]
        rep_genome_status = int(line.split("\t")[19])
        if rep_genome_status == 1:
            rep_genomes.add(genome_id)

print(f"making db_red {input_database}")
make_db_red(input_database, rep_genomes, output_database)
