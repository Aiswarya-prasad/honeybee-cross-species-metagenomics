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

def make_db_red(input_fasta, output_fasta):
    records = r for r in SeqIO.parse(input_fasta, "fasta")
    ## INCOMPLETE
    print(f"adding {number_red} records to reduced database")

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--input', metavar="input", required=True, help="Database of concatenated mags", action="store")
requiredNamed.add_argument('--output',metavar="output",required=True, help="Info about MAG status with MAG name in first column and 1 0r 0 denating if it is ref id in the third (tsv format)", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
args = parser.parse_args()

input_database = args.input
output_database = args.output

print(f"making db_red {db_red}")
concat_contigs(input_database, output_database)
