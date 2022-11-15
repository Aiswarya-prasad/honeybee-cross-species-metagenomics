#!/usr/bin/env python3
import sys
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

def parse_MAG_name_from_geneid(header):
    """
    read header with gene id in ffn files and get MAG name
    eg. MAG_C1.5_18 from >MAG_C1.5_18_01186 hypothetical protein
    it removes the final number after the underscore which goes from 1 to n
    where n is the gene id numbered by prokka
    """
    parsed_header = "_".join(header.split(" ")[0].split(">")[1].split("_")[:-1])
    return(parsed_header)

def get_magOTU(mag, ref_info):
    with open(ref_info, "r") as ref_info_fh:
        for line in ref_info_fh:
            line = line.strip()
            if line.startswith("ID"):
                header = line+"\n"
                continue
            genome_id = line.split("\t")[0]
            magOTU = line.split("\t")[11]
            if mag == genome_id:
                return(magOTU)

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--ref_info',metavar="ref_info",required=True, help="Info about MAG status with MAG name in first column and other information accroding to checkpoint output", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
requiredNamed.add_argument('--input_seq_dir', metavar="input_seq_dir", required=True, help="Directory containing genes, aligned and back-translated genes", action="store")

args = parser.parse_args()

ref_info = args.ref_info
input_seq_dir = args.input_seq_dir
