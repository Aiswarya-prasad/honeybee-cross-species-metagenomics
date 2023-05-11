#!/usr/bin/env python3
import sys
import os.path
import argparse
from Bio import SeqIO

"""
This script renames scaffolds in a fasta file made by spades.
It takes the existing scaffold name (inlcuding node number, coverage and length put in by spades)
and appends it to the name of the sample separated by an underscore.
These are already filtered scaffolds done by the script parse_spades_metagenome.py
The script filt_orfs.py also does a similar renaming in the resulting ffn file.
"""

#Usage: python3 scripts/rename_scaffolds.py --scaffolds_in data/scaffolds/contigs.fasta --scaffolds_out data/scaffolds/contigs_renamed.fasta --sample sample

#Parse input options
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--scaffolds_in', metavar='scaffolds_in', required=True, help="scaffolds from spades", action="store")
requiredNamed.add_argument('--scaffolds_out', metavar="scaffolds_out", required=True, help="reanamed scaffolds with sample name added", action="store")
requiredNamed.add_argument('--sample', metavar="sample", required=True, help="sample name to add to scaffold names", action="store")
args = parser.parse_args()

scaffolds_in = args.scaffolds_in
scaffolds_out = args.scaffolds_out
sample = args.sample

seq_records_filt = []

for seq_record in SeqIO.parse(scaffolds_in, "fasta"):
    header = seq_record.id
    new_header = sample+"_"+header
    seq_record.id = new_header
    seq_record.description = ''
    seq_records_filt.append(seq_record)

SeqIO.write(seq_records_filt, scaffolds_out, "fasta")