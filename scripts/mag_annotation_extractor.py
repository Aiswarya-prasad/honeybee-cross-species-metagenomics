#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

"""
not using this code at the moment
"""


parser = argparse.ArgumentParser(description="")

# Step 1: Identify MAG contigs
mag_contigs_file = "path/to/mag_contigs.txt"  # File containing MAG contig IDs
with open(mag_contigs_file, 'r') as f:
    mag_contigs = set([line.strip() for line in f])

# Step 2: Parse Prodigal output to collect MAG-specific outputs
prodigal_output_dir = "path/to/prodigal_output"  # Directory containing Prodigal output files
mag_faa_file = "path/to/mag_specific.faa"  # Output MAG-specific FAA file
mag_gff_file = "path/to/mag_specific.gff"  # Output MAG-specific GFF file
mag_ffn_file = "path/to/mag_specific.ffn"  # Output MAG-specific FFN file

with open(mag_faa_file, 'w') as faa_out, open(mag_gff_file, 'w') as gff_out, open(mag_ffn_file, 'w') as ffn_out:
    for file_name in os.listdir(prodigal_output_dir):
        if file_name.endswith(".faa"):
            faa_file = os.path.join(prodigal_output_dir, file_name)
            sample_name = file_name.split(".")[0]
            with open(faa_file, 'r') as faa_in:
                for line in faa_in:
                    if line.startswith(">"):
                        contig_id = line.strip().lstrip('>')
                        if contig_id in mag_contigs:
                            faa_out.write(line)
                            for i in range(1, 4):
                                faa_out.write(faa_in.readline())
                    else:
                        continue

        elif file_name.endswith(".gff"):
            gff_file = os.path.join(prodigal_output_dir, file_name)
            sample_name = file_name.split(".")[0]
            with open(gff_file, 'r') as gff_in:
                for line in gff_in:
                    contig_id = line.split('\t')[0]
                    if contig_id in mag_contigs:
                        gff_out.write(line)

        elif file_name.endswith(".ffn"):
            ffn_file = os.path.join(prodigal_output_dir, file_name)
            sample_name = file_name.split(".")[0]
            with open(ffn_file, 'r') as ffn_in:
                for line in ffn_in:
                    if line.startswith(">"):
                        contig_id = line.strip().lstrip('>')
                        if contig_id in mag_contigs:
                            ffn_out.write(line)
                            for i in range(1, 4):
                                ffn_out.write(ffn_in.readline())
                    else:
                        continue
