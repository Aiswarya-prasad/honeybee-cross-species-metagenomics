#!/usr/bin/python3

import os
import sys
import argparse
import gzip
import pandas as pd
import numpy as np
from Bio import SeqIO

# def parse_args():
#     parser = argparse.ArgumentParser(description='Subset GTDB-Tk MSA to a list of genomes')
#     parser.add_argument('-i', '--input', type=str, required=True,
#                         help='Input GTDB-Tk MSA file')
#     parser.add_argument('-o', '--output', type=str, required=True,
#                         help='Output GTDB-Tk MSA file')
#     return parser.parse_args()
# for now writing for note_book style execution

metadata = "results/09_MAGs_collection/All_mags_sub/All_mags_sub_metadata_summary.tsv"
gtdb_ref_meta = "config/gtdb_metadata.tsv"
gtdb_msa = "results/09_MAGs_collection/All_mags_sub/gtdb_output/align/All_mags_sub.bac120.msa.fasta.gz"

with open(metadata, "r") as f:
    mag_metadata = pd.read_csv(f, sep="\t", header=0)

# read in gtdb metadata
with open(gtdb_ref_meta, "r") as f:
    gtdb_metadata = pd.read_csv(f, sep="\t", header=0)

# high and medium quality mags
high_quality = mag_metadata.loc[mag_metadata["Quality"] == "high"]
medium_quality = mag_metadata.loc[mag_metadata["Quality"] == "medium"]
mags_df = pd.concat([high_quality, medium_quality])
mags = mags_df['ID'].tolist()
# only genera of interest (i.e. those with at least two or more high or medium quality mag)
genera_of_interest = list(mags_df['Genus'].value_counts().loc[lambda x: x > 2].index.tolist())
genera_of_interest.pop(genera_of_interest.index("g__"))

genus_refs = {}
for genus in genera_of_interest:
    genus_refs[genus] = gtdb_metadata.loc[gtdb_metadata["Genus"] == genus]["Accession"].tolist()

genus_mags = {}
for genus in genera_of_interest:
    genus_mags[genus] = mags_df.loc[mags_df["Genus"] == genus]["ID"].tolist()

# subset gtdb msa to only include genera of interest
with gzip.open(gtdb_msa, "rt") as f:
    records = SeqIO.parse(f, "fasta")
    mag_recs = []
    ref_recs = []
    for record in records:
        if record.id in mags:
            mag_recs.append(record)
        elif record.id in genus_refs:
            ref_recs.append(record)
        else:
            continue

for genus in genera_of_interest:
    with open(f"results/11_phylogenies/gtdb_trees/{genus}/align/{genus}.msa.fasta", "w") as f:
        records_to_write = []
        for record in mag_recs:
            if record.id in genus_mags[genus]:
                records_to_write.append(record)
        for record in ref_recs:
            if record.id in genus_refs[genus]:
                records_to_write.append(record)
        SeqIO.write(records_to_write, f, "fasta")
    