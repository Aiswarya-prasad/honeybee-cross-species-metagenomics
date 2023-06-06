#!/usr/bin/env python3

"""
        python scripts/make_mag_rep_database.py \
                --collect_mags_marker {input.collect_mags_marker} \
                --mag_metadata_summary {input.mag_metadata_summary} \
                --mag_rep_database {output.mag_rep_database}
"""

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='Make MAG representative database')
parser.add_argument('--collect_mags_marker', help='collect_mags_marker')
parser.add_argument('--mag_metadata_summary', help='mag_metadata_summary')
parser.add_argument('--mag_rep_database', help='mag_rep_database')
args = parser.parse_args()

collect_mags_marker = args.collect_mags_marker
mag_metadata_summary = args.mag_metadata_summary
mag_rep_database = args.mag_rep_database

collect_mags_marker = "results/09_MAGs_collection/All_mags_sub/MAGs/collect_mags.done",
mags_dir = collect_mags_marker.split("/collect_mags.done")[0]
mag_metadata_summary = "results/09_MAGs_collection/All_mags_sub/MAGs/mag_metadata_summary.tsv"
mag_rep_database = "results/09_MAGs_collection/All_mags_sub/MAGs_rep_db/mag_rep_database.fa"


mag_info = pd.read_csv(mag_metadata_summary, sep="\t")
ref_mags = mag_info[mag_info['reference'].astype(str) == '1']['ID']

if os.path.exists(mag_rep_database):
    os.remove(mag_rep_database)

with open(mag_rep_database, "a") as db_fh:
    for mag_name in ref_mags:
        fasta_file = os.path.join(os.getcwd(), mags_dir, mag_name+".fa")
        for record in SeqIO.parse(fasta_file, "fasta"):
            SeqIO.write(record, db_fh, "fasta")