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
print(collect_mags_marker)
mags_dir = collect_mags_marker.split("/collect_mags.done")[0]
mag_metadata_summary = args.mag_metadata_summary
mag_rep_database = args.mag_rep_database


mag_info = pd.read_csv(mag_metadata_summary, sep="\t")
ref_mags = mag_info[mag_info['Representative'].astype(str) == '1']['ID']

print(f"Writing MAG representative database with {len(ref_mags)} MAGs")

db_records = []
with open(mag_rep_database, "a") as db_fh:
    for mag_name in ref_mags:
        fasta_file = os.path.join(os.getcwd(), mags_dir, mag_name+".fa")
        for record in SeqIO.parse(fasta_file, "fasta"):
            db_records.append(record)

with open(mag_rep_database, "w+") as db_fh:
    SeqIO.write(db_records, db_fh, "fasta")

print(f"Done!")
