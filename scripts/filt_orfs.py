#!/usr/bin/env python3
import os
import sys
import re
import glob
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

"""
This script reads the ffn file made by prodigal and makes a new one where the name of the contig
starts with sample name prefixed with '_' as a separator and only keeps orfs that are marked have 
a flag of 00 and are more than 300 bp long.
The corresponding faa records can be parsed out of the file in 06_metagenomicORFs/{sample}/{sample}.faa
The corresponding gff records can be parsed out of the file in 06_metagenomicORFs/{sample}/{sample}.faa
or can be written by parsing the ffn file using biopython (https://biopython.org/wiki/GFF_Parsing).
In addition, filtering steps are now added to exclude the genes that sit on contigs that are 
annotated by whokaryote as belonging to a eukaryote
# taxonomy based filtering can be done later after the gene catalog is made

usage example:
    python scripts/filt_orfs.py --who {whokaryote_result} 
            --tax {tax_result}  --taxtype {params.taxtype} \
            --ffn_in {output.scaffolds_ffn} --ffn_out {output.orfs} \
            --sample {wildcards.sample} --log {output.filt_log}
"""

# Usage: python3 scripts/filt_orfs.py --ffn_in data/ffn/contigs.ffn --ffn_out data/ffn/contigs_filt.ffn --sample sample --log data/ffn/contigs_filt.log

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--ffn_in', metavar="ffn_in", required=True, help="provided ffn input", action="store")
requiredNamed.add_argument('--ffn_out',metavar="ffn_out",required=True, help="filtered ffn output", action="store")
requiredNamed.add_argument('--sample',metavar="ffn_out",required=True, help="sample name", action="store")
requiredNamed.add_argument('--who',metavar="ffn_out",required=True, help="output tsv from whokaryote", action="store")
requiredNamed.add_argument('--tax',metavar="ffn_out",required=True, help="taxonomic profiling of contigs", action="store")
requiredNamed.add_argument('--taxtype',metavar="ffn_out",required=True, help="taxonomic profiling of contigs - name of tool / type of input", action="store")
requiredNamed.add_argument('--log',metavar="ffn_out",required=True, help="filtered ffn output", action="store")

args = parser.parse_args()
ffn_in = args.ffn_in
ffn_out = args.ffn_out
sample = args.sample
who = args.who
tax = args.tax
taxtype = args.taxtype
log = args.log

# for testing only
# who="results/05_assembly/contig_fates/whokaryote/D2-4/whokaryote_predictions_S.tsv"
# tax="results/05_assembly/contig_fates/kaiju/nr/D2-4_fullnames.txt"
# taxtype="kaiju"
# ffn_in="results/06_metagenomicORFs/D2-4/prodigal_out/D2-4.ffn"
# ffn_out="results/06_metagenomicORFs/D2-4/filt_orfs/D2-4.ffn"
# sample="D2-4"
# log="results/06_metagenomicORFs/D2-4/orfs_filt_sumary.log"
# 

count_records = 0
lost_records_length_partial = 0
lost_records_eukaryote = 0
count_records_filt = 0
unknown_who_and_tax = 0

seq_records_filt = []

contigs_karyote = {}

with open(who, "r") as who_fh:
    for line in who_fh:
        line = line.strip()
        if line.startswith("NODE"):
            contig = line.split("\t")[0]
            whokaryote = line.split("\t")[1]
            contigs_karyote[contig] = whokaryote

tax_df = pd.DataFrame(columns=["classification", "contig", "taxonomy"])

with open(tax, "r") as tax_fh:
    rows_to_concat = []
    for line in tax_fh:
        line = line.strip()
        data = line.split("\t")
        classification = data[0]
        contig = data[1]
        taxonomy = ''
        if classification == 'C':
            try:
                taxonomy = data[7]
            except IndexError:
                if data[2] == '1':
                    taxonomy = 'root'
        elif classification == 'U':
            taxonomy = 'unknown'

        # Append the data as a dictionary to the list
        rows_to_concat.append({"classification": classification, "contig": contig, "taxonomy": taxonomy})

# Concatenate the list of dictionaries to create the DataFrame
tax_df = pd.concat([tax_df, pd.DataFrame(rows_to_concat)], ignore_index=True)

# contig_name = "NODE_169933_length_1001_cov_4.932347"
def tax_info_ok(contig_name):
    if taxtype == "kaiju":
        if tax_df.loc[tax_df["contig"] == contig_name]["classification"].values == 'U':
            return False
        if tax_df.loc[tax_df["contig"] == contig_name]["classification"].values == 'C':
            return True
        print(f"could not determine tax info for contig {contig_name}")
    else:
        print(f"taxtype {taxtype} not supported yet. No filtering for this step.")

def isprokaryote(contig_name):
    global unknown_who_and_tax
    if contig_name not in contigs_karyote:
        print(f"contig {contig_name} not found in whokaryote results")
        print(f"checking tax info for contig {contig_name}")
        try:
            return tax_info_ok(contig_name)
        except:
            unknown_who_and_tax += 1
            return False
    else:
        if contigs_karyote[contig_name] == "prokaryote":
            return True
        else:
            return False

for seq_record in SeqIO.parse(ffn_in, "fasta"):
    count_records += 1
    header = seq_record.id
    length = len(seq_record)
    partial_flag = seq_record.description.split(";")[1].split("=")[1]
    contig_name = "_".join(header.split("_")[1:-1])
    if length > 300 and partial_flag == "00":
        if isprokaryote(contig_name):
            count_records_filt += 1
            seq_record.id = header
            seq_record.description = ''
            seq_records_filt.append(seq_record)
        else:
            lost_records_eukaryote += 1
    else:
        lost_records_length_partial += 1

SeqIO.write(seq_records_filt, ffn_out, "fasta")

with open(log, "w") as log_fh:
    log_fh.write(f"sample\tcount_records\trecords_filt_length\trecords_filt_eukaryote\trecords_filt\tfraction\n")
    print(f"sample\tcount_records\trecords_filt_length\trecords_filt_eukaryote\trecords_filt\tfraction\n")
    records_filt_length = count_records - lost_records_length_partial
    # this adds back the ones that were removed because they were unknown
    # even though the unknown ones were removed anyway
    # this decrease is reflected in count_records_filt anyway
    records_filt_eukaryote = records_filt_length - lost_records_eukaryote + unknown_who_and_tax
    fraction = count_records_filt/count_records*100
    log_fh.write(f"{sample}\t{count_records}\t{records_filt_length}\t{records_filt_eukaryote}\t{count_records_filt}\t{fraction}\n")
    print(f"{sample}\t{count_records}\t{records_filt_length}\t{records_filt_eukaryote}\t{count_records_filt}\t{fraction}\n")
