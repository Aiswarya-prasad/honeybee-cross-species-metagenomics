#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Collect prodigal output based on checkm output
gff file ID will be renamed to sample_ID eg. 1_1 -> D4-1_1_1
faa file ID will be renamed to sample_ID eg. 1_1 -> D4-1_1_1
fna files is the same as the input mag file and can be accessed from there
ffn file is made by getting the sequence from the fa file based on the gff file
Usage:
python scripts/rename_collect_prodigal_mags.py \
            --mag {wildcards.mag} \
            --mag_fa {input.mag_fa} \
            --prodigal_faa {input.prodigal_checkm_faa} \
            --prodigal_gff {input.prodigal_checkm_gff} \
            --outdir {params.outdir} \
            --outfile_fna {output.renamed_fna} \
            --outfile_ffn {output.renamed_ffn} \
            --outfile_faa {output.renamed_faa} \
            --outfile_gff {output.renamed_gff}
"""

parser = argparse.ArgumentParser(description='Rename MAGs')
parser.add_argument('--mag', type=str, help='MAG name')
parser.add_argument('--mag_fa', type=str, help='MAG fasta file')
parser.add_argument('--prodigal_faa', type=str, help='Prodigal faa file')
parser.add_argument('--prodigal_gff', type=str, help='Prodigal gff file')
parser.add_argument('--outdir', type=str, help='Output directory')
parser.add_argument('--outfile_ffn', type=str, help='Output ffn file')
parser.add_argument('--outfile_faa', type=str, help='Output faa file')
parser.add_argument('--outfile_gff', type=str, help='Output gff file')
args = parser.parse_args()

mag = args.mag
mag_fa = args.mag_fa
prodigal_faa = args.prodigal_faa
prodigal_gff = args.prodigal_gff
outdir = args.outdir
outfile_ffn = args.outfile_ffn
outfile_faa = args.outfile_faa
outfile_gff = args.outfile_gff

mag = "F3-1_1"
mag_fa = "results/09_MAGs_collection/All_mags_sub/MAGs/F3-1_1.fa"
prodigal_faa = "results/09_MAGs_collection/All_mags_sub/prodigal_output/from_checkm/F3-1_1.faa"
prodigal_gff = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed/F3-1_1.gff"
outdir = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed"
outfile_ffn = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed/F3-1_1.ffn" 
outfile_faa = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed/F3-1_1.faa"
outfile_gff = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed/F3-1_1.gff"

