#!/usr/bin/env python3
import os
import sys
import argparse
import itertools
import pandas as pd
from collections import Counter

"""
Read the orthogroups genecounts file and create a cog file for motupan
The file should either be a JSON-file encoding a dictionary where 
the keys are the genome names and the values are lists of traits/genes
(example in example_files/example_genome2cog.json). 
Or a TAB-separated file, where the first column is the genome name and 
followed by TAB-separated trait/gene-names (example in example_files/example_genome2cog.tsv).
For this script,
the output will be a tsv file with genome name as the first column followed
by a tab-seperated list of OG names next to it.

example usage:
python3 prepare_motupan_OG_file.py --input Orthogroups.GeneCount.tsv --output motupan_cog.tsv
"""

parser  = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--input', metavar="input", required=True, help="Orthofinder output listing genes inside each orthogroup", action="store")
requiredNamed.add_argument('--output_core', metavar='output', required=True, help="output file to be given to motupan", action="store")
args = parser.parse_args()

input_orthofile = args.input
output_cogfile = args.output

input_orthofile = "/...<project_dir_path>.../211018_Medgenome_india_samples/database/MAGs_database_Orthofinder/g__Lactobacillus/OrthoFinder/Results_g__Lactobacillus/Orthogroups/Orthogroups.GeneCount.tsv"
output_cogfile = "/...<project_dir_path>.../211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/motupan_cog.tsv" 

ortho_df = pd.read_csv(input_orthofile, sep = "\t")
ortho_df = ortho_df.set_index('Orthogroup')
ortho_df.index.names = [None]
cog_df = {}
for genome in ortho_df.keys():
    if genome == "Total":
        continue
    cog_df[genome] = [x for x in ortho_df.index[ortho_df[genome] > 0]]

with open(output_cogfile, "w") as out_fh:
    for genome in cog_df:
        out_list = "\t".join(cog_df[genome])
        out_fh.write(f"{genome}\t{out_list}\n")