#!/usr/bin/env python3
import os
import sys
import argparse
import itertools
import pandas as pd
from collections import Counter

"""
This script parses the motupan output and generates from it
a table of genecounts that can be used by the visualization
scripts to plot a heatmap
The aim of this visualization is to see how well the inference
of motupan matches with the occurance (and prevalence) of the gene
How often does the earlier cut-off of half-core in MAGs falsely
choses a gene as "core" and what does motupan make of it
format required:
First column, OG (name of orthogroup) - this column can be unnamed
each subsequent column has genome name and the number of 
    genes of that OG in that genome
Usage:
python3 genecounts_from_motupan.py --input /scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/mOTUpan.tsv --output_genecounts /scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/mOTUpan_gene_counts.tsv --output_cores "/scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/g__g__Lactobacillus_single_ortho.txt"
"""

parser  = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--input', metavar="input", required=True, help="Output of mOTUpan in tsv format", action="store")
requiredNamed.add_argument('--output_core', metavar='output', required=True, help="File to write list of single copy core genes to in the same format as is obtained by previously implemented script from Orthofinder", action="store")
requiredNamed.add_argument('--output_genecounts', metavar='output', required=True, help="File to write output genecounts to", action="store")
args = parser.parse_args()

input_orthofile = args.input
output_orthofile = args.output

input_orthofile = "/scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/mOTUpan.tsv"
output_genecounts = "/scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/mOTUpan_gene_counts.tsv"
output_cores = "/scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/g__g__Lactobacillus_single_ortho.txt"

def genome_from_gene_name(gene):
    genome = "_".join(gene.split("_")[:-1])
    return(genome)

OG_df = pd.DataFrame()
header = False
with open(input_orthofile, "r") as in_file:
    for line in in_file:
        if line:
            line = line.strip()
            if line.startswith("#"):
                if line.startswith("#genomes="):
                    line = line.split("#genomes=")[1]
                    genomes_list = [x.split(":")[0] for x in line.split(";")]
                    OG_df = pd.DataFrame(index = genomes_list)
                continue
            if header:
                pass
            else:
                header = line
                continue
            split_line = line.split("\t")
            if len(split_line) < 6:
                continue
            OG = split_line[0]
            freqs = Counter([genome_from_gene_name(x) for x in split_line[6].split(";")])
            OG_df[OG] = pd.Series(freqs)
            OG_df[OG] = OG_df[OG].fillna(0)
        
OG_df.T.to_csv(output_genecounts, sep = "\t")
        
# resume here to write 

# Also make a script to get single copy core genes from motupan annotations!