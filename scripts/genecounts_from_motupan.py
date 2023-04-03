#!/usr/bin/env python3
import sys

"""
This script parses the motupan output and generates from it
a table of genecounts that can be used by the visualization
scripts to plot a heatmap
The aim of this visualization is to see how well the inference
of motupan matches with the occurance (and prevalence) of the gene
How often does the earlier cut-off of half-core in MAGs falsely
choses a gene as "core" and what does motupan make of it
format required:
First column, OG (name of orthogroup)
each subsequent column has genome name and the number of 
    genes of that OG in that genome
Final column is called Total - sum of all the number of genes
    in that OG
Usage:
python3 genecounts_from_motupan.py --input /scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/mOTUpan.tsv --output /scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/mOTUpan_gene_counts.tsv
"""

parser  = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--input', metavar="input", required=True, help="Output of mOTUpan in tsv format", action="store")
requiredNamed.add_argument('--output', metavar='output', required=True, help="File to write output genecounts to", action="store")
args = parser.parse_args()

input_orthofile = args.input
output_orthofile = args.output

# resume here