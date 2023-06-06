#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

"""
collect all the gene counts from the depthfiles in the directory provided
and create a matrix of gene counts mxn where n is the number of genes and m is the number of samples
each row is a gene and each column is a sample
"""


parser = argparse.ArgumentParser(description="create a matrix of gene counts from depthfiles")
parser.add_argument("-d", "--depthfiles_dir", help="directory containing depthfiles", required=True)
parser.add_argument("-s", "--suffix", help="pattern appended after sample name to depth file", required=True)
parser.add_argument("-l", "--samples_list", nargs = "+", help="optional argument to specify samples", required=False)
parser.add_argument("-g", "--gene_lengths", help="file to write gene lengths to", required=True)
parser.add_argument("-c", "--gene_catalog", help="gene catalog fasta that all samples were mapped to", required=True)
parser.add_argument("-o", "--outfile", help="output file name", required=True)
args = parser.parse_args()

def make_gene_lengths_file(catalog, lengths_file):
    with open(catalog, "r") as catalog_file:
        with open(lengths_file, "w+") as lengths_fh:
            for record in SeqIO.parse(catalog_file, "fasta"):
                lengths_fh.write(f"{record.id}\t{len(record.seq)}\n")
    



depths_dir = args.depthfiles_dir
suffix = args.suffix
if args.samples:
    samples = args.samples
    depth_files = [depths_dir + sample + suffix for sample in samples]
else:
    depth_files = glob.glob(depths_dir + "*" + suffix)
    samples = [os.path.basename(depth_file).split(suffix)[0] for depth_file in depth_files]
outfile = args.outfile
gene_catalog = args.gene_catalog
gene_lengths = args.gene_lengths


depths_dir = "results/08_gene_content/01_profiling/"
suffix = "_mapped.depth"
depth_files = glob.glob(depths_dir + "*" + suffix)
samples = [os.path.basename(depth_file).split(suffix)[0] for depth_file in depth_files]
outfile = "results/08_gene_content/03_gene_counts/all_genes_matrix.txt"
gene_catalog = "results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta"
gene_lengths = "results/08_gene_content/03_gene_counts/gene_lengths.txt"
os.makedirs(os.path.dirname(gene_lengths), exist_ok=True)

make_gene_lengths_file(gene_catalog, gene_lengths)

gene_lengths_dict = {}
with open(gene_lengths, 'r') as f:
    for line in f:
        gene, length = line.strip().split('\t')
        gene_lengths_dict[gene] = int(length)

outfile_matrices = {sample: f"results/08_gene_content/03_gene_counts/by_sample/{sample}_gene_counts_matrix.txt" for sample in samples}


gene_counts = pd.DataFrame(index=gene_lengths_dict.keys(), columns=samples).fillna(0)
genes = set()

# Process depth files
# for depth_file, sample in zip(depth_files, samples):
#     with open(depth_file, 'r') as file:
#         for line in file:
#             fields = line.strip().split('\t')
#             gene = fields[0]
            
#             position = int(fields[1])
#             depth = int(fields[2])

#             genes.add(gene)
#             if gene not in gene_counts[sample]:
#                 gene_counts[sample][gene] = []

#             # Adjust the depth by dividing by gene length
#             adjusted_depth = depth / gene_lengths_dict.get(gene, 1)
#             gene_counts[sample][gene].append(adjusted_depth)

# # Create the gene count matrix
# gene_count_matrix = pd.DataFrame.from_dict(gene_counts, orient='columns').fillna(0)

# # Save the gene count matrix to a file
# gene_count_matrix.to_csv(outfile, sep='\t')

# # Optional: Calculate the mean or median adjusted gene counts
# gene_counts_mean = gene_count_matrix.mean()
# gene_counts_median = gene_count_matrix.median()

# # Print the mean and median adjusted gene counts
# for sample in samples:
#     print(f"Sample: {sample}")
#     print("Mean adjusted gene counts:")
#     print(gene_counts_mean[sample])
#     print("Median adjusted gene counts:")
#     print(gene_counts_median[sample])
#     print()