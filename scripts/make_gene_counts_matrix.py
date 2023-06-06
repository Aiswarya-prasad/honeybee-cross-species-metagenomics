#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import numpy as np
import collections
from Bio import SeqIO

"""
collect all the gene counts from the coveragefiles in the directory provided
and create a matrix of gene counts mxn where n is the number of genes and m is the number of samples
each row is a gene and each column is a sample
"""


parser = argparse.ArgumentParser(description="create a matrix of gene counts from coveragefiles")
parser.add_argument("-d", "--coveragefiles_dir", help="directory containing coveragefiles", required=True)
parser.add_argument("-s", "--suffix", help="pattern appended after sample name to coverage file", required=True)
parser.add_argument("-l", "--samples_list", nargs = "+", help="optional argument to specify samples", required=False)
parser.add_argument("-g", "--gene_lengths", help="file to write gene lengths to", required=True)
parser.add_argument("-f", "--filtering_info", help="file to write filtering results to", required=True)
parser.add_argument("-c", "--gene_catalog", help="gene catalog fasta that all samples were mapped to", required=True)
parser.add_argument("-o", "--outfile", help="output file name", required=True)
args = parser.parse_args()

def make_gene_lengths_file(catalog, lengths_file):
    with open(catalog, "r") as catalog_file:
        with open(lengths_file, "w+") as lengths_fh:
            for record in SeqIO.parse(catalog_file, "fasta"):
                lengths_fh.write(f"{record.id}\t{len(record.seq)}\n")
    
coverages_dir = args.coveragefiles_dir
suffix = args.suffix
if args.samples:
    samples = args.samples
    coverage_files = [coverages_dir + sample + suffix for sample in samples]
else:
    coverage_files = glob.glob(coverages_dir + "*" + suffix)
    samples = [os.path.basename(coverage_file).split(suffix)[0] for coverage_file in coverage_files]
outfile = args.outfile
outdir = os.path.dirname(outfile)
gene_catalog = args.gene_catalog
gene_lengths = args.gene_lengths
filtering_info = args.filtering_info


coverages_dir = "results/08_gene_content/01_profiling/"
suffix = "_mapped.coverage"
coverage_files = glob.glob(coverages_dir + "*" + suffix)
samples = [os.path.basename(coverage_file).split(suffix)[0] for coverage_file in coverage_files]
outfile = "results/08_gene_content/03_gene_counts/all_genes_matrix.txt"
outdir = os.path.dirname(outfile)
gene_catalog = "results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta"
gene_lengths = "results/08_gene_content/03_gene_counts/gene_lengths.txt"
filtering_info = "results/08_gene_content/03_gene_counts/gene_filtering_info.txt"
os.makedirs(os.path.dirname(gene_lengths), exist_ok=True)

# make_gene_lengths_file(gene_catalog, gene_lengths)

gene_lengths_dict = {}
with open(gene_lengths, 'r') as f:
    for line in f:
        gene, length = line.strip().split('\t')
        gene_lengths_dict[gene] = int(length)

outfile_matrices = {sample: f"{outdir}/by_sample/{sample}_gene_counts_matrix.txt" for sample in samples}

print(f"read {len(gene_lengths_dict.keys())} genes from {gene_catalog}")
gene_counts = pd.DataFrame(index=gene_lengths_dict.keys(), columns=samples)

stats_df = pd.DataFrame(index=samples, columns=["total_genes", "filter_coverage_breadth_50", "filter_mean_mapQ_30", "final_mapped"])

for coverage_file, sample in zip(coverage_files, samples):
    sample_index = samples.index(sample)
    print(f"Processing {sample} ({sample_index+1}/{len(samples)})")
    
    # Process the coverage file and extract statistics
    coverage_data = pd.read_csv(coverage_file, sep="\t", header=None, names=["rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"])
    
    # Compute required statistics
    total_genes = len(gene_lengths_dict)
    filter_coverage_breadth_50 = (coverage_data["covbases"] > 0).mean() * 100
    filter_mean_mapQ_30 = (coverage_data["meanmapq"] >= 30).mean() * 100
    final_mapped = ((coverage_data["covbases"] > 0) & (coverage_data["meanmapq"] >= 30)).sum()
    
    # Store statistics in the DataFrame
    stats_df.loc[sample] = [total_genes, filter_coverage_breadth_50, filter_mean_mapQ_30, final_mapped]

# Print the statistics DataFrame
print(stats_df)