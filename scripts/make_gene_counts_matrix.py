#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import numpy as np
import collections
from Bio import SeqIO
from statistics import geometric_mean

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
os.makedirs(f"{outdir}/by_sample/", exist_ok=True)

print(f"read {len(gene_lengths_dict.keys())} genes from {gene_catalog}")

gene_counts = pd.DataFrame(index=gene_lengths_dict.keys(), columns=samples)
gene_counts = pd.DataFrame(index=, columns=samples)

stats_df = pd.DataFrame(index=samples, columns=["total_genes", "filter_coverage_breadth_50", "filter_mean_mapQ_30", "filter_coverage_breadth_70", "final_mapped_mapQ_30_breadth_50"])

for coverage_file, sample in zip(coverage_files, samples):
    sample_index = samples.index(sample)
    print(f"Processing {sample} ({sample_index+1}/{len(samples)})")
    
    # Process the coverage file and extract statistics
    coverage_data = pd.read_csv(coverage_file, sep="\t")
    
    # Compute required statistics
    total_genes = len(coverage_data)
    filter_coverage_breadth_70 = len(coverage_data[coverage_data["covbases"] > 70])
    filter_coverage_breadth_50 = len(coverage_data[coverage_data["covbases"] > 50])
    filter_mean_mapQ_30 = len(coverage_data[coverage_data["meanmapq"] >= 30])
    final_mapped = len(coverage_data[(coverage_data["covbases"] > 50) & (coverage_data["meanmapq"] >= 30)])
    
    # Store statistics in the DataFrame
    stats_df.loc[sample, "total_genes"] = total_genes
    stats_df.loc[sample, "filter_coverage_breadth_70"] = filter_coverage_breadth_70
    stats_df.loc[sample, "filter_coverage_breadth_50"] = filter_coverage_breadth_50
    stats_df.loc[sample, "filter_mean_mapQ_30"] = filter_mean_mapQ_30
    stats_df.loc[sample, "final_mapped_mapQ_30_breadth_50"] = final_mapped
    
    # make a matrix of gene counts considering only genes that passed the filtering criteria
    coverage_data_filt = coverage_data[(coverage_data["covbases"] > 50) & (coverage_data["meanmapq"] >= 30)]
    size_factor = np.median(coverage_data_filt["meandepth"] / np.exp(np.mean(np.log(coverage_data_filt["meandepth"]))))
    
    
    # gene_counts_data = coverage_data.set_index("#rname")["meandepth"]
    gene_counts_data = coverage_data.set_index("#rname")["meandepth"] * size_factor
    # Update gene counts DataFrame with the mapped counts
    gene_counts[sample] = gene_counts_data.reindex(gene_counts.index)
    # write gene_counts[sample] to outfile_matrices[sample]
    gene_counts[sample].to_csv(outfile_matrices[sample], sep="\t", header=False)
# replace nan with 0 in gene_counts
gene_counts.to_csv(outfile, sep="\t", header=True)

stats_df.to_csv(filtering_info, sep="\t", header=True)

metadata="results/09_MAGs_collection/All_mags_sub/All_mags_sub_metadata_summary.tsv"

outfile_kos = "results/08_gene_content/03_gene_counts/ko_family_genes_matrix.txt"

# collect dram annotation deatails for genes in the gene catalog
annotations_dict = {name: "NA" for name in gene_counts.index}
dram_annotations = glob.glob("results/08_gene_content/02_DRAM_annotations/*/annotations.tsv")
for annotation_file in dram_annotations:
    with open(annotation_file, "r") as fh:
        header = fh.readline()
        kegg_ind = header.split("\t").index("kegg_id")
        for line in fh:
            print(line)
            name = "".join(line.split("\t")[0].split("_orfs_")[1])
            if name in annotations_dict.keys():
                kegg_id_text = line.split("\t")[kegg_ind]
                if kegg_id_text:
                    annotations_dict[name] = kegg_id_text
# correct
# make a matrix by summing up meandepth of filtered genes that belong to the same KO family
gene_counts_ko = pd.DataFrame(index=list(set(annotations_dict.values())), columns=gene_counts.columns)
gene_counts_ko = gene_counts_ko.fillna(0)
for gene in gene_counts.index:
    ko = annotations_dict[gene]
    if ko != "NA":
        gene_counts_ko.loc[ko] += gene_counts.loc[gene]
# print gene_counts_ko not nan
gene_counts_ko = gene_counts_ko[gene_counts_ko.notna().any(axis=1)]

gene_counts_ko.to_csv(outfile_kos, sep="\t", header=True)

outfile_GHs = "results/08_gene_content/03_gene_counts/GH_family_genes_matrix.txt"

annotations_dict_gh = {name: "NA" for name in gene_counts.index}
dram_annotations = glob.glob("results/08_gene_content/02_DRAM_annotations/*/annotations.tsv")
for annotation_file in dram_annotations:
    with open(annotation_file, "r") as fh:
        header = fh.readline()
        cazy_ind = header.split("\t").index("cazy_hits")
        for line in fh:
            name = "".join(line.split("\t")[0].split("_orfs_")[1])
            if name in annotations_dict_gh.keys():
                cazy_id_text = line.split("\t")[cazy_ind]
                if cazy_id_text:
                    annotations_dict_gh[name] = cazy_id_text
# annotations_dict_gh which is not na
# extract GH families in the format GHxx from annotations_dict_gh
def get_gh_family(x):
    if "[GH" in x:
        num = x.split("[GH")[1].split("]")[0]
        return f"GH{num}"
    else:
        return x
# make clean annotation gh with get_gh_family
annotations_dict_gh_clean = {key: get_gh_family(x) if "[GH" in x else "NA" for key, x in annotations_dict_gh.items()}

gene_counts_GH = pd.DataFrame(index=list(set(annotations_dict_gh_clean.values())), columns=gene_counts.columns)
gene_counts_GH = gene_counts_GH.fillna(0)
for gene in gene_counts.index:
    GH = annotations_dict_gh_clean[gene]
    if GH != "NA":
        gene_counts_GH.loc[GH] += gene_counts.loc[gene]
# print gene_counts_GH not nan
gene_counts_GH = gene_counts_GH[gene_counts_GH.notna().any(axis=1)]
# print gene_counts_GH not 0
gene_counts_GH = gene_counts_GH[(gene_counts_GH.T != 0).any()]
gene_counts_GH.to_csv(outfile_GHs, sep="\t", header=True)
