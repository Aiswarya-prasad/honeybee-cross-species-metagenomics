#!/usr/bin/env python3

import os
import sys
import re
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from itertools import groupby


"""
python3 scripts/parse_cluster_file.py --cluster_file {input.cdhit_clustering} --cluster_out {output.cdhit_clusters}
"""

parser = argparse.ArgumentParser(description="parse cdhit cluster file")
parser.add_argument("-c", "--cluster_file", help="cdhit cluster file", required=True)
parser.add_argument("-o", "--cluster_out", help="output file name", required=True)
args = parser.parse_args()

cluster_file = args.cluster_file
cluster_out = args.cluster_out


cluster_file = "results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta.clstr"
cluster_out = "results/08_gene_content/00_cdhit_clustering/cluster_host_affiliations.tsv"
                
# create empty dataframe to store number of genes and number of hosts and a column
# indicating if present in given host (m, c, d, f, a) for each cluster

with open(cluster_file, "r") as cfh:
        cluster_groups = [
            list(group)
            for key, group in groupby(cfh, lambda line: line.startswith(">Cluster"))
            if not key
        ]

clusters_df = pd.DataFrame(columns=['cluster', 'num_genes', 'num_hosts', 'M', 'C', 'D', 'F', 'A', 'mean_length', 'standard_deviation', 'min_length', 'max_length', 'length_representative', 'representative_gene'])
total_clusters = len(cluster_groups)
for i, genes in enumerate(cluster_groups, start=1):
    print(f"Processing cluster {i} of {total_clusters}", end="\r")
    length_dict = {key.split("...")[0].split(">")[1]: int(value.split('nt')[0]) for key, value in zip([g.split()[2] for g in genes], [g.split()[1] for g in genes]) }
    ids = length_dict.keys()
    rep_id = [g.split()[2].split("...")[0].split(">")[1] for g in genes if g.strip().endswith("*")].pop()
    clusters_df.loc[i, 'cluster'] = i
    clusters_df.loc[i, 'num_genes'] = len(ids)
    hosts = [id.split("_NODE")[0][0] for id in ids]
    clusters_df.loc[i, 'num_hosts'] = len(set(hosts))
    clusters_df.loc[i, 'M'] = 1 if hosts.count('M') > 0 else 0
    clusters_df.loc[i, 'C'] = 1 if hosts.count('C') > 0 else 0
    clusters_df.loc[i, 'D'] = 1 if hosts.count('D') > 0 else 0
    clusters_df.loc[i, 'F'] = 1 if hosts.count('F') > 0 else 0
    clusters_df.loc[i, 'A'] = 1 if hosts.count('A') > 0 else 0
    clusters_df.loc[i, 'mean_length'] = np.mean(list(length_dict.values()))
    clusters_df.loc[i, 'standard_deviation'] = np.std(list(length_dict.values()))
    clusters_df.loc[i, 'min_length'] = np.min(list(length_dict.values()))
    clusters_df.loc[i, 'max_length'] = np.max(list(length_dict.values()))
    clusters_df.loc[i, 'length_representative'] = length_dict[rep_id]
    clusters_df.loc[i, 'representative_gene'] = rep_id
    
clusters_df.to_csv(cluster_out, sep="\t", index=False)