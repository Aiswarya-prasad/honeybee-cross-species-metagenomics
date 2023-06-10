#!/usr/bin/env python3

import os
import sys
import re
import argparse
import pickle
import glob
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

print(f"starting to process {cluster_file}")

# cluster_file="results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta.clstr"
# cluster_out="results/08_gene_content/00_cdhit_clustering/cluster_info.tsv"
                
pickle_file = cluster_file + ".pickle"

if os.path.exists(pickle_file):
    print(f"loading pickle file from {pickle_file}")
    with open(pickle_file, "rb") as pkl_fh:
        cluster_groups = pickle.load(pkl_fh)
else:
    with open(cluster_file, "r") as cfh:
            cluster_groups = [
                list(group)
                for key, group in groupby(cfh, lambda line: line.startswith(">Cluster"))
                if not key
            ]
    with open(pickle_file, "wb") as pkl_fh:
        pickle.dump(cluster_groups, pkl_fh)

print("cluster information loaded")

progress_marker = cluster_out + ".progress"

if os.path.exists(progress_marker) and os.path.exists(cluster_out):
    with open(progress_marker, "r") as pm_fh:
        already_processed = int(pm_fh.read())
else:
    already_processed = 0

def save_progress(progress_marker, i):
    with open(progress_marker, "w+") as pm_fh:
        pm_fh.write(str(i))

total_clusters = len(cluster_groups)
def get_gh_family(x):
    if "[GH" in x:
        num = x.split("[GH")[1].split("]")[0]
        return f"GH{num}"
    else:
        return x

print("reading annotations gh")
pickle_file_gh = os.path.join(os.path.dirname(cluster_out), "annotations_gh.pickle")
if os.path.exists(pickle_file_gh):
    print(f"loading pickle file from {pickle_file_gh}")
    with open(pickle_file_gh, "rb") as pkl_fh:
        annotations_dict_gh = pickle.load(pkl_fh)
else:
    annotations_dict_gh = {}
    dram_annotations = glob.glob("results/08_gene_content/02_DRAM_annotations/*/annotations.tsv")
    for annotation_file in dram_annotations:
        with open(annotation_file, "r") as fh:
            header = fh.readline()
            cazy_ind = header.split("\t").index("cazy_hits")
            for line in fh:
                name = "".join(line.split("\t")[0].split("_orfs_")[1])
                cazy_id_text = line.split("\t")[cazy_ind]
                if cazy_id_text and "[GH" in cazy_id_text:
                    annotations_dict_gh[name] = get_gh_family(cazy_id_text)
                else:
                    annotations_dict_gh[name] = "NA"
    with open(pickle_file_gh, "wb") as pkl_fh:
        pickle.dump(annotations_dict_gh, pkl_fh)

print("reading annotations ko")
pickle_file_ko = os.path.join(os.path.dirname(cluster_out), "annotations_ko.pickle")
if os.path.exists(pickle_file_ko):
    print(f"loading pickle file from {pickle_file_ko}")
    with open(pickle_file_ko, "rb") as pkl_fh:
        annotations_dict_ko = pickle.load(pkl_fh)
else:
    annotations_dict_ko = {}
    dram_annotations = glob.glob("results/08_gene_content/02_DRAM_annotations/*/annotations.tsv")
    for annotation_file in dram_annotations:
        with open(annotation_file, "r") as fh:
            header = fh.readline()
            kegg_ind = header.split("\t").index("kegg_id")
            for line in fh:
                name = "".join(line.split("\t")[0].split("_orfs_")[1])
                kegg_id_text = line.split("\t")[kegg_ind]
                if kegg_id_text:
                    annotations_dict_ko[name] = kegg_id_text
                else:
                    annotations_dict_ko[name] = "NA"
    with open(pickle_file_ko, "wb") as pkl_fh:
        pickle.dump(annotations_dict_ko, pkl_fh)

print("Starting to write clusters info to file")

with open(cluster_out, 'a+') as out_fh:
    if already_processed > total_clusters:
        sys.exit("ERROR: already processed more clusters than total clusters")
    if already_processed == 0:
        header_line = "\t".join(['cluster', 'num_genes', 'num_hosts', 'hosts', 'GH_list', 'KO_list', 'M', 'C', 'D', 'F', 'A', 'mean_length', 'standard_deviation', 'min_length', 'max_length', 'length_representative', 'representative_gene', 'representative_KO', 'representative_GH', 'num_GHs', 'num_KOs'])
        out_fh.write(f"{header_line}\n")
    if already_processed == total_clusters:
        sys.exit("DONE: already processed all clusters!")
    else:
        print(f"{already_processed} clusters already processed. Will begin/continue from cluster {already_processed + 1}")
    for i, genes in enumerate(cluster_groups, start=1):
        if i <= already_processed:
            continue
        print(f"Processing cluster {i} of {total_clusters}", end="\r")
        cluster_info = {}
        length_dict = {key.split("...")[0].split(">")[1]: int(value.split('nt')[0]) for key, value in zip([g.split()[2] for g in genes], [g.split()[1] for g in genes]) }
        ids = length_dict.keys()
        rep_id = [g.split()[2].split("...")[0].split(">")[1] for g in genes if g.strip().endswith("*")].pop()
        cluster_info['cluster'] = i
        cluster_info['num_genes'] = len(ids)
        hosts = [id.split("_NODE")[0][0] for id in ids]
        cluster_info['num_hosts'] = len(set(hosts))
        cluster_info['hosts'] = ",".join(set(hosts))
        ghs = [annotations_dict_gh[id] for id in ids]
        kos = [annotations_dict_ko[id] for id in ids]
        cluster_info['GHs'] = ",".join(set(ghs))
        cluster_info['KOs'] = ",".join(set(kos))
        cluster_info['M'] = 1 if hosts.count('M') > 0 else 0
        cluster_info['C'] = 1 if hosts.count('C') > 0 else 0
        cluster_info['D'] = 1 if hosts.count('D') > 0 else 0
        cluster_info['F'] = 1 if hosts.count('F') > 0 else 0
        cluster_info['A'] = 1 if hosts.count('A') > 0 else 0
        cluster_info['mean_length'] = np.mean(list(length_dict.values()))
        cluster_info['standard_deviation'] = np.std(list(length_dict.values()))
        cluster_info['min_length'] = np.min(list(length_dict.values()))
        cluster_info['max_length'] = np.max(list(length_dict.values()))
        cluster_info['length_representative'] = length_dict[rep_id]
        cluster_info['representative_gene'] = rep_id
        cluster_info['representative_KO'] = annotations_dict_ko[rep_id]
        cluster_info['representative_GH'] = annotations_dict_gh[rep_id]
        cluster_info['num_GHs'] = len(set(ghs) - set(["NA"]))
        cluster_info['num_KOs'] = len(set(kos)  - set(["NA"]))
        outstring = "\t".join([str(x) for x in cluster_info.values()])
        # write to outfile by line
        out_fh.write(f"{outstring}\n")
        if i % 10000 == 0:
            save_progress(progress_marker, i)
    save_progress(progress_marker, i)