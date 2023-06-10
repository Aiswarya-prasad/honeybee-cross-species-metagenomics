#!/usr/bin/env python3

import os
from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse

# parser = argparse.ArgumentParser()
# requiredNamed = parser.add_argument_group('required arguments')
# requiredNamed.add_argument('--out_',metavar="db_metafile",required=True, help="File from checkpoint output with all mag details", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
# requiredNamed.add_argument('--ortho', metavar='orthofinder_file', required=True, help="Filtered single-copy ortholog gene file (orthofinder format)", action="store")
# requiredNamed.add_argument('--beddir', metavar="bedfiles_dir", required=True, help="Directory containing bed-files for genomes in db", action="store")
# requiredNamed.add_argument('--bamfile', metavar="bamfile_notlistanymore", help="List of bam-files, one line per file", action="store")
# requiredNamed.add_argument('--group', metavar="group", required=True, help="Name of sample being evaluated", action="store")
# requiredNamed.add_argument('--sample', metavar="sample", required=True, help="Name of genus / group being evaluated", action="store")
# requiredNamed.add_argument('--outfile', metavar="out_file", required=True, help="Output file containing coverage of genefamilies in each magOTU", action="store")
# args = parser.parse_args()

# database_metafile = args.info
# ortho_file = args.ortho
# bedfiles_dir = args.beddir

samples = ["M1.1", "M1.2", "M1.3", "M1.4", "M1.5",
              "C1.1", "C1.2", "C1.3", "C1.4", "C1.5",
              "C2.1", "C2.2", "C2.3", "C2.4", "C2.5",
              "C3.1", "C3.2", "C3.3", "C3.4", "C3.5",
              "D1.1","D1.2","D1.3","D1.4","D1.5",
              "D2.1","D2.2","D2.3","D2.4","D2.5",
              "D3.1","D3.2","D3.3","D3.4","D3.5",
              "F1.1","F1.2","F1.3","F1.4","F1.5",
              "F2.1","F2.2","F2.3","F2.4","F2.5",
              "F3.1","F3.2","F3.3","F3.4","F3.5"
              ]
contig_fates_dict = {}
for sample in samples:
    contig_fates = "06_MAG_binning/contig_fates/"+sample+"_contig_fates.csv"
    with open(contig_fates, "r") as fates_fh:
        header = fates_fh.readline()
        for line in fates_fh:
            line = line.strip()
            sample = line.split(",")[0]
            contig_name = line.split(",")[1]
            bin_name = line.split(",")[6]
            if sample not in contig_fates_dict.keys():
                contig_fates_dict[sample] = {contig_name: bin_name}
            else:
                contig_fates_dict[sample][contig_name] = bin_name

MAGs_info = pd.read_csv("06_MAG_binning/MAGs_GenomeInfo_auto.tsv", sep ="\t")
All_Groups = [grp for grp in MAGs_info.Group_auto.unique() if grp == grp]
All_Groups.append("NA")

genes_in_genus_annotated = {}
genes_in_binned = {}
genes_in_unbinned = {}
genes_in_total = {}

out_num_genes = "/scratch/aprasad/211018_Medgenome_india_samples/Figures/Number_genes_binned.csv"

with open(out_num_genes, "w") as out_num_fh:
    out_num_fh.write(f"Sample, Num_genes, Num_genes_binned, Num_genes_unbinned, Perc_genes_binned")
    print(f"Sample, Num_genes, Num_genes_binned, Num_genes_unbinned, Perc_genes_binned")

    for sample in samples:
        genes_in_binned[sample] = 0
        genes_in_unbinned[sample] = 0
        genes_in_total[sample] = 0
        ORFs_ffn = "12_species_validation/metagenomic_orfs/"+sample+"_orfs.ffn"
        records = SeqIO.parse(ORFs_ffn, "fasta")
        for record in records:
            gene_id = record.id
            contig_name = "_".join(record.id.split(sample+"_")[1].split("_")[:-1])
            genes_in_total[sample] += 1
            bin_name = contig_fates_dict[sample][contig_name]
            if bin_name == "NA":
                genes_in_unbinned[sample] += 1
            else:
                genes_in_binned[sample] += 1
                if bin_name in All_Groups:
                    genes_in_genus_annotated[sample] += 1
    for sample in samples:
        out_num_fh.write(f"{sample}, {genes_in_total[sample]}, {genes_in_binned[sample]}, {genes_in_unbinned[sample]}, {genes_in_binned[sample]/genes_in_total[sample]*100}")
        print(f"{sample}, {genes_in_total[sample]}, {genes_in_binned[sample]}, {genes_in_unbinned[sample]}, {genes_in_binned[sample]/genes_in_total[sample]*100}")

# num_genes_by_group = {}
# num_orfs_in_mag = {}
# for mag in MAGs_info["ID"]:
#     num_orfs_in_mag[mag] = 0
# for sample in samples:
#     num_genes_by_group[sample] = {}
#     for group in All_Groups:
#         num_genes_by_group[sample][group] = 0
#     ORFs_ffn = "12_species_validation/metagenomic_orfs/"+sample+"_orfs.ffn"
#     for record in SeqIO.parse(ORFs_ffn, "fasta"):
#         gene_id = record.id
#         contig_name = "_".join(record.id.split(sample+"_")[1].split("_")[:-1])
#         bin_name = contig_fates_dict[sample][contig_name]
#         if bin_name == "NA":
#             pass
#         else:
#             num_orfs_in_mag[bin_name] += 1
#             Genus = MAGs_info.query("ID == @bin_name")["Group_auto"].values[0]
#             if Genus != Genus:
#                 Genus = "NA"
#             num_genes_by_group[sample][Genus] += 1

# for sample in samples:
#     print(f"{num_genes_by_group[sample]}")
