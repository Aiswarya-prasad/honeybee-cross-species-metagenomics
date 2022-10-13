#!/usr/bin/env python3
import os
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from statistics import mean
import argparse

#Usage: python3 calc_perc_id_orthologs.py Locustag_Phylotype_SDP.txt OG0000400_aln_trim.fasta

def get_genome_group_dict(path):
    """
    Returns dictionary of genomes and groups with each value being a the group
    correspoinding to the key which is the genome
    """
    g_list_dict = {}
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            cluster = line.split("\t")[11]
            group = line.split("\t")[18]
            # only include groups of interest!
            if group == "g__":
                group = "g__"+cluster
            g_list_dict[genome] = group
    return(g_list_dict)

#Function for getting the elements in list1 which are not contained in list2
def Diff(list1, list2):
    return (list(list(set(list1)-set(list2)) + list(set(list2)-set(list1))))

#Function for calculating percentage identity in alignment with two sequences (skipping positions with gaps in either seq)
def pairwise_perc_id(aln):
    nb_seqs = len(aln)
    nb_col = aln.get_alignment_length()
    nb_col_ungapped = 0
    nb_matches = 0
    i = 0
    gap_char = '-'
    selected_columns = list()
    while (i < nb_col):
        aln_col = aln[:,i]
        nb_gaps = aln_col.count(gap_char)
        if nb_gaps == 0:
            nb_col_ungapped += 1
            aln_col_str = aln_col.upper()
            if aln_col_str[0] == aln_col_str[1]:
                nb_matches += 1
        i += 1
    pairwise_perc_id = round(nb_matches/nb_col_ungapped, 3)
    return(pairwise_perc_id)

#Function for generating pairwise alignments between two groups of genomes, and calculating the percentage identity in each alignment. Will return list of perc-ids for all pairwise combinations between the two groups
def group_perc_id(in_group, out_group):
    group_perc_ids = list()
    for genome_in in in_group:
        for genome_out in out_group:
            aln_seqrecords = list()
            aln_seqrecords.append(aln[genome_aln_index[genome_in]])
            aln_seqrecords.append(aln[genome_aln_index[genome_out]])
            pairwise_aln = MultipleSeqAlignment(aln_seqrecords)
            perc_id = pairwise_perc_id(pairwise_aln)
            group_perc_ids.append(perc_id)
    return(group_perc_ids)

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('--meta', action='store', help='Tab-delim file with genome-id in tab1 and cluster affiliation - parsed assuming that it is the checkpoint file. Change script if it is not.')
requiredNamed.add_argument('--trim_file', action='store', help='Core gene family alignment-file in fasta-format trimmed - description to be updated')
requiredNamed.add_argument('--outfile', action='store', help='Output file path - description to be updated')
args = parser.parse_args()

genome_db_meta = args.meta
aln_in = args.trim_file
outfile = args.outfile
#Open the SDP-affiliation file and save data in dict
genome_group = get_genome_group_dict(genome_db_meta)
#Open the alignment file, get ordered list of genome-ids, and save their index positions (in order to subset the alignment with genome-id downstream)
split_filename = aln_in.split('_')
OG_id = split_filename[0]
aln = AlignIO.read(aln_in, "fasta")
genome_ids = list()
genome_aln_index = dict()
aln_index = 0
for record in aln:
    recored_unparsed = record.id
    record_parsed = "_".join(recored_unparsed.split(" ")[0].split("_")[:-1])
    genome_ids.append(record_parsed)
    genome_aln_index[record_parsed] = aln_index
    aln_index += 1

#Identify the SDP-groups contained within the alignment, save in dict
group_genomes = dict()
for genome_id in genome_ids:
    group = genome_group[genome_id]
    if group not in group_genomes:
        group_genomes[group] = list()
    group_genomes[group].append(genome_id)

#Check the number of SDPs contained within the alignment. If more than one, continue by calculating alignment percentage identity stats across SDPs. If only one SDP, exit script.
nb_groups = len(group_genomes)
if nb_groups == 1:
    with open(outfile, "a") as fh_outfile:
        line_out = [OG_id, str(0), str(0), str(0)]
        line_out_str = "\t".join(line_out)
        fh_outfile.write(line_out_str + '\n')
        exit()

#Compare the genomes in each SDP to all other genomes in the alignment: calculate percentage identity for all pairwise combinations. Calculate the max, min, and mean values, print to file
with open(outfile, "a") as fh_outfile:
    all_perc_ids = list()
    for group in group_genomes.keys():
        in_genomes = group_genomes[group]
        out_genomes = Diff(genome_ids, in_genomes)
        group_perc_ids = group_perc_id(in_genomes, out_genomes)
        all_perc_ids.extend(group_perc_ids)
    mean_perc_id = round(mean(all_perc_ids),3)
    min_perc_id = min(all_perc_ids)
    max_perc_id = max(all_perc_ids)
    line_out = [OG_id, str(mean_perc_id), str(min_perc_id), str(max_perc_id)]
    line_out_str = "\t".join(line_out)
    fh_outfile.write(line_out_str + '\n')
