#!/usr/bin/env python3
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from statistics import mean

#Usage: python3 calc_perc_id_orthologs.py Locustag_Phylotype_SDP.txt OG0000400_aln_trim.fasta

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

if (len(sys.argv) < 3):
    print('Two input-files are required for this script:')
    print('1. Tab-delim file with genome-id in tab1 and SDP-affiliation in tab 3')
    print('2. A core gene family alignment-file in fasta-format')
    print("Exiting script!")
    exit()

#Open the SDP-affiliation file and save data in dict
group_file = sys.argv[1]
with open(group_file, "r") as fh_group_file:
    genome_group = dict()
for line in fh_group_file:
    line = line.strip()
    split_line = line.split('\t')
    genome_id = "_".join(split_line[:-1])
    group = split_line[2]
    genome_group[genome_id] = group

#Open the alignment file, get ordered list of genome-ids, and save their index positions (in order to subset the alignment with genome-id downstream)
aln_in = sys.argv[2]
split_filename = aln_in.split('_')
OG_id = split_filename[0]
aln = AlignIO.read(aln_in, "fasta")
genome_ids = list()
genome_aln_index = dict()
aln_index = 0
for record in aln:
    genome_ids.append(record.id)
    genome_aln_index[record.id] = aln_index
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
    fh_outfile = open(sys.argv[3], 'a')
    line_out = [OG_id, str(0), str(0), str(0)]
    line_out_str = "\t".join(line_out)
    fh_outfile.write(line_out_str + '\n')
    fh_outfile.close()
    exit()

#Compare the genomes in each SDP to all other genomes in the alignment: calculate percentage identity for all pairwise combinations. Calculate the max, min, and mean values, print to file
fh_outfile = open(sys.argv[3], 'a')
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
fh_outfile.close()
