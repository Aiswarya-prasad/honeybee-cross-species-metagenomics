#!/usr/bin/env python3

"""
Edited version of GSL's SAGE_FetchAndConcat2.py to work with snakemake
# parse table and extract genes ID per genomes
# dict to store genomes names as key and  dict key= all values sequences as value
#%%
"""

# modules

import os
from os import path
import sys
import shutil
import numpy as np
import pandas as pd
import argparse

#%%
#functions

# fasta = OG
# Path = PathAlignedOG
# pipeNames

def fasta_to_dict(fasta, Path, pipeNames) :
    dict_OG = {}
    seq = ''
    pathToFaa = Path + '/' + fasta
    with open(pathToFaa) as f:
        OG_ID = LongName(f.readline(), pipeNames)
        for line in f :
            if line.startswith(">") :
                print(line)
                dict_OG[OG_ID] = list(seq)
                OG_ID = LongName(line, pipeNames)
                seq = ''
            else :
                seq += line.rstrip()
            print(OG_ID)
            dict_OG[OG_ID] = list(seq)
    return(dict_OG)


def LongName(string, pipeNames='False') :
    """
    LongName function needs to be edited according to the format of the annotation output
    Currently it accounts for '_' in sample names and assumes that the words in the annotation
    output header are separated by ' '
    eg. MAG_D3.2_9_00823 Fructose import permease protein FruG_3
    and returns MAG_D3.2_9
    """
    if len(string.rstrip().split('_'))>2 and pipeNames==False :
        ID = '_'.join(string.rstrip().split(' ')[0].split('_')[0:-1])
        return(ID)
    elif pipeNames=='True' :
        ID='_'.join(string.rstrip().split('|')[0].split('_')[0:-1])
        return(ID)
    else :
        ID = string.rstrip().split('_')[0]
        return(ID)

def insert_newlines(string, every=64) :
    return ('\n'.join(string[i:i + every] for i in range(0, len(string), every)))

def get_g_dict_for_groups(path):
    """
    Returns dictionary of genomes and groups with each value being a list of
    genomes corresponding to a given group
    """
    g_list_dict = {}
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
            group = line.split("\t")[18]
            if group not in g_list_dict.keys():
                g_list_dict[group] = []
            g_list_dict[group].append(genome)
    return(g_list_dict)

parser = argparse.ArgumentParser()
parser.add_argument('--aligned_dir', action='store', help='name of directory where input alignments are')
parser.add_argument('--pruned_dir', action='store', help='name of directory to write outputs to')
parser.add_argument('--pipe_names', action='store', help='True if names have pipes')
parser.add_argument('--genomes_list', action='store', help='True if names have pipes')
parser.add_argument('--group', action='store', help='Group wildcard passed from snakemake')
args = parser.parse_args()

PathAlignedOG=args.aligned_dir
OutDir=args.pruned_dir
WriteAll=True
pipeNames=args.pipe_names
if pipeNames == 'False':
    pipeNames = False
genomes_file=args.genomes_list
group=args.group

# PathAlignedOG=snakemake.input.aligned_dir
# OutDir=snakemake.output.pruned_dir
# WriteAll=True
# pipeNames=snakemake.params.pipe_names
# if pipeNames == 'False':
#     pipeNames = False
# genomes_file=snakemake.input.genomes_list

Genomes = get_g_dict_for_groups(genomes_file)[group]
Genomes = [">"+genome for genome in Genomes]
print(f"{Genomes}")

if not os.path.exists(OutDir):
    os.makedirs(OutDir)

OGList = [item for item in os.listdir(PathAlignedOG) if not item.startswith(".") and item.endswith("_aligned.faa")]
AlignedFiles=[OG for OG in OGList if os.path.isfile(PathAlignedOG+'/'+OG)]
MainDictionary={}
for OG in AlignedFiles :
    OGDictionary = fasta_to_dict(OG, PathAlignedOG, pipeNames)
    MainDictionary[OG]=OGDictionary

# Genomes = ['A-1-12', 'A12', 'A2', 'Aw-20', 'E1', 'ESL0196', 'ESL0304', 'ESL0323', 'ESL0324', 'MS1-3', 'N-S3', 'N-S4', 'N-W7', 'N9', 'PEB0171', 'wkB2', 'wkB339', 'wkB237', 'wkB273', 'wkB298', 'ESL0253', 'HK3', 'Nev4-2', 'Occ4-2', 'wkB12', 'DSM_17694', 'DSM_2580', 'MWU2920', 'ATCC_15551', 'HLGZ1', 'DSM_17474', 'NCTC_13336', 'ATCC23332', 'DSM_16848', 'NCTC10283', 'MAG_C1.1_17', 'MAG_C1.1_2', 'MAG_C1.2_15', 'MAG_C1.3_13', 'MAG_C1.3_2', 'MAG_C1.4_17', 'MAG_C1.5_6', 'MAG_C1.5_9', 'MAG_C2.1_5', 'MAG_C2.1_7', 'MAG_C2.2_12', 'MAG_C2.3_19', 'MAG_C2.3_5', 'MAG_C2.4_11', 'MAG_C2.4_13', 'MAG_C2.5_18', 'MAG_C3.1_1', 'MAG_C3.1_3', 'MAG_C3.2_13', 'MAG_C3.2_9', 'MAG_C3.3_14', 'MAG_C3.3_6', 'MAG_C3.4_21', 'MAG_C3.5_2', 'MAG_D1.4_34', 'MAG_D2.3_21', 'MAG_M1.1_3', 'MAG_M1.1_7', 'MAG_M1.2_7', 'MAG_M1.2_9', 'MAG_M1.3_12', 'MAG_M1.3_3', 'MAG_M1.4_18', 'MAG_M1.5_11', 'MAG_M1.5_6']
# OG = "OG0000997_aligned.faa"
# PathAlignedOG = "."
# pipeNames = False
# alnOG = OG

#%%
# 1 : MainDictionary key(OG number) : value(OG_ID) (dictionary)
# 2 : OGDictionary = key(genome) : value(sequence)
allPrunedSequences = {}
FIRST='True'
gap = {'-'}
# ProteinTable = pd.DataFrame(AlignedFiles)
# ProteinTable.insert(0, 'Genome',AlignedFiles , True)
# i=0

for alnOG in sorted(MainDictionary.keys()) :
    # i+=1
    print(f"###..........Processing {alnOG}..........###")
    seq_table = pd.DataFrame.from_dict(MainDictionary[alnOG], orient='columns')
    gap_count = seq_table.isin(gap).sum(1)/seq_table.shape[1]
    """
    # to allow MAGs to be kept as some may have too many gaps... So no pruning
    prunedSeq_table = seq_table[gap_count<0.5]
    pruned_seq = prunedSeq_table.to_dict('list')
    """
    pruned_seq = seq_table.to_dict('list')
    pruned_seq = {key : ''.join(pruned_seq[key]) for key in pruned_seq.keys()}
    """
    check if all genomes are present in this OG
    if all genomes are not present, add as many gaps as the alignment length
    """
    if len(Genomes) > len(pruned_seq.keys()):
        print(f"{alnOG} only has {len(pruned_seq.keys())} out of {len(Genomes)} genomes")
        missing_genome_names = [item.split(">")[1] for item in Genomes if item not in pruned_seq.keys()]
        print(f"missing genomes are {missing_genome_names} adding columns for these")
        missing_genomes = [item for item in Genomes if item not in pruned_seq.keys()]
        alignment_length = [len(MainDictionary[alnOG][item]) for item in MainDictionary[alnOG].keys()][0]
        for genome in missing_genomes:
            pruned_seq[genome] = list(gap)[0]*alignment_length
    else:
        print(f"{alnOG} is in all the {len(Genomes)} genomes")
    if WriteAll==True:
        pruned_aln = OutDir + '/' + alnOG.split(".")[0] + '_pruned.fasta'
        with open(pruned_aln, 'w') as pruned_alignment :
            for key, value in pruned_seq.items() :
                print(key + '\n')
                print(insert_newlines(value) + '\n')
                pruned_alignment.write(key + '\n')
                pruned_alignment.write(insert_newlines(value) + '\n')


    for key, value in pruned_seq.items() :
        print(key)
        if FIRST=='True' :
            allPrunedSequences[key] = (insert_newlines(value) + '\n')

        else :
            print('else')
            allPrunedSequences[key] += insert_newlines(value) + '\n'
    FIRST='False'
#%%
full_aln =  OutDir + '/' + 'CoreGeneAlignment.fasta'
with open(full_aln, 'w') as pruned_alignment :
    for key, value in allPrunedSequences.items() :
        GenomeName = key
        pruned_alignment.write(GenomeName + '\n')
        pruned_alignment.write(value)
   # list_allseq.append(list_seq)
