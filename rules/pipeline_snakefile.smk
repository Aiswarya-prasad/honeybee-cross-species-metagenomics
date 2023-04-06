#!/usr/bin/env python

"""
name: honeybee-MAGs-pipeline
description: Snakefile to be used to launch the pipeline and run qc, assembly, binning and some downstream steps
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
dependencies:
    - trim-qc
"""

import os
import itertools

configfile: "config/config.yaml"
# read config information into local variables to improve readability of rules
LOCAL_BACKUP = config["LocalBackup"]
SAMPLES_INDIA = config["SAMPLES_INDIA"]
SAMPLES_KE = config["SAMPLES_KE"]
# SAMPLES_KE = []
SAMPLES = SAMPLES_INDIA
PROJECT_IDENTIFIER = config["ProjectIdentifier"]
BACKUP_PATH = config["BackupPath"]
DBs = config["GENOME_DBs"]
ADAPTERS = config["Adapters"]
PROJECT_PATH = config["ProjectPath"]
GROUPS = ["g__Bombilactobacillus",
          "g__Lactobacillus",
          "g__Bifidobacterium",
          "g__Gilliamella",
          "g__Frischella",
          "g__Snodgrassella",
          "g__Bartonella",
          "g__Enterobacter",
          # "g__",
          "g__Pectinatus",
          "g__Apibacter",
          "g__Dysgonomonas",
          "g__Spiroplasma",
          # "g__Zymobacter",
          "g__Entomomonas",
          "g__Saezia",
          "g__Parolsenella",
          "g__WRHT01",
          "g__Commensalibacter",
          "g__Apilactobacillus",
          "g__Bombella"]
if LOCAL_BACKUP:
    localrules: backup

#  some helper functions used across rules

def get_g_dict_for_defined_groups(path):
    """
    Returns dictionary of genomes and groups with each value being a list of
    genomes corresponding to a given group
    used for the phylogenies because trees are only made for requested groups
    """
    g_list_dict = {}
    for group in GROUPS:
        g_list_dict[group] = []
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
            # only include groups of interest!
            if group not in GROUPS:
                continue
            g_list_dict[group].append(genome)
    return(g_list_dict)


def get_g_dict_for_groups_from_data(path):
    """
    Returns dictionary of genomes and groups with each value being a list of
    genomes corresponding to a given group
    used to get all MAGs from checkpoint output for all downstream steps
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
            if group in g_list_dict.keys():
                g_list_dict[group].append(genome)
            else:
                g_list_dict[group] = [genome]
    return(g_list_dict)

def get_rep_genomes_dict(path):
    """
    Returns dictionary of genomes and groups with each value being the a dict
    of magOTUs and corresponding representative genome from checkpoint output
    group1:
        cluster1: rep_genome
        cluster2: rep_genome
        cluster3: rep_genome
    group2:
        cluster4: rep_genome
    remember:
        import itertools
        people = {1: {'name': 'John', 'age': '27', 'sex': 'Male'},
         2: {'name': 'Marie', 'age': '22', 'sex': 'Female'}}
        list(itertools.chain.from_iterable(list(y.values()) for y in people.values()))
    """
    path = "06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv"
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
            rep_status = int(line.split("\t")[19])
            if rep_status != 1:
                continue
            if group == "g__":
                group = "g__"+cluster
            if group not in g_list_dict.keys():
                g_list_dict[group] = {cluster: genome}
            else:
                g_list_dict[group].update({cluster: genome})
    return(g_list_dict)

def convertToMb(string):
    """
    This function can convert text in the form
    xxG to mb
    If it does not end with G, it returns the string
    It does not handle other cases of invalid input
    """
    if string.endswith("G"):
        number = int(string.split("G")[0])
        return(number*1000)
    else:
        return(string)

def convertToSec(string):
    """
    This function can convert text in the form
    D-hh:mm:ss to seconds
    D - # days
    hh # hours
    mm # mins
    ss # secs
    """
    days = string.split("-")[0]
    hrs = string.split("-")[1].split(":")[0]
    min = string.split("-")[1].split(":")[1]
    sec = string.split("-")[1].split(":")[2]
    total = int(sec)
    total = total + 60*int(min)
    total = total + 60*60*int(hrs)
    total = total + 24*60*60*int(days)
    return(total)

def num_genomes_in_group(group, path):
    """
    This function finds the number of genomes in group
    so make_tree is only run for those with at least 3
    """
    if group in GROUPS:
        return len(get_g_dict_for_defined_groups(path)[group])
    else:
        return len(get_g_dict_for_groups_from_data(path)[group])
# path = "06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv"
# group_provided = "g__Gilliamella"
def num_magOTUs_in_group(group_provided, path):
    """
    Returns number of magOTUs within a given group
    used to get only groups for which running the "SDP_validation"
    part of the pipeline makes sense
    """
    magOTU_group_dict = {}
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            cluster = line.split("\t")[11]
            group = line.split("\t")[18]
            if group in magOTU_group_dict.keys():
                if cluster not in magOTU_group_dict[group]:
                    magOTU_group_dict[group].append(cluster)
            else:
                magOTU_group_dict[group] = [cluster]
    if group_provided in magOTU_group_dict.keys():
        return(len(magOTU_group_dict[group_provided]))
    else:
        return(0)

def get_MAGs_list(path):
    """
    read list of MAGs from checkpoint output
    (06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv) and return the list of Mags
    of a given sample either as a list per sample, as a complete list including
    all samples or a dictonary
    """
    g_list = []
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
            # sample = "_".join(genome.split("MAG_")[1].split("_")[:-1])
            g_list.append(genome)
        return(g_list)

# based on targets (?), include rules

include: "rules/time-qc.smk"
include: "rules/motus-profiling.smk"
include: "rules/assemble-qc.smk"

# to do next, 
# binning
# drep and gtdb ..
# rename MAGs (with completeness and genus in the name?) and checkpoint table
# annotation
# phylogenies
# core coverage
# sdp validation
# popcogent