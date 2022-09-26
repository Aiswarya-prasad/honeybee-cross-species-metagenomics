#!/usr/bin/env python3
import sys
import os

# adapted from Kirsten's get_singlecp_orthologs.py

def is_scp_core(gene_list, genomes, MAGs):
    is_core = 0
    is_single = True
    num_genomes = len(genomes)
    num_MAGs = len(MAGs)
    match_count_genomes = dict()
    match_count_MAGs = dict()
    genomes_seen = set()
    for gene in gene_list:
        split_gene = gene.split('_')
        genome_id = "_".join(split_gene[:-1])
        if genome_id in genomes_seen:
            is_single = False
        else:
            genomes_seen.add(genome_id)
        if genome_id not in MAGs:
            match_count_genomes[genome_id] = match_count_genomes.get(genome_id,0) + 1
        else:
            match_count_MAGs[genome_id] = match_count_MAGs.get(genome_id,0) + 1
    nb_matches_MAGs = len(match_count_MAGs)
    nb_matches_genomes = len(match_count_genomes)
    if nb_matches_genomes == num_genomes and nb_matches_MAGs >= num_MAGs/2 and is_single:
        is_core = 1
    return(is_core)

# if it is single (in all not just in MAGs)
# and it is present in all MAGs, it is counted
def is_scp_MAG(gene_list, genomes, MAGs):
    is_core_MAG = 0
    is_single = True
    num_genomes = len(genomes)
    num_MAGs = len(MAGs)
    match_count_genomes = dict()
    match_count_MAGs = dict()
    genomes_seen = set()
    for gene in gene_list:
        split_gene = gene.split('_')
        genome_id = "_".join(split_gene[:-1])
        if genome_id in MAGs:
            if genome_id in genomes_seen:
                is_single = False
            else:
                genomes_seen.add(genome_id)
            match_count_MAGs[genome_id] = match_count_MAGs.get(genome_id,0) + 1
    nb_matches_MAGs = len(match_count_MAGs)
    if nb_matches_MAGs == num_MAGs and is_single:
        is_core_MAG = 1
    return(is_core_MAG)

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

#Get genome prefixes and nb of unique genomes from orthofile
orthofile=snakemake.input.ortho_file
genomes_file=snakemake.input.genomes_list
group=snakemake.params.group


genomesAndMAGs = get_g_dict_for_groups(genomes_file)[group]
MAGs = [genome for genome in genomesAndMAGs if "MAG_" in genome]
genomes = [genome for genome in genomesAndMAGs if "MAG_" not in genome]

#Go through each ortholog-family, print to file if single-copy
outfile = snakemake.output.single_ortho
outfile_MAG = snakemake.output.single_ortho_MAGs
with open(orthofile, "r") as fh_orthofile:
    with open(outfile, "w") as fh_outfile:
        with open(outfile_MAG, "w") as fh_outfile_MAG:
            for og in fh_orthofile:
                og = og.strip()
                og_split = og.split(" ")
                og_split.pop(0)
                core_status = is_scp_core(og_split, genomes, MAGs)
                core_mag_status = is_scp_MAG(og_split, genomes, MAGs)
                if core_status == 1:
                    fh_outfile.write(f"{og}\n")
                if core_mag_status == 1:
                    fh_outfile_MAG.write(f"{og}\n")
