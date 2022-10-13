#!/usr/bin/env python3
import sys
import os


# adapted from Kirsten's get_singlecp_orthologs.py
def is_MAG(genome_name):
    split = genome_name.split("_")
    if len(split) == 3 and split[0] == "MAG":
        return(True)
    else:
        return(False)

# specify if the og is present in only MAGs, only isolates or both
# can be MAGs, Isolates or Both
def present_in(gene_list, genomes, MAGs):
    Isolates_present = True
    if len(genomes) == 0:
        Isolates_present = False
    in_MAGs = False
    in_Isolates = False
    for gene in gene_list:
        split_gene = gene.split('_')
        genome_id = "_".join(split_gene[:-1])
        if genome_id in MAGs:
            in_MAGs = True
        if Isolates_present:
            if genome_id in genomes:
                in_Isolates = True
        if in_Isolates and in_MAGs:
            return("Both")
    if in_Isolates:
        return("Isolates")
    else:
        if in_MAGs:
            return("MAGs")

# specify if the og is present in single copy only in MAGs, only in isolates or both
# can be MAGs, Isolates, Both or NotSingle
def single_copy_in(gene_list, genomes, MAGs):
    Isolates_present = True
    if len(genomes) == 0:
        Isolates_present = False
    single_in_MAGs = True
    if Isolates_present:
        single_in_Isolates = True
    else:
        single_in_Isolates = False
    genomes_seen_MAGs = set()
    if Isolates_present:
        genomes_seen_Isolates = set()
    for gene in gene_list:
        split_gene = gene.split('_')
        genome_id = "_".join(split_gene[:-1])
        if genome_id in MAGs:
            if genome_id in genomes_seen_MAGs:
                single_in_MAGs = False
                continue
            else:
                genomes_seen_MAGs.add(genome_id)
        if Isolates_present:
            if genome_id in genomes:
                if genome_id in genomes_seen_Isolates:
                    single_in_Isolates = False
                    continue
                else:
                    genomes_seen_Isolates.add(genome_id)
    if single_in_MAGs and single_in_Isolates:
        return("Both")
    else:
        if single_in_MAGs:
            return("MAGs")
        if single_in_Isolates:
            return("Isolates")
        return("NotSingle")

# specify if the og is present in single copy and is core (present in all) only in MAGs, only in isolates or both
# can be MAGs, Isolates, Both or NotCore
def core_in(gene_list, genomes, MAGs):
    num_MAGs = len(MAGs)
    num_genomes = len(genomes)
    Isolates_present = True
    if num_genomes == 0:
        Isolates_present = False
    is_core_MAGs = False
    is_core_Isolates = False
    genomes_seen_MAGs = set()
    if Isolates_present:
        genomes_seen_Isolates = set()
    for gene in gene_list:
        split_gene = gene.split('_')
        genome_id = "_".join(split_gene[:-1])
        if genome_id in MAGs:
            if genome_id in genomes_seen_MAGs:
                continue
            else:
                genomes_seen_MAGs.add(genome_id)
        if Isolates_present:
            if genome_id in genomes:
                if genome_id in genomes_seen_Isolates:
                    continue
                else:
                    genomes_seen_Isolates.add(genome_id)
        if len(genomes_seen_MAGs) == num_MAGs:
            is_core_MAGs = True
        if Isolates_present:
            if len(genomes_seen_Isolates) == num_genomes:
                is_core_Isolates = True
    if is_core_MAGs and is_core_Isolates:
        return("Both")
    else:
        if is_core_MAGs:
            return("MAGs")
        else:
            if is_core_Isolates:
                return("Isolates")
            else:
                return("NotCore")

def core_in_half_the_MAGs(gene_list, MAGs):
    num_MAGs = len(MAGs)
    is_core_in_half_the_MAGs = False
    genomes_seen_MAGs = set()
    for gene in gene_list:
        split_gene = gene.split('_')
        genome_id = "_".join(split_gene[:-1])
        if genome_id in MAGs:
            if genome_id in genomes_seen_MAGs:
                continue
            else:
                genomes_seen_MAGs.add(genome_id)
        if len(genomes_seen_MAGs) >= num_MAGs/2:
            is_core_in_half_the_MAGs = True
    return(is_core_in_half_the_MAGs)

def get_og_status(single, core, half_core_MAGs, Isolates_present):
    I_copy = "x"
    I_core = "x"
    M_copy = "x"
    M_core = "x"
    if core == "Both":
        I_core = "Core"
        M_core = "Core"
    if core == "MAGs":
        I_core = "-"
        M_core = "Core"
    if core == "Isolates":
        I_core = "Core"
        M_core = "-"
    if core == "NotCore":
        I_core = "-"
        M_core = "-"
    if single == "Both":
        I_copy = "Scp"
        M_copy = "Scp"
    if single == "MAGs":
        I_copy = "-"
        M_copy = "Scp"
    if single == "Isolates":
        I_copy = "Scp"
        M_copy = "-"
    if single == "NotSingle":
        I_copy = "-"
        M_copy = "-"
    if M_core == "-":
        if half_core_MAGs:
            M_core = "Core0.5"
    if Isolates_present:
        return(f"Isolates:{I_copy}{I_core} ; MAGs:{M_copy}{M_core}")
    else:
        return(f"; MAGs:{M_copy}{M_core}")

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
            cluster = line.split("\t")[11]
            group = line.split("\t")[18]
            # only include groups of interest!
            if group == "g__":
                group = "g__"+cluster
            if group not in g_list_dict.keys():
                g_list_dict[group] = []
            g_list_dict[group].append(genome)
    return(g_list_dict)

#Get genome prefixes and nb of unique genomes from orthofile
orthofile=snakemake.input.ortho_file
orthofile_filt=snakemake.input.ortho_file_filt
perc_id=snakemake.input.perc_id
genomes_file=snakemake.input.genomes_list
group=snakemake.params.group

print('reading single_ortho file')
genome_ids = dict()
OG_fams = dict()
with open(orthofile, "r") as fh_ortho_in:
    for line in fh_ortho_in
        line = line.strip()
        split_line = line.split()
        OG_id = split_line.pop(0)[:-1]
        OG_fams[OG_id] = split_line
        genome_count = genome_count_line(split_line)
        for genome in genome_count:
            genome_ids[genome] = 1

outfile = snakemake.output.summary_orthogroups_filt
with open(orthofile, "r") as fh_orthofile:
    with open(outfile, "w") as fh_outfile:
        fh_outfile.write(f"group,og_name,present_in,single_copy_in,core_in,core_in_half_the_MAGs,status\n")
        for og in fh_orthofile:
            og = fh_orthofile.readline()
            og = og.strip()
            og_name = og.split(":")[0]
            og_split = og.split(" ")
            og_split.pop(0)
            present = present_in(og_split, genomes, MAGs)
            single = single_copy_in(og_split, genomes, MAGs)
            core = core_in(og_split, genomes, MAGs)
            half_core_MAGs = core_in_half_the_MAGs(og_split, MAGs)
            status = get_og_status(single, core, half_core_MAGs, Isolates_present)
            fh_outfile.write(f"{group},{og_name},{present},{single},{core},{half_core_MAGs},{status}\n")
