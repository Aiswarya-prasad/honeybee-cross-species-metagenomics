#!/usr/bin/env python3
import sys
import os

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
# orthofile=os.path.join(os.getcwd(), "database/MAGs_database_Orthofinder/g__Gilliamella/OrthoFinder/g__Gilliamella_single_ortho.txt")
orthofile=snakemake.input.ortho_file
# orthofile=os.path.join(os.getcwd(), "database/MAGs_database_Orthofinder/g__Gilliamella/OrthoFinder/g__Gilliamella_single_ortho_filt.txt")
orthofile_filt=snakemake.input.ortho_file_filt
perc_id=snakemake.input.perc_id
genomes_file=snakemake.input.genomes_list
ref_info=snakemake.input.ref_info
group=snakemake.params.group

outfile = snakemake.output.summary_orthogroups_filt
all_genomes = set(get_g_dict_for_groups(genomes_file)[group])
total_genomes = len(all_genomes)

passed_ogs = set()
with open(orthofile_filt, "r") as fh_orthofile_filt:
    for line in fh_orthofile_filt:
        line = line.strip()
        og_name = line.split(":")[0]
        passed_ogs.add(og_name)

percid_dict = {}
with open(perc_id, "r") as percid_fh:
    for line in percid_fh:
        line = line.strip()
        og_name = line.split("\t")[0]
        # report max perc_id
        percid_dict[og_name] = line.split("\t")[3]

# ref_info = "06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv"
rep_genomes = set()
rep_genome_cluster_dict = dict()
with open(ref_info, "r") as ref_info_fh:
    for line in ref_info_fh:
        line = line.strip()
        split_line = line.split("\t")
        genome_id = split_line[0]
        rep_status = split_line[19]
        cluster = split_line[11]
        if rep_status == "1":
            rep_genomes.add(genome_id)
            rep_genome_cluster_dict[genome_id] = cluster

og_genomes_dict = {}
with open(outfile, "w") as outfile_fh:
    outfile_fh.write(f"group,og_name,present_in,missing_in,perc_id,present_in_ref,passed_filter,missing_genomes,represented_cluster\n")
    with open(orthofile, "r") as fh_orthofile:
        for line in fh_orthofile:
            line = line.strip()
            og_name = line.split(":")[0]
            og_split = line.split(" ")
            og_split.pop(0)
            genomes_covered = set(["_".join(x.split("_")[:-1]) for x in og_split])
            genomes_missing = all_genomes - genomes_covered
            og_genomes_dict[og_name] = genomes_missing
            present_in = len(genomes_covered)
            missing_in = len(genomes_missing)
            if rep_genomes & genomes_covered:
                present_in_ref = "Y"
            else:
                present_in_ref = "N"
            represented_cluster = ";".join([rep_genome_cluster_dict[x] for x in genomes_covered if x in rep_genomes])
            print(f"Checking orthogroup: {og_name}")
            if og_name not in percid_dict.keys():
                print(f"{og_name} missing in perc_id.txt")
                perc_id = "NA"
            else:
                perc_id = percid_dict[og_name]
            passed_filter = "Y" if og_name in passed_ogs else "N"
            missing_genomes = ";".join(genomes_missing)
            outfile_fh.write(f"{group},{og_name},{present_in},{missing_in},{perc_id},{present_in_ref},{passed_filter},{missing_genomes},{represented_cluster}\n")
