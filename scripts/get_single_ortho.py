#!/usr/bin/env python3
import sys

# adapted from Kirsten's get_singlecp_orthologs.py

def get_genome_prefixes(orthofile):
    fh_orthofile = open(orthofile)
    genome_prefixes = dict()
    for line in fh_orthofile:
        line = line.strip()
        split_line = line.split(" ")
        split_line.pop(0)
        for gene in split_line:
            split_gene = gene.split('_')
            genome_id = split_gene[0]
            genome_prefixes[genome_id] = 1
    fh_orthofile.close()
    return(genome_prefixes)

def is_scp_core(gene_list, group_size):
    is_core = 0
    match_count = dict()
    nb_genes = len(gene_list)
    for gene in gene_list:
        split_gene = gene.split('_')
        genome_id = split_gene[0]
        match_count[genome_id] = match_count.get(genome_id,0) + 1
    nb_matches = len(match_count)
    if nb_matches == group_size and nb_genes == group_size:
        is_core = 1
    return(is_core)

#Get genome prefixes and nb of unique genomes from orthofile
orthofile=snakemake.input.ortho_file
genome_prefixes = get_genome_prefixes(orthofile)
group_size = len(genome_prefixes)

#Go through each ortholog-family, print to file if single-copy
fh_orthofile = open(orthofile)
outfile = snakemake.output.ortho_single
fh_outfile = open(outfile, 'w')
for og in fh_orthofile:
    og = og.strip()
    og_split = og.split(" ")
    og_split.pop(0)
    core_status = is_scp_core(og_split, group_size)
    if core_status == 1:
        fh_outfile.write(f"{og}\n")
fh_orthofile.close()
fh_outfile.close()
