#!/usr/bin/env python

"""
name: mag_phylogenies
description: xxx
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - xxx
scripts:
    - x*.py
targets:
    - ...
"""
# orthofinder and motupan

rule make_phylo_metadata:
    input:
        mag_metadata = "results/09_MAGs_collection/All_mags_sub/All_mags_sub_metadata_summary.tsv",
        Isolate_metadata = "results/11_phylogenies/00_prepare_genomes/Isolate_metadata.tsv",
    output:
        all_metadata = "results/11_phylogenies/00_prepare_genomes/all_genomes_metadata.tsv",

    

rule prepare_isolate_genomes:

rule prepare_faa
rule run_orthofinder_phylo
rule summarise_orthogroups
rule get_single_ortho_phylo
rule extract_orthologs_phylo
rule align_orthologs
rule prune_and_concat
rule make_tree
