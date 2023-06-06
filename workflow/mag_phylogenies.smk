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
        mag_metadata = lambda wildcards: checkpoints.mag_metadata.get().output.mag_metadata,
        Isolate_metadata = "config/Outgroup_isolate_genomes.tsv",
    output:
        phylo_metadata = "results/11_phylogenies/phylo_genomes_metadata.tsv",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_phylo_metadata",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("8G")
    threads: 2
    log: "results/11_phylogenies/phylo_genomes_metadata.log"
    # benchmark:
    # conda: 
    shell:
        """
        python3 make_phylo_metadata.py \
            --mag_metadata {input.mag_metadata} \
            --Isolate_metadata {input.Isolate_metadata} \
            --phylo_metadata {output.phylo_metadata}
        """

def get_genome_path(genome, mag_summary):
    with open(metadata, "r") as f:
        header = f.readline()
        ID_ind = header.split("\t").index("ID")
        IDs = [line.split("\t")[ID_ind] for line in f.readlines()]
        if genome not in IDs:
            print(f"{genome} not found in {metadata}, will be downloaded from NCBI")   
            with open("results/11_phylogenies/assembly_summary.txt", 'r') as f:
        for line in f:
            if line.startswith('#') and 'ftp_path' in line:
                ind_link = line.strip().split('\t').index('ftp_path')
            else:
                line_split = line.strip().split('\t')
                if line_split[0] == genome:
                    ftp_link = line_split[ind_link]
                    ftp_id=''.join(ftp_link.split("/")[-1])
                    fna_path = f'{ftp_link}/{ftp_id}_genomic.fna.gz'
                    faa_path = f'{ftp_link}/{ftp_id}_protein.faa.gz'
                    type = "link"
        else:
            fna_path = f"results/09_MAGs_collection/All_mags_sub/MAGs/{genome}.fa"
            faa_path = f"results/09_MAGs_collection/All_mags_sub/prodigal_output/from_checkm/{genome}.faa"
            type = "path"
        genome_path_dict = {"fna": fna_path, "faa": faa_path, "type": type}
        return genome_path_dict
            

rule download_assembly_summary:
    output:
        assembly_summary = "results/11_phylogenies/assembly_summary.txt"
    params:
        input_type = get_genome_path(wildcards.genome, checkpoints.mag_summary.get().output.mag_summary)["type"],
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="download_assembly_summary",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("8G")
    threads: 2
    log: "results/11_phylogenies/download_assembly_summary.log"
    shell:
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -O {output.assembly_summary}
        """

rule prepare_genomes_faa:
    input:
        phylo_genomes_metadata = "results/11_phylogenies/phylo_genomes_metadata.tsv",
        assembly_summary = "results/11_phylogenies/assembly_summary.txt",
        genome_faa = get_genome_path(wildcards.genome, checkpoints.mag_summary.get().output.mag_summary)["faa"],
    output:
        genome_faa = lambda wildcards: "results/11_phylogenies/00_genomes_faa/{genome}.faa",
    params:
        input_type = get_genome_path(wildcards.genome, checkpoints.mag_summary.get().output.mag_summary)["type"],
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{genome}_prepare_genomes_faa",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("8G")
    threads: 2
    log: "results/11_phylogenies/00_genomes_faa/{genome}_prepare_genomes.log"
    benchmark: "results/11_phylogenies/00_genomes_faa/{genome}_prepare_genomes.benchmark.json"
    # conda: 
    run:
        if {params.input_type} == "link":
            shell("wget {input.genome_faa} -O {output.genome_faa}.gz")
            shell("gunzip {output.genome_faa}.gz")
        else:
            shell("cp {input.genome_faa} {output.genome_faa}")

rule prepare_faa_orthofinder:
    input:
        genome_faa = "results/11_phylogenies/00_genomes_faa/{genome}.faa"
    output:
        genome_faa = "results/11_phylogenies/01_orthofinder/{group}/{genome}.faa",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_phylo_metadata",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("8G")
    threads: 2
    log: "results/11_phylogenies/01_orthofinder/prepare_faa_orthofinder.log"
    run:
        i = 1
        with open({output.genome_faa}, "w") as fh:
            with open(f"{input.genome_faa}", "r") as fh:
                for line in fh:
                    if line.startswith(">"):
                        old_header = line.split(">")[1]
                        new_header = f">{wildcards.genome}_{(four_digit(i))} {old_header}"
                        i += 1
                        fh.write(new_header)
                    else:
                        fh.write(line)


rule run_orthofinder_phylo
rule summarise_orthogroups
rule get_single_ortho_phylo
rule extract_orthologs_phylo
rule align_orthologs
rule prune_and_concat
rule make_tree
