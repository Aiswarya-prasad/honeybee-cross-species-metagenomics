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

checkpoint make_phylo_metadata:
    input:
        mag_metadata = lambda wildcards: checkpoints.mag_metadata_summary.get().output.metadata,
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
        mem_mb = convertToMb("16G")
    threads: 2
    log: "results/11_phylogenies/phylo_genomes_metadata.log"
    # benchmark:
    # conda: 
    shell:
        """
        python3 scripts/make_phylo_metadata.py \
            --mag_metadata {input.mag_metadata} \
            --Isolate_metadata {input.Isolate_metadata} \
            --phylo_metadata {output.phylo_metadata}
        """

rule download_assembly_summary:
    output:
        assembly_summary = "results/11_phylogenies/assembly_summary.txt"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="download_assembly_summary",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
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
        genome_file_faa = lambda wildcards: get_genome_path(wildcards.genome, checkpoints.mag_metadata_summary.get().output.metadata, "faa")["path"],
        genome_file_ffn = lambda wildcards: get_genome_path(wildcards.genome, checkpoints.mag_metadata_summary.get().output.metadata, "ffn")["path"],
    output:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.faa",
        genome_ffn = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.ffn",
    params:
        input_type = lambda wildcards: get_genome_path(wildcards.genome, checkpoints.mag_metadata_summary.get().output.metadata, "faa")["type"],
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{genome}_prepare_genomes_faa",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 2
    log: "results/11_phylogenies/00_genomes/{genome}/{genome}_prepare_genomes.log"
    benchmark: "results/11_phylogenies/00_genomes/{genome}/{genome}_prepare_genomes.benchmark.json"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        genome_ffn={output.genome_ffn}
        genome_gff=${{genome_ffn/.ffn/.gff}}
        if [ {params.input_type} == "fna" ]
        then
            echo "running prodigal as fna was returned"
            # remember genoe_file_faa is actually fna here!
            prodigal -i {input.genome_file_faa} -a {output.genome_faa} -o ${{genome_gff}} -f gff
            bedtools getfasta -fi {input.genome_file_faa} -bed ${{genome_gff}} -fo ${{genome_ffn}}
        else
            cp {input.genome_file_faa} {output.genome_faa}
            echo "copied faa file to final destination"
            cp {input.genome_file_ffn} {output.genome_ffn}
            echo "copied ffn file to final destination"
        fi
        """

rule rename_faa_and_ffn:
    input:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.faa",
        genome_ffn = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.ffn",
    output:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}.faa",
        genome_ffn = "results/11_phylogenies/00_genomes/{genome}/{genome}.ffn",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_phylo_metadata",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 4
    log: "results/11_phylogenies/00_genomes/{genome}/{genome}_rename_faa_and_ffn.log"
    run:
        i = 1
        with open(output.genome_faa, "w") as fh,  open(f"{input.genome_faa}") as f:
            for line in f:
                if line.startswith(">"):
                    fh.write(f">{wildcards.genome}_{i}\n")
                    i += 1
                else:
                    fh.write(line)
        j = 1
        with open(output.genome_ffn, "w") as fh,  open(f"{input.genome_ffn}") as f:
            for line in f:
                if line.startswith(">"):
                    fh.write(f">{wildcards.genome}_{j}\n")
                    j += 1
                else:
                    fh.write(line)
        if i != j:
            # consider adding additional checks to make sure the faa and fna record names are correctly matched
            sys.exit("Number of genes in faa and ffn files are not equal for {wildcards.genome}")


rule collect_faa_orthofinder:
    input:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}.faa"
    output:
        genome_faa = "results/11_phylogenies/01_orthofinder_input/{genus}/{genome}.faa",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 2
    log: "results/11_phylogenies/01_orthofinder_input/{genus}_collect_input_logs/{genome}_collect.log"
    benchmark: "results/11_phylogenies/01_orthofinder_input/{genus}_collect_input_logs/{genome}_collect.benchmark"
    run:
        shell("cp {input.genome_faa} {output.genome_faa}")


rule run_orthofinder:
    input:
        faa_files = lambda wildcards: expand("results/11_phylogenies/01_orthofinder_input/{{genus}}/{genome}.faa", genome=get_mags_for_genus_phylogeny(wildcards.genus, checkpoints.make_phylo_metadata.get().output.phylo_metadata)),
    output:
        trees = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Species_Tree/SpeciesTree_rooted.txt",
        orthofile = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Orthogroups/Orthogroups.txt"
    params:
        outdir = "results/11_phylogenies/02_orthofinder_results/{genus}",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("1-2:10:00") if wildcards.genus == "g__Lactobacillus" or wildcards.genus in ["g__Bifidobacterium", "g__Lactobacillus"] else convertToSec("0-15:10:00"),
    resources:
        mem_mb = convertToMb("200G")
    threads: 16
    log: "results/11_phylogenies/02_orthofinder_results/{genus}_run_orthofinder.log"
    benchmark: "results/11_phylogenies/02_orthofinder_results/{genus}_run_orthofinder.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        # https://github.com/davidemms/OrthoFinder#running-orthofinder
        rm -rf {params.outdir}
        orthofinder -ot -t {threads} -n {wildcards.genus} \
        -f $(dirname {input.faa_files[0]}) -o {params.outdir} \
        -M msa -T fasttree \
        2>&1 | tee -a {log}
        """

rule run_orthofinder_iqtree:
    input:
        trees = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Species_Tree/SpeciesTree_rooted.txt",
        orthofile = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Orthogroups/Orthogroups.txt"
    output:
        trees = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}_iqtree/Species_Tree/SpeciesTree_rooted.txt"
    params:
        prev_outdir = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}",
        new_outdir = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}_iqtree",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("1-2:10:00") if wildcards.genus == "g__Lactobacillus" or wildcards.genus in ["g__Bifidobacterium", "g__Lactobacillus"] else convertToSec("0-15:10:00"),
    resources:
        mem_mb = convertToMb("200G")
    threads: 16
    log: "results/11_phylogenies/02_orthofinder_results/{genus}_iqtree_run_orthofinder.log"
    benchmark: "results/11_phylogenies/02_orthofinder_results/{genus}_iqtree_run_orthofinder.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        rm -rf {params.new_outdir}
        orthofinder -fg {params.prev_outdir} \
            -M msa -T iqtree -ot -t {threads} \
            -n {wildcards.genus}_iqtree 2>&1 | tee {log}
        """
