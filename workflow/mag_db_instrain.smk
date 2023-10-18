#!/usr/bin/env python

"""
name: mag_db_instrain
description: Takes binning results from metabat2 and summarises it then runs checkm, gtdbtk, drep on mags and creates a filtered mag database for instrain (non-redundant) and redundant for mapping to
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - checkm_evaluation
scripts:
    - *.py
targets:
    - ...
"""

rule make_mag_rep_database:
    input:
        collect_mags_marker = "results/09_MAGs_collection/MAGs/collect_mags.done",
        mag_metadata_summary = lambda wildcards: checkpoints.mag_metadata_summary.get().output.metadata,
        rep_mags = lambda wildcards: expand("results/09_MAGs_collection/MAGs/{mag}.fa", mag = get_rep_mags(checkpoints.mag_metadata_summary.get().output.metadata)),
    output:
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:10:00"),
    threads: 4
    log: "results/10_instrain/00_prepare_mags/mag_rep_database.log"
    benchmark: "results/10_instrain/00_prepare_mags/mag_rep_database.benchmark"
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python scripts/make_mag_rep_database.py \
                --collect_mags_marker {input.collect_mags_marker} \
                --mag_metadata_summary {input.mag_metadata_summary} \
                --mag_rep_database {output.mag_rep_database}
        """

rule get_genes_mag_rep_database:
    input:
        # ffn_files = lambda wildcards: ["results/09_MAGs_collection/prodigal_output/renamed/{mag}/{mag}.ffn" for mag in get_rep_mags(checkpoints.mag_metadata_summary.get().output.metadata)]
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa"
    output:
        instrain_genes_faa = "results/10_instrain/00_prepare_mags/mag_rep_database_genes.faa",
        instrain_genes_fna = "results/10_instrain/00_prepare_mags/mag_rep_database_genes.fna"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:10:00"),
    threads: 4
    log: "results/10_instrain/00_prepare_mags/get_genes_mag_rep_database.log"
    benchmark: "results/10_instrain/00_prepare_mags/get_genes_mag_rep_database.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        prodigal -i {input.mag_rep_database} -d {output.instrain_genes_fna} -a {output.instrain_genes_faa} -p meta &> {log}
        """
        # python3 mag_annotation_extractor.py

rule bwa_index_rep_db:
    input:
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa"
    output:
        bwa_index = multiext("results/10_instrain/00_prepare_mags/mag_rep_database.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:10:00"),
    resources:
        mem_mb = convertToMb("8G")
    threads: 8
    log: "results/10_instrain/00_prepare_mags/bwa_index.log"
    benchmark: "results/10_instrain/00_prepare_mags/bwa_index.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa index {input.mag_rep_database} &> {log}
        """

# rule map_to_rep_MAGs:
#     input:
#         reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
#         reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
#         bwa_index = multiext("results/10_instrain/00_prepare_mags/mag_rep_database.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
#         mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa"
#     output:
#         bam = "results/10_instrain/01_mapping/{sample}/{sample}.bam",
#         flagstat = "results/10_instrain/01_mapping/{sample}/{sample}_flagstat.tsv",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         jobname="{sample}_map_to_rep_MAGs",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-07:10:00"),
#     resources:
#         mem_mb = convertToMb("50G")
#     threads: 4
#     log: "results/10_instrain/01_mapping/{sample}/{sample}_map_to_rep_MAGs.log"
#     benchmark: "results/10_instrain/01_mapping/{sample}/{sample}_map_to_rep_MAGs.benchmark"
#     conda: "../config/envs/mapping-env.yaml"
#     shell:
#         """
#         bwa mem -t {threads} {input.mag_rep_database} {input.reads1} {input.reads2} | samtools view -bh - | samtools sort - > {output.bam}
#         samtools flagstat {output.bam} > {output.flagstat}
#         """

rule bowtie_index:
    input:
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa"
    output:
        bowtie_index = multiext("results/10_instrain/00_prepare_mags/mag_rep_database.fa", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="bowtie_index",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/10_instrain/00_prepare_mags/bowtie_index.log"
    benchmark: "results/10_instrain/00_prepare_mags/bowtie_index.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2-build {input.mag_rep_database} {input.mag_rep_database} &> {log}
        """

rule map_to_rep_MAGs_bowtie2:
    input:
        reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
        reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
        bwa_index = multiext("results/10_instrain/00_prepare_mags/mag_rep_database.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa"
    output:
        bam = "results/10_instrain/01_mapping/{sample}/{sample}_bowtie.bam",
        flagstat = "results/10_instrain/01_mapping/{sample}/{sample}_bowtie_flagstat.tsv",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_map_to_rep_MAGs_bowtie",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/10_instrain/01_mapping/{sample}/{sample}_map_to_rep_MAGs_bowtie.log"
    benchmark: "results/10_instrain/01_mapping/{sample}/{sample}_map_to_rep_MAGs_bowtie.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2 -X 1000 -x {input.mag_rep_database} -1 {input.reads1} -2 {input.reads2} | samtools view -bh - | samtools sort - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """


# A .text file with two columns separated by tabs, 
# where the first column is the name of a scaffold 
# and the second column is the name of the bin / genome the scaffold belongs to.
rule make_scaffold_to_bin_file:
    input:
        rep_mags_db = "results/10_instrain/00_prepare_mags/mag_rep_database.fa",
        rep_mags = lambda wildcards: expand("results/09_MAGs_collection/MAGs/{mag}.fa", mag = get_rep_mags(checkpoints.mag_metadata_summary.get().output.metadata)),
    output:
        scaffold_to_bin_file = "results/10_instrain/00_prepare_mags/scaffold_to_bin_file.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_scaffold_to_bin_file",
        account="pengel_spirit",
        runtime_s=convertToSec("0-01:10:00"),
    resources:
        mem_mb = convertToMb("8G")
    threads: 4
    log: "results/10_instrain/00_prepare_mags/make_scaffold_to_bin_file.log"
    benchmark: "results/10_instrain/00_prepare_mags/make_scaffold_to_bin_file.benchmark"
    run:
        with open(output.scaffold_to_bin_file, "w") as f:
            for mag in input.rep_mags:
                with open(mag, "r") as m:
                    for line in m:
                        mag_name = os.path.basename(mag).split(".")[0]
                        if line.startswith(">"):
                            scaffold = line.strip().split(">")[1]
                            f.write(f"{scaffold}\t{mag_name}\n")

rule instrain_profile:
    input:
        bam = "results/10_instrain/01_mapping/{sample}/{sample}_bowtie.bam",
        # bam = "results/10_instrain/01_mapping/{sample}/{sample}.bam",
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa",
        scaffold_to_bin_file = "results/10_instrain/00_prepare_mags/scaffold_to_bin_file.tsv",
        instrain_genes_file = "results/10_instrain/00_prepare_mags/mag_rep_database_genes.fna"
    output:
        # gene_info = "results/10_instrain/02_instrain_profile/{sample}/output/{sample}_gene_info.tsv",
        # scaffold_info = "results/10_instrain/02_instrain_profile/{sample}/output/{sample}_scaffold_info.tsv",
        # mapping_info = "results/10_instrain/02_instrain_profile/{sample}/output/{sample}_mapping_info.tsv",
        # genome_info = "results/10_instrain/02_instrain_profile/{sample}/output/{sample}_genome_info.tsv",
        # SNVs = "results/10_instrain/02_instrain_profile/{sample}/output/{sample}_SNVs.tsv.gz",
        marker = touch("results/10_instrain/02_instrain_profile/{sample}_profile.done/")
    params:
        outdir = "results/10_instrain/02_instrain_profile/{sample}",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_instrain_profile",
        account="pengel_spirit",
        runtime_s=convertToSec("0-22:10:00"),
    resources:
        mem_mb = convertToMb("250G")
    threads: 8
    conda: "../config/envs/instrain_env.yaml"
    log: "results/10_instrain/02_instrain_profile/{sample}_instrain_profile.log"
    benchmark: "results/10_instrain/02_instrain_profile/{sample}_instrain_profile.benchmark"
    shell:
        """
        inStrain profile {input.bam} {input.mag_rep_database} -o {params.outdir} \
                -p {threads} -g {input.instrain_genes_file} \
                --max_insert_relative 5 \
                -s {input.scaffold_to_bin_file}
        touch {output.marker}
        """

rule instrain_profile_plot:
    input:
        marker = "results/10_instrain/02_instrain_profile/{sample}_profile.done/"
    output:
        done = touch("results/10_instrain/04_instrain_plot_marker/{sample}_profile_plots.done")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_instrain_profile_plot",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("100G")
    threads: 16
    log: "results/10_instrain/04_instrain_plot_marker/{sample}_profile_plots.log"
    benchmark: "results/10_instrain/04_instrain_plot_marker/{sample}_profile_plots.benchmark"
    conda: "../config/envs/instrain_env.yaml"
    shell:
        """
        profile={input.marker}
        profile=${{profile/_profile.done/}}
        inStrain plot -i ${{profile}} -pl a -p {threads}
        touch {output.done}
        """

# Do instrain compare for all samples but one genome at a time! - update the code
# and run it later. Also needs profile to be done in database mode.. So maybe new rule for that

# get_rep_mags

rule instrain_profile_db_mode:
    input:
        bam = "results/10_instrain/01_mapping/{sample}/{sample}_bowtie.bam",
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa",
        scaffold_to_bin_file = "results/10_instrain/00_prepare_mags/scaffold_to_bin_file.tsv",
        instrain_genes_file = "results/10_instrain/00_prepare_mags/mag_rep_database_genes.fna",
    output:
        marker = touch("results/10_instrain/02_instrain_profile_db_mode/{sample}_profile_db_mode.done/")
    params:
        outdir = "results/10_instrain/02_instrain_profile_db_mode/{sample}",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_instrain_profile_db_mode",
        account="pengel_spirit",
        runtime_s=convertToSec("0-22:10:00"),
    resources:
        mem_mb = convertToMb("250G")
    threads: 8
    conda: "../config/envs/instrain_env.yaml"
    log: "results/10_instrain/02_instrain_profile_db_mode/{sample}_instrain_profile_db_mode.log"
    benchmark: "results/10_instrain/02_instrain_profile_db_mode/{sample}_instrain_profile_db_mode.benchmark"
    shell:
        """
        inStrain profile {input.bam} {input.mag_rep_database} -o {params.outdir} \
                -p {threads} -g {input.instrain_genes_file} \
                --max_insert_relative 5 --database_mode \
                -s {input.scaffold_to_bin_file}
        touch {output.marker}
        """

rule instrain_compare:
    input:
        markers = expand("results/10_instrain/02_instrain_profile_db_mode/{sample}_profile_db_mode.done/",  sample=SAMPLES_sub),
        # repeat with all samples of interest later
        scaffold_to_bin_file = "results/10_instrain/00_prepare_mags/scaffold_to_bin_file.tsv"
    output:
        compare_marker = touch("results/10_instrain/03_instrain_compare/log_files/{mag}_compared.done"),
    params:
        outdir = lambda wildcards: f"results/10_instrain/03_instrain_compare/{get_species_from_rep_mag(wildcards.mag)}",
        profiles = [f"results/10_instrain/02_instrain_profile_db_mode/{sample}" for sample in SAMPLES],
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00"),
    resources:
        mem_mb = convertToMb("200G")
    conda: "../config/envs/instrain_env.yaml"
    log: "results/10_instrain/03_instrain_compare/log_files/instrain_compare_{mag}.log"
    benchmark: "results/10_instrain/03_instrain_compare/log_files/instrain_compare_{mag}.benchmark"
    # log: lambda wildcards: f"results/10_instrain/03_instrain_compare/log_files/instrain_compare_{get_species_from_rep_mag(wildcards.mag)}.log"
    # benchmark: lambda wildcards: f"results/10_instrain/03_instrain_compare/log_files/instrain_compare_{get_species_from_rep_mag(wildcards.mag)}.benchmark"
    threads: 8
    shell:
        """
        inStrain compare -i {params.profiles} -s {input.scaffold_to_bin_file} \
                -p {threads} -o {params.outdir} \
                --database_mode --genome {wildcards.mag} || touch {output.compare_marker}.singleton
        touch {output.compare_marker}
        """

rule aggregate_compare:
    input:
        markers = lambda wildcards: [f"results/10_instrain/03_instrain_compare/log_files/{mag}_compared.done" for mag in get_rep_mags(checkpoints.mag_metadata_summary.get().output.metadata)]
    output:
        aggregate_marker = touch("results/10_instrain/03_instrain_compare/all_compared.done"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("2G")
    log: "results/10_instrain/03_instrain_compare/aggregate_compare.log"
    benchmark: "results/10_instrain/03_instrain_compare/aggregate_compare.benchmark"