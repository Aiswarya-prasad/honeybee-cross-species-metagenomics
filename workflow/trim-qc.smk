#!/usr/bin/env python

"""
name: trim-qc
description: trimming and qc before and after trimming and trimmed reads concatenated (only concat reads kept)
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - raw_qc
        + runs fastqc on raw reads
    - make_adapters
        + custom script used when adaptor sequences are not standard (currently not in use)
    - trim
        + runs trimmomatic on raw reads
    - trim_qc
        + runs fastqc on trimmed reads
scripts:
    - write_adapters.py (currently not in use)
targets:
    - trimmed files
    - concatenated files after trimming
    - raw and trimmmed file fastqc results
"""

rule raw_qc:
    input:
        reads = lambda wildcards: get_input_file(raw_paths_dict, wildcards.sample, wildcards.read, wildcards.lane, wildcards.run)[0],
    output:
        html = "results/00_rawreads/fastqc/{sample}_{lane}_{read}_{run}_fastqc.html",
        zip = "results/00_rawreads/fastqc/{sample}_{lane}_{read}_{run}_fastqc.zip"
    params:
        outdir = "results/00_rawreads/fastqc/",
        outname = lambda wildcards: f"{wildcards.sample}_{wildcards.lane}_{wildcards.read}_{wildcards.run}",
        mailto = "aiswarya.prasad@unil.ch",
        mailtype = "BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname = lambda wildcards: f"{wildcards.sample}_{wildcards.read}_qc",
        account = "pengel_spirit",
        runtime_s = convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "results/00_rawreads/fastqc/{sample}_{lane}_{read}_{run}_qc.log"
    benchmark: "results/00_rawreads/fastqc/{sample}_{lane}_{read}_{run}_qc.benchmark"
    threads: 2
    conda: "../config/envs/trim-qc-env.yaml"
    shell:
        """
        zcat {input.reads} | fastqc -t {threads} -o {params.outdir} stdin:{params.outname}
        """

"""
# rule make_adapters:
#     input:
#         "../config/index_table.csv"
#     output:
#         "../config/Adapters-PE.fa"
#     script:
#         "../scripts/write_adapters.py"
"""

rule trim:
    input:
        reads1 = lambda wildcards: get_input_file(raw_paths_dict, wildcards.sample, "R1", wildcards.lane, wildcards.run)[0],
        reads2 = lambda wildcards: get_input_file(raw_paths_dict, wildcards.sample, "R2", wildcards.lane, wildcards.run)[0],
        adapter = "config/Adapters-PE.fa"
        # adapter=rules.make_adapters.output
    output:
        reads1 = "results/00_trimmedreads/{sample}_{lane}_R1_{run}_trim.fastq.gz",
        reads2 = "results/00_trimmedreads/{sample}_{lane}_R2_{run}_trim.fastq.gz",
        reads1_unpaired = temp("results/00_trimmedreads/{sample}_{lane}_R1_{run}_unpaired.fastq.gz"),
        reads2_unpaired = temp("results/00_trimmedreads/{sample}_{lane}_R2_{run}_unpaired.fastq.gz")
    params:
        adapters = lambda wildcards: ADAPTERS["SAMPLES_INDIA"] if wildcards.sample in SAMPLES_INDIA else (ADAPTERS["SAMPLES_KE"] if wildcards.sample in SAMPLES_KE else ADAPTERS["default"]),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_{lane}_{run}_trim",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/00_trimmedreads/{sample}_{lane}_{run}_trim.log"
    benchmark: "results/00_trimmedreads/{sample}_{lane}_{run}_trim.benchmark"
    conda: "../config/envs/trim-qc-env.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.reads1} {input.reads2} {output.reads1} {output.reads1_unpaired} {output.reads2} {output.reads2_unpaired} ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:28 TRAILING:28  MINLEN:60 &> {log}
        """

rule trim_qc:
    input:
        reads = "results/00_trimmedreads/{sample}_{lane}_{read}_{run}_trim.fastq.gz",
    output:
        html="results/00_trimmedreads/fastqc/{sample}_{lane}_{read}_{run}_trim_fastqc.html",
        zip="results/00_trimmedreads/fastqc/{sample}_{lane}_{read}_{run}_trim_fastqc.zip"
    threads: 2
    conda: "../config/envs/trim-qc-env.yaml"
    params:
        outdir="results/00_trimmedreads/fastqc/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_{lane}_{read}_{run}_trim_qc",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "results/00_trimmedreads/fastqc/{sample}_{lane}_{read}_{run}_trim_qc.log"
    benchmark: "results/00_trimmedreads/fastqc/{sample}_{lane}_{read}_{run}_trim_qc.benchmark"
    resources:
        mem_mb = 8000
    shell:
        """
        fastqc -t {threads} {input.reads} -o {params.outdir} &> {log}
        """

rule concatenate_reads:
    input:
        reads = lambda wildcards: [f"results/00_trimmedreads/{x.split('.fastq.gz')[0]}_trim.fastq.gz" for x in get_list_of_values(get_renamed_input_files({key: raw_paths_dict[key] for key in raw_paths_dict.keys() if key == wildcards.sample})) if f"_{wildcards.read}" in x]
    output:
        concat_reads = "results/01_trimmedconcatreads/{sample}_{read}.fastq.gz"
    threads: 2
    conda: "../config/envs/trim-qc-env.yaml"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_{read}_concat",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "results/01_trimmedconcatreads/{sample}_{read}_concat.log"
    benchmark: "results/01_trimmedconcatreads/{sample}_{read}_concat.benchmark"
    resources:
        mem_mb = 8000
    shell:
        """
        cat {input.reads} > {output.concat_reads}
        """