#!/usr/bin/env python

"""
name: trim-qc
description: trimming and qc before and after
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
    - trimmed_files
    - qc_raw = expand("fastqc/raw/{sample}_{read}_fastqc.html", sample=SAMPLES+SAMPLES_KE, read=config["READS"]),
    - qc_trimmed = expand("fastqc/trim/{sample}_{read}_trim_fastqc.html", sample=SAMPLES+SAMPLES_KE, read=config["READS"]),
"""

rule raw_qc:
    input:
        reads=ancient(os.path.join("00_RawData", "{sample}_{read}.fastq.gz")),
    output:
        html="fastqc/raw/{sample}_{read}_fastqc.html",
        zip="fastqc/raw/{sample}_{read}_fastqc.zip"
    params:
        outdir="fastqc/raw/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname=lambda wildcards: wildcards.sample+"_"+wildcards.read+"_qc",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "fastqc/raw/{sample}_{read}_qc.log"
    benchmark: "fastqc/raw/{sample}_{read}_qc.benchmark"
    threads: 2
    conda: "envs/trim-qc-env.yaml"
    shell:
        """
        fastqc -t {threads} {input.reads} -o {params.outdir} &> {log}
        """

# """
# rule make_adapters:
#     input:
#         "config/index_table.csv"
#     output:
#         "config/Adapters-PE.fa"
#     script:
#         "scripts/write_adapters.py"
# """

rule trim:
    input:
        reads1=os.path.join("00_RawData", "{sample}_R1.fastq.gz"),
        reads2=os.path.join("00_RawData", "{sample}_R2.fastq.gz"),
        adapter="config/Adapters-PE.fa"
        # adapter=rules.make_adapters.output
    output:
        reads1 = "01_Trimmed/{sample}_R1_trim.fastq.gz",
        reads2 = "01_Trimmed/{sample}_R2_trim.fastq.gz",
        reads1_unpaired = temp("01_Trimmed/{sample}_R1.unpaired.fastq.gz"),
        reads2_unpaired = temp("01_Trimmed/{sample}_R2.unpaired.fastq.gz")
    params:
        # add adapter definition to config later
        adapters = lambda wildcards: ADAPTERS["SAMPLES_INDIA"] if wildcards.sample in SAMPLES_INDIA else (ADAPTERS["SAMPLES_KE"] if wildcards.sample in SAMPLES_KE else ADAPTERS["default"]),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_trim",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "01_Trimmed/{sample}_trim.log"
    benchmark: "01_Trimmed/{sample}_trim.benchmark"
    conda: "envs/trim-qc-env.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.reads1} {input.reads2} {output.reads1} {output.reads1_unpaired} {output.reads2} {output.reads2_unpaired} ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:28 TRAILING:28  MINLEN:60 &> {log}
        """

rule trim_qc:
    input:
        reads="01_Trimmed/{sample}_{read}_trim.fastq.gz",
    output:
        html="fastqc/trim/{sample}_{read}_trim_fastqc.html",
        zip="fastqc/trim/{sample}_{read}_trim_fastqc.zip"
    threads: 2
    conda: "envs/trim-qc-env.yaml"
    params:
        outdir="fastqc/trim/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_{read}_trim_qc",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "fastqc/trim/{sample}_{read}_trim_qc.log"
    benchmark: "fastqc/trim/{sample}_{read}_trim_qc.benchmark"
    resources:
        mem_mb = 8000
    shell:
        """
        fastqc -t {threads} {input.reads} -o {params.outdir} &> {log}
        """