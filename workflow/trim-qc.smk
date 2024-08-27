#!/usr/bin/env python

"""
name: trim-qc
description: trimming and qc before and after trimming and trimmed reads concatenated (only concat reads kept) for further analysis
             also includes a rules to map reads to host and then to MAGs in different orders
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - raw_qc:
        + runs fastqc on raw reads
    - trim:
        + trims reads using trimmomatic
    - trim_qc:
        + runs fastqc on trimmed reads
    - concatenate_reads:
        + concatenates trimmed reads
    - re_pair_reads:
        + re-pairs reads (only used for some that htre errors)
    - get_host_unmapd:
        + gets unmapped reads from host mapping
    - host_unmapd_map_to_mags_rep:
        + maps unmapped reads from host mapping to representative MAGs
    - get_host_unmapd_map_mags_unmapped:
        + gets unmapped reads from host mapping to representative MAGs
    - get_unmapd_to_rep_MAGs:
        + gets unmapped reads from representative MAGs mapping
    - unmapd_to_rep_MAGs_host_mapping:
        + maps unmapped reads from representative MAGs mapping to host
    - get_unmapd_to_rep_MAGs_unmaped_host:
        + gets unmapped reads from representative MAGs mapping to host
    - rule unmapd_rep_MAGs_map_all_MAGs_hq:
        + maps unmapped reads from representative MAGs mapping to all MAGs
scripts:
    - write_adapters.py (currently not in use as output hence made is available)
"""

rule raw_qc:
    input:
        reads = lambda wildcards: get_input_file(raw_paths_dict, wildcards.sample, wildcards.read, wildcards.lane)[0],
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


# rest of the pipeline starts from clean reads so commenting out this rule for now
# rule re_pair_reads:
#     input:
#         reads1 = "results/01_trimmedconcatreads/{sample}_R1.fastq.gz",
#         reads2 = "results/01_trimmedconcatreads/{sample}_R2.fastq.gz"
#     output:
#         reads1 = "results/01_cleanreads/{sample}_R1_repaired.fastq.gz",
#         reads2 = "results/01_cleanreads/{sample}_R2_repaired.fastq.gz",
#         singletons = "results/01_cleanreads/{sample}_singletons.fastq.gz"
#     params:
#         java_mem="170", # 70 worked for most but 11/150 samples needed more - should not be more than 85% total memory requested
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         jobname="{sample}_re-paired",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-23:00:00"),
#     resources:
#         mem_mb = convertToMb("400G")
#     threads: 2 # 4 worked for most but 11/150 samples needed more
#     log: "results/01_cleanreads/{sample}_repaired.log"
#     benchmark: "results/01_cleanreads/{sample}_repaired.benchmark"
#     conda: "../config/envs/mapping-env.yaml"
#     shell:
#         """
#         repair.sh -Xmx{params.java_mem}g threads={threads} in1={input.reads1} in2={input.reads2} out1={output.reads1} out2={output.reads2} outs={output.singletons} repair &> {log}
#         """

# map to host and then the unmapped ones to the MAG database
rule get_host_unmapd:
    input:
        bam = "results/03_host_mapping/{sample}_hostfiltered.bam", # bam with unmapped reads from host mapping
    output:
        reads1 = temp("results/04_MapToDBs/{sample}/{sample}_host_unmapd_R1.fastq.gz"),
        reads2 = temp("results/04_MapToDBs/{sample}/{sample}_host_unmapd_R2.fastq.gz"),
        bam_unmapped = temp("results/04_MapToDBs/{sample}/{sample}_host_unmapd.bam"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_get_host_unmapd",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/04_MapToDBs/{sample}/{sample}_get_host_unmapd.log"
    benchmark: "results/04_MapToDBs/{sample}/{sample}_get_host_unmapd.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        samtools view -b {input.bam} | samtools sort -n - > {output.bam_unmapped}
        outreads1={output.reads1}
        outreads2={output.reads2}
        bedtools bamtofastq -i {input.bam} -fq ${{outreads1/.gz/}} -fq2 ${{outreads2/.gz/}} &> {log}
        gzip ${{outreads1/.gz/}}
        gzip ${{outreads2/.gz/}}
        """


rule host_unmapd_map_to_mags_rep:
    input:
        bwa_index = multiext("results/10_instrain/00_prepare_mags/mag_rep_database.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        reads1_hf = "results/04_MapToDBs/{sample}/{sample}_host_unmapd_R1.fastq.gz",
        reads2_hf = "results/04_MapToDBs/{sample}/{sample}_host_unmapd_R2.fastq.gz",
        mag_rep_database = "results/10_instrain/00_prepare_mags/mag_rep_database.fa"
    output:
        bam = temp("results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_MAGs_rep.bam"),
        flagstat_hf = "results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_MAGs_rep.flagstat"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_host_unmapd_map_to_mags_rep",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_to_mags_rep.log"
    benchmark: "results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_to_mags_rep.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa mem -t 4 {input.mag_rep_database} {input.reads1_hf} {input.reads2_hf} | samtools view -bh - | samtools sort - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat_hf}
        """

rule get_host_unmapd_map_mags_unmapped:
    input:
        bam = "results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_MAGs_rep.bam",
    output:
        reads1 = "results/04_MapToDBs/{sample}/{sample}_R1_host_unmapd_map_MAGs_rep_unmapped.fastq.gz",
        reads2 = "results/04_MapToDBs/{sample}/{sample}_R2_host_unmapd_map_MAGs_rep_unmapped.fastq.gz",
        bam_unmapped = temp("results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_MAGs_rep_unmapped.bam"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_get_host_unmapd_map_mags_unmapped",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_mags_rep.log"
    benchmark: "results/04_MapToDBs/{sample}/{sample}_host_unmapd_map_mags_rep.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        samtools view -b -f 4 {input.bam} | samtools sort -n - > {output.bam_unmapped}
        outreads1={output.reads1}
        outreads2={output.reads2}
        bedtools bamtofastq -i {input.bam} -fq ${{outreads1/.gz/}} -fq2 ${{outreads2/.gz/}} &> {log}
        gzip ${{outreads1/.gz/}}
        gzip ${{outreads2/.gz/}}
        """

rule get_unmapd_to_rep_MAGs:
    input:
        bam = "results/10_instrain/01_mapping/{sample}/{sample}.bam",
    output:
        reads1 = "results/04_MapToDBs/{sample}/{sample}_R1_unmapd_rep_MAGs.fastq.gz",
        reads2 = "results/04_MapToDBs/{sample}/{sample}_R2_unmapd_rep_MAGs.fastq.gz",
        bam_unmapped = temp("results/04_MapToDBs/{sample}/{sample}_unmapd_rep_MAGs_unmapped.bam"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_get_unmapd_to_rep_MAGs",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/04_MapToDBs/{sample}/{sample}_get_unmapd_to_rep_MAGs.log"
    benchmark: "results/04_MapToDBs/{sample}/{sample}_get_unmapd_to_rep_MAGs.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        samtools view -b -f 4 {input.bam} | samtools sort -n - > {output.bam_unmapped}
        outreads1={output.reads1}
        outreads2={output.reads2}
        bedtools bamtofastq -i {input.bam} -fq ${{outreads1/.gz/}} -fq2 ${{outreads2/.gz/}} &> {log}
        gzip ${{outreads1/.gz/}}
        gzip ${{outreads2/.gz/}}
        """


rule unmapd_to_rep_MAGs_host_mapping:
    input:
        reads1 = "results/04_MapToDBs/{sample}/{sample}_R1_unmapd_rep_MAGs.fastq.gz",
        reads2 = "results/04_MapToDBs/{sample}/{sample}_R2_unmapd_rep_MAGs.fastq.gz",
        host_db = "data/host_database/apis_bees_db.fasta",
        host_db_index = multiext("data/host_database/apis_bees_db.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam = temp("results/04_MapToDBs/{sample}/{sample}_unmapd_rep_MAGs_host_mapping.bam"),
        flagstat = "results/04_MapToDBs/{sample}/{sample}_unmapd_rep_MAGs_host_mapping.flagstat",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_unmapd_rep_MAGs_host_mapping",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/04_MapToDBs/{sample}/{sample}_unmapd_rep_MAGs_host_mapping.log"
    benchmark: "results/04_MapToDBs/{sample}/{sample}_unmapd_rep_MAGs_host_mapping.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa mem -t {threads} {input.host_db} {input.reads1} {input.reads2} | samtools view -bh - | samtools sort - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule get_unmapd_to_rep_MAGs_unmaped_host:
    input:
        bam = "results/04_MapToDBs/{sample}/{sample}_unmapd_rep_MAGs_host_mapping.bam",
    output:
        reads1 = "results/04_MapToDBs/{sample}/{sample}_R1_unmapd_rep_MAGs_host_unmapd.fastq.gz",
        reads2 = "results/04_MapToDBs/{sample}/{sample}_R2_unmapd_rep_MAGs_host_unmapd.fastq.gz",
        bam_unmapped = temp("results/04_MapToDBs/{sample}/{sample}_unmapd_rep_MAGs_host_unmapd.bam"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_get_unmapd_to_rep_MAGs",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/04_MapToDBs/{sample}/{sample}_get_unmapd_to_rep_MAGs.log"
    benchmark: "results/04_MapToDBs/{sample}/{sample}_get_unmapd_to_rep_MAGs.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        samtools view -b -f 4 {input.bam} | samtools sort -n - > {output.bam_unmapped}
        outreads1={output.reads1}
        outreads2={output.reads2}
        bedtools bamtofastq -i {input.bam} -fq ${{outreads1/.gz/}} -fq2 ${{outreads2/.gz/}} &> {log}
        gzip ${{outreads1/.gz/}}
        gzip ${{outreads2/.gz/}}
        """

# later also check that rep mags recruit as many reads as "all MAGs" together would

# # map reads unmapped to the rep database now to "all MAGs"
# # to see if they are unmapped becuase they are missing from the rep database
# # other better way to check this?
# rule unmapd_rep_MAGs_map_all_MAGs_hq:
#     input:
#         reads1 = "results/04_MapToDBs/{sample}/{sample}_R1_unmapd_rep_MAGs.fastq.gz",
#         reads2 = "results/04_MapToDBs/{sample}/{sample}_R2_unmapd_rep_MAGs.fastq.gz",
#         all_mags