#!/usr/bin/env python

"""
name: backmapping-binning
description: takes the assembled scaffolds and performs backmapping followed by writing depth files and finally uses those for binning by metabat2
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - rename_scaffolds
    - build_bwa_index
    - backmapping
    - bamQC
scripts:
    - rename_scaffolds.py
targets:
    - build_bwa_index
"""

rule rename_scaffolds:
    input:
        scaffolds = rules.assemble_metagenomes.output.scaffolds,
    output:
        scaffolds = "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample_assembly}_scaffolds.fasta",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="build_bwa_index",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample_assembly}_rename_scaffolds.log"
    benchmark: "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample_assembly}_rename_scaffolds.benchmark"
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python3 scripts/rename_scaffolds.py --scaffolds_in {input.scaffolds} --scaffolds_out {output.scaffolds} --sample {wildcards.sample} &>> {log}
        """

rule build_bwa_index:
    input:
        scaffolds = "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample_assembly}_scaffolds.fasta"
    output:
        bwa_index = multiext("results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample_assembly}_scaffolds.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="build_bwa_index",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample_assembly}_build_bwa_index.log"
    benchmark: "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample_assembly}_build_bwa_index.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa index {input.scaffolds} &>> {log}
        """

rule backmapping:
    input:
        scaffolds = lambda wildcards: f"results/07_MAG_binng_QC/00_assembled_scaffolds/{wildcards.sample_assembly}/{wildcards.sample_assembly}_scaffolds.fasta",
        bwa_index = lambda wildcards: multiext(f"results/07_MAG_binng_QC/00_assembled_scaffolds/{wildcards.sample_assembly}/{wildcards.sample_assembly}_scaffolds.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
        reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
        singletons = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_singletons.fastq.gz"
    output:
        flagstat = "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}/{sample}_mapped_flagstat.txt",
        bam = temp("results/07_MAG_binng_QC/01_backmapping/{sample_assembly}/{sample}.bam") # only unmapped reads excluded
    params:
        match_length = 50,
        edit_distance = 5, # methods in microbiomics recommends 95 perc identity
        # since reads are 150 bp long, 5 mismatches is 3.3% mismatch which is almost as instrain recommends
        # even less chances of strains mismapping
        filter_script = "scripts/filter_bam.py",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample_assembly}_{sample}_backmapping",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}/{sample}_backmapping.log"
    benchmark: "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}/{sample}_backmapping.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa mem -a -t {threads} {input.scaffolds} {input.reads1} {input.reads2} \
        | samtools view -F 4 -h - |  python3 {params.filter_script} -e 5 -m 50 | samtools sort -O bam -@ {threads} > {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        """

rule make_depthfile:
    input:
        bams = expand("results/07_MAG_binng_QC/01_backmapping/{{sample_assembly}}/{sample}.bam", sample=SAMPLES_INDIA+SAMPLES_MY)
    output:
        depthfile = "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}_depths/{sample_assembly}_depthfile.txt"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample_assembly}_make_depthfile",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}_depths/{sample_assembly}_make_depthfile.log"
    benchmark: "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}_depths/{sample_assembly}_make_depthfile.benchmark"
    conda: "../config/envs/mags-env.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths –outputDepth {output.depthfile} {input.bams}
        """

rule run_metabat2:
    input:
        scaffolds = lambda wildcards: f"results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample}_scaffolds.fasta",
        depthfile = "results/07_MAG_binng_QC/01_backmapping/{sample}_depths/{sample}_depthfile.txt",
    output:
        bins = directory("results/07_MAG_binng_QC/02_bins/{sample}/")
    params:
        clustersize = 200000, # 200kb metabat2 default
        maxEdges = 500, # 500 metabat2 default
        mincontiglen = 2000, # 2000 metabat2 default
        minCV = 1, # 1 metabat2 default
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_make_depthfile",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/07_MAG_binng_QC/02_bins/{sample}/{sample}_run_metabat2.log"
    benchmark: "results/07_MAG_binng_QC/02_bins/{sample}/{sample}_run_metabat2.benchmark"
    conda: "../config/envs/mags-env.yaml"
    shell:
        """
        metabat2 -i {input.scaffolds} -a {input.depthfile} \
        -o {output.bins}/{wildcards.sample}_ --minContig {params.mincontiglen} \
        --maxEdges {params.maxEdges} -x {params.minCV} --minClsSize {params.clustersize} --saveCls -v
        --unbinned
        """

# write rule to summarize contig fates and stats for each bin and assembly (check if this is the right way to do it) 
# rule summarize_contig_fates:
#     input:
#         bins = directory("results/07_MAG_binng_QC/02_bins/{sample}/") 
#     output:
#         contig_fates = "results/07_MAG_binng_QC/02_bins/{sample}_contig_fates/{sample}_contig_fates.txt",
#         contig_stats = "results/07_MAG_binng_QC/02_bins/{sample}_contig_stats/{sample}_contig_stats.txt"
#     params:
#         mailto="

# rules for checkm, gtdbtk, and drep

# summarize magOTUs and rename headers as needed

# annotate mags

# orthofinder and motupan

# make non-redundant mag database for instrain
# ! compare corecov and instrain results from the previous analysis !
# make redundant magdatabase for core coverage estimation

