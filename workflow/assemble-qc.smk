#!/usr/bin/env python

"""
name: assemble-qc
description: assembly of metagenomes, mapping and other steps for qc
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - assemble_metagenomes
        + assemble sample using default settings
    - map_to_assembly
        + map reads to assemble (to be used for qc of assembly and later to extract gene coverages - so do not delete bam file until then)
    - summarize_mapping_assembly
        + use python script to summarize assmebly length and number of reads mapped for each sample
    - prodigal_get_orfs
        + identify ORFs in all filtered contigs and remove partial and short ORFs
    - dram_annotate_orfs
        + annotate ORFs with dram
scripts:
    - parse_spades_metagenome.py
        + filters scaffolds based on provided length and coverage
    - assembly_summary.py
        + summarizes assemebly length, number of contigs before and after filtering and number of reads mapped for each sample
targets:
    - summary_assembly
"""

# samples that failed in the first run (250):
LARGE_SAMPLES = ["A2-2", "A2-3", "A3-4", "A4-4", "A6-4", "D1-2", "D2-1", "D2-2", "D2-4", "D2-5", "D3-2", "D9-5", "F2-5", "F3-4", "F3-5", "F4-1", "F7-5", "F8-2", "F8-4"]
# samples that failed in the first run (450 - 800) ("A6-4", "D1-2", "D3-2" were completed with the higher RAM):
LARGE_SAMPLES = ["A2-2", "A2-3", "A3-4", "A4-4", "D2-1", "D2-2", "D2-4", "D2-5", "D9-5", "F2-5", "F3-4", "F3-5", "F4-1", "F7-5", "F8-2", "F8-4"]

rule run_bbnorm:
    input:
        reads1 = "results/01_trimmedconcatreads/{sample}_R1.fastq.gz",
        reads2 = "results/01_trimmedconcatreads/{sample}_R2.fastq.gz",
    output:
        reads1 = "results/05_assembly/bbnormed_reads/{sample}_R1.fastq.gz",
        reads2 = "results/05_assembly/bbnormed_reads/{sample}_R2.fastq.gz",
        hist = "results/05_assembly/bbnormed_reads/{sample}.hist",
        peaks = "results/05_assembly/bbnormed_reads/{sample}.peaks",
    params:
        java_mem="3",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_run_bbnorm",
        account="pengel_spirit",
        runtime_s=convertToSec("0-12:00:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/05_assembly/bbnormed_reads/{sample}_bbnorm.log"
    benchmark: "results/05_assembly/bbnormed_reads/{sample}_bbnorm.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bbnorm.sh -Xmx{params.java_mem}g threads={threads} in1={input.reads1} in2={input.reads2} out1={output.reads1} out2={output.reads2} target=40 mindepth=0 hist={output.hist} peaks={output.peaks} &> {log}
        """
        

rule assemble_metagenomes:
    input:
        reads1 = lambda wildcards: f"results/01_trimmedconcatreads/{wildcards.sample}_R1.fastq.gz" if wildcards.sample not in LARGE_SAMPLES else f"results/05_assembly/bbnormed_reads/{wildcards.sample}_R1.fastq.gz",
        reads2 = lambda wildcards: f"results/01_trimmedconcatreads/{wildcards.sample}_R2.fastq.gz" if wildcards.sample not in LARGE_SAMPLES else f"results/05_assembly/bbnormed_reads/{wildcards.sample}_R2.fastq.gz",
    output:
        scaffolds_unparsed = temp("results/05_assembly/all_reads_assemblies/{sample}_scaffolds_unparsed.fasta"),
        scaffolds = "results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta",
        graph = "results/05_assembly/all_reads_assemblies/{sample}_assembly_graph.fastg",
        spades_log = "results/05_assembly/all_reads_assemblies/{sample}_spades.log",
    params:
        outdir = lambda wildcards: "results/05_assembly/all_reads_assemblies/"+wildcards.sample,
        length_t = 1000,
        cov_t = 1,
        memory_limit = lambda wildcards, resources: "800",
        # memory_limit = lambda wildcards, resources: "850" if resources.attempt > 1 else "250",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards, resources: convertToSec("3-00:00:00"),
        # runtime_s=lambda wildcards, resources: convertToSec("1-20:00:00") if resources.attempt > 1 else convertToSec("0-20:00:00"),
    resources:
        attempt = lambda wildcards, attempt: attempt,
        mem_mb = lambda wildcards, attempt: convertToMb("800")
        # mem_mb = lambda wildcards, attempt: convertToMb("800G") if attempt > 1 else convertToMb("200G")
    retries: 0
    # retries: 2
    threads: 4
    log: "results/05_assembly/all_reads_assemblies/{sample}_assemble_metagenomes.log"
    benchmark: "results/05_assembly/all_reads_assemblies/{sample}_assemble_metagenomes.benchmark"
    conda: "../config/envs/spades-env.yaml"
    shell:
        """
        if [ -f \"{params.outdir}/scaffolds.fasta\" ]; then
          echo {params.outdir}/scaffolds.fasta\" exists copying and cleaning\" | tee {log}
        else
          if [ -d \"{params.outdir}/\" ]; then
            echo {params.outdir}\" exists resuming spades\" | tee -a {log}
            spades.py -m {params.memory_limit} -o {params.outdir} --restart-from last || true &>> {log}
          else
            echo \"{params.outdir} not found. Starting new spades run.\" | tee -a {log}
            spades.py -m {params.memory_limit} --meta -1 {input.reads1} -2 {input.reads2} -t {threads} -o {params.outdir} || true &>> {log}
          fi
        fi
        if [ -f \"{params.outdir}/scaffolds.fasta\" ]; then
          cp {params.outdir}/scaffolds.fasta {output.scaffolds_unparsed}
        else
          echo \"assembly may have failed for \"{wildcards.sample} | tee -a {log}
          echo \"touching file \"{output.scaffolds_unparsed} | tee -a {log}
          touch {output.scaffolds_unparsed}
        fi
        python3 scripts/parse_spades_metagenome.py -i {output.scaffolds_unparsed} -o {output.scaffolds} -l {params.length_t} -c {params.cov_t} &>> {log}
        cp {params.outdir}/assembly_graph.fastg {output.graph} &>> {log}
        cp {params.outdir}/spades.log {output.spades_log} &>> {log}
        rm -rf {params.outdir} &>> {log}
        """

rule bamQC:
    input:
        bam = "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}/{sample_assembly}.bam", # only unmapped reads excluded
    output:
        outdir = directory("results/07_MAG_binng_QC/01_backmapping/qualimap_results/{sample_assembly}/")
    params:
        java_mem="300G",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="bamQC",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:30:00"),
    log: "results/07_MAG_binng_QC/01_backmapping/qualimap_results/{sample_assembly}_bamQC.log"
    benchmark: "results/07_MAG_binng_QC/01_backmapping/qualimap_results/{sample_assembly}_bamQC.benchmark"
    threads: 8
    resources:
        mem_mb = 300000,
    conda: "../config/envs/qualimap-env.yaml"
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir {output.outdir} -outformat html --java-mem-size={params.java_mem}
        """

rule assembly_summary:
    input:
        flagstat = expand("results/07_MAG_binng_QC/01_backmapping/{sample}/{sample}_flagstat.txt", sample = SAMPLES_INDIA+SAMPLES_MY),
        scaffolds = expand("results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta" , sample = SAMPLES_INDIA+SAMPLES_MY),
    output:
        outfile = "results/05_assembly/all_reads_assemblies/assembly_summary.txt",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="assembly_summary",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:00:00"),
    resources:
        mem_mb = 8000,
    threads: 4
    log: "results/05_assembly/all_reads_assemblies/assembly_summary.log"
    benchmark: "results/05_assembly/all_reads_assemblies/assembly_summary.benchmark"
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python3 scripts/assembly_summary.py --flagstat {input.flagstat} --scaffolds {input.scaffolds} --outfile {output.outfile} &> {log}
        """