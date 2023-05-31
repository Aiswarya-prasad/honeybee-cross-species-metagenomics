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
# LARGE_SAMPLES = ["A2-2", "A2-3", "A3-4", "A4-4", "D2-1", "D2-2", "D2-4", "D2-5", "D9-5", "F2-5", "F3-4", "F3-5", "F4-1", "F7-5", "F8-2", "F8-4"]
SAMPLES_TO_HOST_FILTER = ["D2-2", "A2-2", "F2-5"]

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

rule bwa_index_db:
    input:
        host_db = "data/{db_name}/{db}.fasta",
    output:
        host_db_index = multiext("data/{db_name}/{db}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{db_name}_{db}_bwa_index_db",
        account="pengel_spirit",
        runtime_s=convertToSec("0-01:00:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/03_host_mapping/{db_name}_{db}_bwa_index_db.log"
    benchmark: "results/03_host_mapping/{db_name}_{db}_bwa_index_db.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa index {input.host_db} &> {log}
        """


rule map_reads_to_host:
    input:
        reads1 = "results/01_cleanreads/{sample}_R1_repaired.fastq.gz",
        reads2 = "results/01_cleanreads/{sample}_R2_repaired.fastq.gz",
        host_db = "data/host_database/apis_bees_db.fasta",
        host_db_index = multiext("data/host_database/apis_bees_db.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam = temp("results/03_host_mapping/{sample}.bam"),
        bam_hostfiltered = temp("results/03_host_mapping/{sample}_hostfiltered.bam"),
        flagstat = "results/03_host_mapping/{sample}_flagstat.txt",
        coverage = "results/03_host_mapping/{sample}_coverage.tsv",
        hist = "results/03_host_mapping/{sample}_coverage_histogram.txt",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_map_reads_to_host",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:00:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/03_host_mapping/{sample}_map_reads_to_host.log"
    benchmark: "results/03_host_mapping/{sample}_map_reads_to_host.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa mem -M -t {threads} {input.host_db} {input.reads1} {input.reads2} \
        | samtools view -bh | samtools sort - > {output.bam}
        samtools coverage {output.bam} > {output.coverage}
        samtools flagstat -O tsv {output.bam} > {output.flagstat}
        samtools coverage {output.bam} -m > {output.hist}
        samtools view -bh -f4 {output.bam} | samtools sort - > {output.bam_hostfiltered}
        """

rule get_non_host_reads:
    input:
        bam_hostfiltered = "results/03_host_mapping/{sample}_hostfiltered.bam"
    output:
        reads1 = "results/05_assembly/hostfiltered_reads/{sample}_R1.fastq.gz",
        reads2 = "results/05_assembly/hostfiltered_reads/{sample}_R2.fastq.gz",
    params:
        java_mem="16",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_get_non_host_reads",
        account="pengel_spirit",
        runtime_s=convertToSec("0-12:00:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/05_assembly/hostfiltered_reads/{sample}_get_non_host_reads.log"
    benchmark: "results/05_assembly/hostfiltered_reads/{sample}_get_non_host_reads.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        samtools index {input.bam_hostfiltered}
        # used the below command for A2-2 alone
        # picard -Xmx{params.java_mem}g SamToFastq I={input.bam_hostfiltered} F={output.reads1} F2={output.reads2} VALIDATION_STRINGENCY=SILENT &> {log}
        outreads1={output.reads1}
        outreads2={output.reads2}
        bedtools bamtofastq -i {input.bam_hostfiltered} -fq ${{outreads1/.gz/}} -fq2 ${{outreads2/.gz/}} &> {log}
        gzip ${{outreads1/.gz/}}
        gzip ${{outreads2/.gz/}}
        """
        
def reads_for_assembly(sample, R_n):
    if sample in LARGE_SAMPLES:
        if sample in SAMPLES_TO_HOST_FILTER:
            return f"results/05_assembly/hostfiltered_reads/{sample}_R{R_n}.fastq.gz"
        else:
            return f"results/05_assembly/bbnormed_reads/{sample}_R{R_n}.fastq.gz"
    else:
        return f"results/01_trimmedconcatreads/{sample}_R{R_n}.fastq.gz"

rule assemble_metagenomes:
    input:
        reads1 = lambda wildcards: reads_for_assembly(wildcards.sample, 1),
        reads2 = lambda wildcards: reads_for_assembly(wildcards.sample, 2),
    output:
        scaffolds_unparsed = temp("results/05_assembly/all_reads_assemblies/{sample}_scaffolds_unparsed.fasta"),
        scaffolds = "results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta",
        graph = "results/05_assembly/all_reads_assemblies/{sample}_assembly_graph.fastg",
        spades_log = "results/05_assembly/all_reads_assemblies/{sample}_spades.log",
    params:
        outdir = lambda wildcards: "results/05_assembly/all_reads_assemblies/"+wildcards.sample,
        length_t = 1000,
        cov_t = 1,
        memory_limit = lambda wildcards, resources: "500",
        # memory_limit = lambda wildcards, resources: "850" if resources.attempt > 1 else "250",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards, resources: convertToSec("3-00:00:00"),
        # runtime_s=lambda wildcards, resources: convertToSec("1-20:00:00") if resources.attempt > 1 else convertToSec("0-20:00:00"),
    resources:
        attempt = lambda wildcards, attempt: attempt,
        mem_mb = lambda wildcards, attempt: convertToMb("500G")
        # mem_mb = lambda wildcards, attempt: convertToMb("800G") if attempt > 1 else convertToMb("200G")
    retries: 0
    # retries: 2
    threads: 5
    # threads: 4
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

# headers and IDs in the gff file are renamed to include sample name
# corresponding ffn file will have the header ID can be obtailed as 
# example: header=">C4-3_NODE_1_length_1160554_cov_128.405656_131"
# def get_ID_from_header(header):    
#     header_text = header.split(">")[1]
#     ID_1 = header_text.split("_NODE_")[0]
#     ID_2 = header_text.split("_NODE_")[1].split("_")[0]
#     ID_3 = header_text.split("_NODE_")[1].split("_")[-1]
#     ID = "_".join([ID_1,ID_2,ID_3])
#     return ID

rule rename_gff_headers:
    input:
        gff = "results/06_metagenomicORFs/{sample_assembly}/{sample_assembly}.gff"
    output:
        gff = "results/06_metagenomicORFs/{sample_assembly}.gff"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="rename_gff",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:30:00"),
    log: "results/06_metagenomicORFs/{sample_assembly}_rename_gff_headers.log"
    benchmark: "results/06_metagenomicORFs/{sample_assembly}_rename_gff_headers.benchmark"
    threads: 4
    resources:
        mem_mb = 3000,
    shell:
        """
        cat {input.gff} | sed -e 's/ID=/ID={wildcards.sample_assembly}_/g' | sed -e 's/NODE/{wildcards.sample_assembly}_NODE/g' > {output.gff}
        """

rule bamQC:
    input:
        bam = "results/07_MAG_binng_QC/01_backmapping/{sample_assembly}/{sample_assembly}.bam", # only unmapped reads excluded
        gff = "results/06_metagenomicORFs/{sample_assembly}.gff"
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
        mem_mb = convertToMb("300G"),
    conda: "../config/envs/qualimap-env.yaml"
    priority: 10
    shell:
        """

        qualimap bamqc -bam {input} -outdir {output.outdir} -outformat html \
        --java-mem-size={params.java_mem} --paint-chromosome-limits \
        --feature-file {input.gff} --outside-stats
        """

rule assembly_summary:
    input:
        flagstat = expand("results/07_MAG_binng_QC/01_backmapping/{sample}/{sample}_mapped_flagstat.txt", sample = SAMPLES_INDIA+SAMPLES_MY),
        scaffolds = expand("results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta" , sample = SAMPLES_INDIA+SAMPLES_MY),
    output:
        outfile = "results/05_assembly/assembly_summary.txt",
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

