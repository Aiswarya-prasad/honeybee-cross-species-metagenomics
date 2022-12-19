#!/usr/bin/env python

"""
name: assemble-qc
description: assembly of metagenomes, mapping and other steps for qc followed by ORF annotation and extraction of coverage
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
        + summarizes assmebly length, number of contigs before and after filtering and number of reads mapped for each sample
targets:
    - summary_assembly = "05_Assembly/scaffolds_mapping/Assembly_mapping_summary.csv"
    - dram_annotation_orfs = expand("06_MetagenomicORFs/dram_annotations/{sample}/annotations.tsv", sample=SAMPLES)
"""

rule assemble_metagenomes:
    input:
        reads1 = rules.trim.output.reads1,
        reads2 = rules.trim.output.reads2,
    output:
        scaffolds_unparsed = temp("05_Assembly/trimmed_reads/{sample}_scaffolds_unparsed.fasta"),
        scaffolds = "05_Assembly/trimmed_reads/{sample}_scaffolds.fasta",
        graph = "05_Assembly/trimmed_reads/{sample}_assembly_graph.fastg",
        spades_log = "05_Assembly/trimmed_reads/{sample}_spades.log",
    params:
        outdir = lambda wildcards: "05_Assembly/trimmed_reads/"+wildcards.sample,
        length_t = 1000,
        cov_t = 1,
        memory_limit = lambda wildcards, attempt: "650" if attempt > 1 else "250",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("1-20:00:00") if attempt > 1 else convertToSec("0-20:00:00"),
    resources:
        mem_mb = lambda wildcards: convertToMb("600G") if attempt > 1 else convertToMb("200G")
    retries: 2
    threads: 4
    log: "05_Assembly/trimmed_reads/{sample}_assemble_metagenomes.log"
    benchmark: "05_Assembly/trimmed_reads/{sample}_assemble_metagenomes.benchmark"
    conda: "envs/spades-env.yaml"
    shell:
        """
        if [ -f \"{params.outdir}/scaffolds.fasta\" ]; then
          echo {params.outdir}/scaffolds.fasta\" exists copying and cleaning\" | tee {log}
        else
          if [ -d \"{params.outdir}/\" ]; then
            echo {params.outdir}\" exists resuming spades\" | tee -a {log}
            spades.py --continue -o {params.outdir} || true &>> {log}
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

rule map_to_assembly:
    input:
        reads1 = rules.trim.output.reads1,
        reads2 = rules.trim.output.reads2,
        scaffolds = rules.assemble_metagenomes.output.scaffolds
    output:
        flagstat = "05_Assembly/scaffolds_mapping/{sample}_assembly_mapping_flagstat.tsv",
        sam = temp("05_Assembly/scaffolds_mapping/{sample}_assembly.sam"),
        bam = "05_Assembly/scaffolds_mapping/{sample}_assembly.bam",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_map_to_assembly",
        account="pengel_spirit",
        runtime_s=convertToSec("1-00:00:00"),
    resources:
        mem_mb = 8000
    threads: 5
    log: "05_Assembly/scaffolds_mapping/{sample}_map_to_assembly.log"
    benchmark: "05_Assembly/scaffolds_mapping/{sample}_map_to_assembly.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        mkdir -p 05_Assembly/scaffolds_mapping &> {log}
        bwa index {input.scaffolds} &>> {log}
        bwa mem -t {threads} {input.scaffolds} {input.reads1} {input.reads2} 1> {output.sam} 2>> {log}
        samtools view -bh {output.sam} | samtools sort - 1> {output.bam} 2>> {log}
        rm 05_Assembly/trimmed_reads/{wildcards.sample}_scaffolds.fasta.amb &>> {log}
        rm 05_Assembly/trimmed_reads/{wildcards.sample}_scaffolds.fasta.ann &>> {log}
        rm 05_Assembly/trimmed_reads/{wildcards.sample}_scaffolds.fasta.bwt &>> {log}
        rm 05_Assembly/trimmed_reads/{wildcards.sample}_scaffolds.fasta.pac &>> {log}
        rm 05_Assembly/trimmed_reads/{wildcards.sample}_scaffolds.fasta.sa &>> {log}
        samtools flagstat -O tsv {output.bam} > {output.flagstat} 2>> {log}
        """

rule summarize_mapping_assembly:
    input:
        scaffolds = expand("05_Assembly/trimmed_reads/{sample}_scaffolds.fasta", sample=SAMPLES),
        scaffolds_unparsed = expand("05_Assembly/trimmed_reads/{sample}_scaffolds_unparsed.fasta", sample=SAMPLES),
        flagstat = expand("05_Assembly/scaffolds_mapping/{sample}_assembly_mapping_flagstat.tsv", sample=SAMPLES)
    output:
        summary_assembly = "05_Assembly/scaffolds_mapping/Assembly_mapping_summary.csv",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="summarize_mapping_assembly",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 5
    log: "05_Assembly/scaffolds_mapping/summarize_mapping_assembly.log"
    benchmark: "05_Assembly/scaffolds_mapping/summarize_mapping_assembly.benchmark"
    shell:
        """
        python3 scripts/assembly_summary.py --scaffolds {input.scaffolds} --scaffolds_unparsed {input.scaffolds_unparsed} --flagstat {input.flagstat} --outfile {output.summary_assembly} &> {log}
        """

rule prodigal_get_orfs:
    input:
        scaffolds = "05_Assembly/trimmed_reads/{sample}_scaffolds.fasta"
    output:
        orfs = "06_MetagenomicORFs/{sample}_orfs.ffn",
        filt_log = "06_MetagenomicORFs/{sample}_orfs_filt_sumary.log",
        scaffolds_ffn = "06_MetagenomicORFs/{sample}/{sample}.ffn",
        scaffolds_faa = "06_MetagenomicORFs/{sample}/{sample}.faa",
        scaffolds_gff = "06_MetagenomicORFs/{sample}/{sample}.gff"
    params:
        outdir = "06_MetagenomicORFs/{sample}/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 16000
    threads: 8
    log: "06_MetagenomicORFs/{sample}_prodigal_get_orfs.log"
    benchmark: "06_MetagenomicORFs/{sample}_prodigal_get_orfs.benchmark"
    conda: "envs/snv-env.yaml"
    shell:
        """
        if [ ! -f {output.scaffolds_ffn} ]; then
            prodigal -i {input.scaffolds} -o {output.scaffolds_gff} -f gff -a {output.scaffolds_faa} -d {output.scaffolds_ffn} -p meta &> {log}
        fi
        python scripts/filt_orfs.py --ffn_in {output.scaffolds_ffn} --ffn_out {output.orfs} --sample {wildcards.sample} --log {output.filt_log}
        """

rule dram_annotate_orfs:
    input:
        scaffolds = "06_MetagenomicORFs/{sample}/{sample}.faa",
        dram_config = "config/dram_config.json"
    output:
        dram_annotations = "06_MetagenomicORFs/dram_annotations/{sample}/annotations.tsv",
    params:
        db_location = "/reference/dram",
        dram_outdir = lambda wildcards: os.path.join("06_MetagenomicORFs/dram_annotations/", wildcards.sample),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("1-20:00:00")
    resources:
        mem_mb = convertToMb("512G")
    threads: 2
    log: "06_MetagenomicORFs/dram_annotations/{sample}/dram_annotate_orfs_{sample}.log"
    benchmark: "06_MetagenomicORFs/dram_annotations/{sample}/dram_annotate_orfs_{sample}.benchmark"
    conda: "envs/mags-env.yaml"
    shell:
        """
        source /etc/profile.d/lmodstacks.sh &>> {log}
        dcsrsoft use old &>> {log}
        export PATH=/dcsrsoft/spack/external/dram/v1.2.4/bin:$PATH &>> {log}
        module load gcc/9.3.0 python &>> {log}
        module load hmmer mmseqs2 prodigal infernal trnascan-se barrnap &>> {log}
        which DRAM.py &>> {log}
        dram_annotations={output.dram_annotations} &>> {log}
        dram_outdir=${{dram_annotations/annotations.tsv}} &>> {log}
        DRAM-setup.py version &>> {log}
        ###
        rm -rf ${{dram_outdir}} &>> {log} # snakemake creates it but DRAM will complain
        DRAM.py annotate_genes -i {input.scaffolds} -o ${{dram_outdir}} --threads {threads} --verbose &>> {log}
        """

# to do next
# cluster genes...?
# rule to get coverage matrix of genes in all orfs