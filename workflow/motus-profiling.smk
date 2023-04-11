#!/usr/bin/env python

"""
name: motu-profiling
description: run motus on trimmed reads for an initial picture of the community profile
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - run_motus
        + runs motus for each sample
    - merge_motus
        + merge all motus output into one file for downstream processing
targets:
    - motus_merged = "results/02_motus_profile/samples_merged.motus", # got to its respective rule and add desired list of samples in its expansion eg. SAMPLES
"""

rule run_motus:
    input:
        reads1 = "results/01_trimmedconcatreads/{sample}_R1.fastq.gz",
        reads2 = "results/01_trimmedconcatreads/{sample}_R2.fastq.gz",
    output:
        motus_temp = "results/02_motus_profile/{sample}.motus"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
        jobname="run_motus_{sample}"
    resources:
        mem_mb = convertToMb("50G")
    threads: 8
    log: "results/02_motus_profile/{sample}_run_motus.log"
    benchmark: "results/02_motus_profile/{sample}_run_motus.benchmark"
    conda: "config/envs/motus-env.yaml"
    shell:
        """
        motus downloadDB &> {log}
        motus profile -f {input.reads1} -r {input.reads2} -n {wildcards.sample} -o {output.motus_temp}  -t {threads} &>> {log}
        """

rule merge_motus:
    input:
        motus_temp = expand("results/02_motus_profile/{sample}.motus", sample=SAMPLES)
    output:
        motus_merged = "results/02_motus_profile/samples_merged.motus"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
        jobname="merge_motus"
    resources:
        mem_mb = convertToMb("50G")
    threads: 8
    log: "results/02_motus_profile/merge_motus.log"
    benchmark: "results/02_motus_profile/merge_motus.benchmark"
    conda: "config/envs/motus-env.yaml"
    shell:
        """
        motus merge -i $(echo \"{input.motus_temp}\" | sed -e 's/ /,/g' ) > {output.motus_merged} 2> {log}
        """
