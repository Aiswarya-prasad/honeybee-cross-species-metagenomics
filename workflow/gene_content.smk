"""
name: gene_content.smk
description: xxx
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - xxx
scripts:
    - 
targets:
    - 
"""

##################################################################
# these lines are all here to make this smkm  self-contained for the moment
# will be changed later
import os
import glob
import yaml
import itertools
from itertools import chain

configfile: "config/config.yaml"
# read config information into local variables to improve readability of rules
if config["LocalBackup"]:
    localrules: backup

SAMPLES_KE = config["SAMPLES_KE"]
SAMPLES_INDIA = config["SAMPLES_INDIA"]
SAMPLES_MY = config["SAMPLES_MY"]
SAMPLES = SAMPLES_KE + SAMPLES_INDIA + SAMPLES_MY
# SAMPLES = config["SAMPLES_KE"]
ADAPTERS = config["Adapters"]

wildcard_constraints:
  sample = '|'.join(SAMPLES),
#   read = '|'.join(["R1", "R2"]),
#   lane = "L*",
  run = "20[0-9]{6,6}"

onstart:
    # this is just for updates - needs to be run before starting the pipeline
    shell("python3 scripts/make_reads_list_file.py")

raw_paths_dict_all = yaml.safe_load(open("config/raw_file_paths.yaml", "r"))
raw_paths_dict = {key: raw_paths_dict_all[key] for key in SAMPLES}

include: "common.smk"

rule all:
    input:
        cdhit_genes = "results/08_gene_content/gene_catalog_cdhit9590.fasta"
##################################################################


rule prodigal_get_orfs:
    input:
        scaffolds = "results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta"
    output:
        orfs = "results/06_metagenomicORFs/{sample}_orfs.ffn",
        filt_log = "results/06_metagenomicORFs/{sample}_orfs_filt_sumary.log",
        scaffolds_ffn = "results/06_metagenomicORFs/{sample}/{sample}.ffn",
        scaffolds_faa = "results/06_metagenomicORFs/{sample}/{sample}.faa",
        scaffolds_gff = "results/06_metagenomicORFs/{sample}/{sample}.gff"
    params:
        outdir = "results/06_metagenomicORFs/{sample}/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 16000
    threads: 8
    log: "results/06_metagenomicORFs/{sample}/{sample}_prodigal_get_orfs.log"
    benchmark: "results/06_metagenomicORFs/{sample}/{sample}_prodigal_get_orfs.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        if [ ! -f {output.scaffolds_ffn} ]; then
            prodigal -i {input.scaffolds} -o {output.scaffolds_gff} -f gff -a {output.scaffolds_faa} -d {output.scaffolds_ffn} -p meta &> {log}
        fi
        python scripts/filt_orfs.py --ffn_in {output.scaffolds_ffn} --ffn_out {output.orfs} --sample {wildcards.sample} --log {output.filt_log}
        """

rule cd_hit_clustering:
    input:
        scaffolds_ffn = expand("results/06_metagenomicORFs/{sample}/{sample}.ffn", sample = SAMPLES_INDIA+SAMPLES_MY),
        scaffolds_faa = expand("results/06_metagenomicORFs/{sample}/{sample}.faa", sample = SAMPLES_INDIA+SAMPLES_MY),
    output:
        gene_catalog_ffn = "results/08_gene_content/gene_catalog_all.ffn",
        gene_catalog_faa = "results/08_gene_content/gene_catalog_all.faa",
        cdhit_genes = "results/08_gene_content/gene_catalog_cdhit9590.fasta",
    params:
        identity_threshold = 0.95,
        memory_limit = 0, # no limit
        seq_id = 0, # if set to 0, then use local sequence identity
        aln_cov = 0.9,
        alg = 1, # 1: most similar cluster, 0: first matching cluster
        strands = 1, # do both +/+ & +/- alignments
        description = 0, # if set to 0, it takes the fasta defline and stops at first space
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-70:00:00"),
    resources:
        mem_mb = convertToMB("512G"),
    threads: 16
    log: "results/08_gene_content/cdhit_clustering.log"
    benchmark: "results/08_gene_content/cdhit_clustering.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        # used scripts/batch-scripts/cdhit_clustering.sh to run this rule earlier
        cat {input.scaffolds_ffn} > {output.gene_catalog_ffn}
        cat {input.scaffolds_faa} > {output.gene_catalog_faa}
        cd-hit-est -i {output.gene_catalog_ffn} -o {output.cdhit_genes} \
            -c {params.identity_threshold} -T {threads} -M {params.memory_limit} \
            -G {params.seq_id} -aS {params.aln_cov} -g {params.alg} \
            -r {params.strands} -d {params.description}
        """

# grep "^>" cdhit9590/gene_catalog_cdhit9590.fasta | \
# cut -f 2 -d ">" | \
# cut -f 1 -d " " > cdhit9590/cdhit9590.headers
# seqtk subseq gene_catalog_all.faa cdhit9590/cdhit9590.headers \
# > cdhit9590/gene_catalog_cdhit9590.faa

# rule dram_annotate_orfs:
#     input:
#         scaffolds = "results/06_metagenomicORFs/{sample}/{sample}.faa",
#         dram_config = "config/dram_config.json"
#     output:
#         dram_annotations = "results/06_metagenomicORFs/dram_annotations/{sample}/annotations.tsv",
#     params:
#         db_location = "/reference/dram",
#         dram_outdir = lambda wildcards: os.path.join("results/06_metagenomicORFs/dram_annotations/", wildcards.sample),
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         account="pengel_spirit",
#         runtime_s=convertToSec("1-20:00:00")
#     resources:
#         mem_mb = convertToMb("512G")
#     threads: 2
#     log: "results/06_metagenomicORFs/dram_annotations/{sample}/dram_annotate_orfs_{sample}.log"
#     benchmark: "results/06_metagenomicORFs/dram_annotations/{sample}/dram_annotate_orfs_{sample}.benchmark"
#     conda: "../config/envs/mags-env.yaml"
#     shell:
#         """
#         source /etc/profile.d/lmodstacks.sh &>> {log}
#         dcsrsoft use old &>> {log}
#         export PATH=/dcsrsoft/spack/external/dram/v1.2.4/bin:$PATH &>> {log}
#         module load gcc/9.3.0 python &>> {log}
#         module load hmmer mmseqs2 prodigal infernal trnascan-se barrnap &>> {log}
#         which DRAM.py &>> {log}
#         dram_annotations={output.dram_annotations} &>> {log}
#         dram_outdir=${{dram_annotations/annotations.tsv}} &>> {log}
#         DRAM-setup.py version &>> {log}
#         ###
#         rm -rf ${{dram_outdir}} &>> {log} # snakemake creates it but DRAM will complain
#         DRAM.py annotate_genes -i {input.scaffolds} -o ${{dram_outdir}} --threads {threads} --verbose &>> {log}
#         """

# to do next
# cluster genes...?
# rule to get coverage matrix of genes in all orfs