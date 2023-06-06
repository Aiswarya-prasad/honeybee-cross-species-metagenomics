"""
name: annotate-orfs
description: ORF annotation and extraction of coverage
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - prodigal_get_orfs
        + identify ORFs in all filtered contigs and remove partial and short ORFs
    - dram_annotate_orfs
        + annotate ORFs with dram
scripts:
    - 
targets:
    - 
"""


rule rename_gff_headers:
    input:
        gff = "results/06_metagenomicORFs/{sample_assembly}/{sample_assembly}.gff"
    output:
        gff = "results/06_metagenomicORFs/{sample_assembly}/{sample_assembly}_renamed.gff"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="rename_gff",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:30:00"),
    log: "results/06_metagenomicORFs/{sample_assembly}/{sample_assembly}_rename_gff_headers.log"
    benchmark: "results/06_metagenomicORFs/{sample_assembly}/{sample_assembly}_rename_gff_headers.benchmark"
    threads: 4
    resources:
        mem_mb = 3000,
    shell:
        """
        cat {input.gff} | sed -e 's/ID=/ID={wildcards.sample_assembly}_/g' | sed -e 's/NODE/{wildcards.sample_assembly}_NODE/g' > {output.gff}
        """

rule rename_faa_headers:
    input:
        faa = "results/06_metagenomicORFs/{sample}/{sample}.faa"
    output:
        faa = "results/06_metagenomicORFs/{sample}/{sample}_renamed.faa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="rename_faa",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:30:00"),
    log: "results/06_metagenomicORFs/{sample}_rename_faa_headers.log"
    benchmark: "results/06_metagenomicORFs/{sample}_rename_faa_headers.benchmark"
    threads: 4
    resources:
        mem_mb = 3000,
    shell:
        """
        cat {input.faa} | sed -e 's/ID=/ID={wildcards.sample}_/g' | sed -e 's/NODE/{wildcards.sample}_NODE/g' > {output.faa}
        """

rule rename_ffn_headers:
    input:
        ffn = "results/06_metagenomicORFs/{sample}/{sample}.ffn"
    output:
        ffn = "results/06_metagenomicORFs/{sample}/{sample}_renamed.ffn"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="rename_ffn",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:30:00"),
    log: "results/06_metagenomicORFs/{sample}_rename_ffn_headers.log"
    benchmark: "results/06_metagenomicORFs/{sample}_rename_ffn_headers.benchmark"
    threads: 4
    resources:
        mem_mb = 3000,
    shell:
        """
        cat {input.ffn} | sed -e 's/ID=/ID={wildcards.sample}_/g' | sed -e 's/NODE/{wildcards.sample}_NODE/g' > {output.ffn}
        """

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
        # takes care of renaming the headers
        python scripts/filt_orfs.py --ffn_in {output.scaffolds_ffn} --ffn_out {output.orfs} --sample {wildcards.sample} --log {output.filt_log}
        """

rule get_filt_orfs_faa:
    input:
        filt_ffn = "results/06_metagenomicORFs/{sample}_orfs.ffn",
        faa = "results/06_metagenomicORFs/{sample}/{sample}_renamed.faa",
    output:
        headers = "results/06_metagenomicORFs/{sample}_orfs_headers.txt",
        filt_faa = "results/06_metagenomicORFs/{sample}_orfs.faa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    log: "results/06_metagenomicORFs/{sample}_get_filt_orfs_faa.log"
    benchmark: "results/06_metagenomicORFs/{sample}_get_filt_orfs_faa.benchmark"
    conda: "../config/envs/genes-env.yaml"
    threads: 2
    shell:
        """
        grep \"^>\" {input.filt_ffn} | cut -f 2 -d \">\" | cut -f 1 -d \" \" > {output.headers}
        seqtk subseq {input.faa} {output.headers} > {output.filt_faa}
        """

rule cdhit_clustering:
    input:
        scaffolds_ffn = expand("results/06_metagenomicORFs/{sample}_orfs.ffn", sample=SAMPLES_INDIA+SAMPLES_MY),
        scaffolds_faa = expand("results/06_metagenomicORFs/{sample}_orfs.faa", sample=SAMPLES_INDIA+SAMPLES_MY),
    output:
        gene_catalog_ffn="results/08_gene_content/gene_catalog_all.ffn",
        gene_catalog_faa="results/08_gene_content/gene_catalog_all.faa",
        cdhit_clustering="results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta.clstr",
        cdhit_genes="results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-70:00:00"),
    resources:
        mem_mb = convertToMb("512G")
    threads: 16
    log: "results/08_gene_content/00_cdhit_clustering/cdhit_clustering.log"
    benchmark: "results/08_gene_content/00_cdhit_clustering/cdhit_clustering.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        cat {input.scaffolds_ffn} > {output.gene_catalog_ffn}
        cat {input.scaffolds_faa} > {output.gene_catalog_faa}
        echo "starting cd-hit-est"
        cd-hit-est -i {output.gene_catalog_ffn} -o {output.cdhit_genes} \
            -c 0.95 -T 64 -M 0 -G 0 -aS 0.9 -g 1 -r 1 -d 0
        """
rule parse_clustering_file:
    input:
        cdhit_clustering="results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta.clstr"
    output:
        cdhit_clusters="results/08_gene_content/00_cdhit_clustering/cluster_host_affiliations.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:00:00")
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python3 scripts/parse_cluster_file.py --cluster_file {input.cdhit_clustering} --cluster_out {output.cdhit_clusters}
        """


rule dram_annotate_orfs:
    input:
        # dram_config = "config/dram_config.json"
        filt_faa = "results/06_metagenomicORFs/{sample}_orfs.faa",
    output:
        dram_annotations = "results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv",
    params:
        # db_location = "/reference/dram",
        dram_outdir = lambda wildcards: os.path.join("results/08_gene_content/02_DRAM_annotations/", wildcards.sample),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:00:00")
    resources:
        mem_mb = convertToMb("100G")
    threads: 8
    log: "results/08_gene_content/02_DRAM_annotations/{sample}_dram_annotate_orfs_{sample}.log"
    benchmark: "results/08_gene_content/02_DRAM_annotations/{sample}_dram_annotate_orfs_{sample}.benchmark"
    conda: "../config/envs/mags-env.yaml"
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
        DRAM.py annotate_genes -i {input.filt_faa} -o ${{dram_outdir}} --threads {threads} --verbose &>> {log}
        """

# rule dram_distill_orfs:
#     input:
#         dram_annotations = "results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv",

rule index_gene_catalog:
    input:
        gene_catalog = "results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta",
    output:
        bwa_index = multiext("results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="build_bwa_index_catalog",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/08_gene_content/00_cdhit_clustering/bwa_index_catalog.log"
    benchmark: "results/08_gene_content/00_cdhit_clustering/bwa_index_catalog.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa index {input.gene_catalog} &> {log}
        """

rule profile_genes:
    input:
        reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
        reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
        gene_catalog = "results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta",
        bwa_index = multiext("results/08_gene_content/00_cdhit_clustering/gene_catalog_cdhit9590.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam = temp("results/08_gene_content/01_profiling/{sample}_mapped.bam"),
        flagstat = "results/08_gene_content/01_profiling/{sample}_mapped.flagstat",
        depth = "results/08_gene_content/01_profiling/{sample}_mapped.depth",
        coverage = "results/08_gene_content/01_profiling/{sample}_mapped.coverage",
        hist = "results/08_gene_content/01_profiling/{sample}_mapped.hist",
    params:
        match_length = 50,
        edit_distance = 5, # methods in microbiomics recommends 95 perc identity
        # since reads are 150 bp long, 5 mismatches is 3.3% mismatch which is almost as instrain recommends
        # even less chances of strains mismapping
        filter_script = "scripts/filter_bam.py",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_profile_genes",
        account="pengel_spirit",
        runtime_s=convertToSec("0-15:10:00"),
    resources:
        mem_mb = convertToMb("40G")
    threads: 4
    log: "results/08_gene_content/01_profiling/{sample}_profile_genes.log"
    benchmark: "results/08_gene_content/01_profiling/{sample}_profile_genes.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa mem -a -t {threads} {input.gene_catalog} {input.reads1} {input.reads2} \
        | samtools view -F 4 -h - |  python3 {params.filter_script} -e 5 -m 50 | samtools sort -O bam -@ {threads} > {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        samtools depth -s -a {output.bam} > {output.depth}
        samtools coverage {output.bam} > {output.coverage}
        samtools coverage -m {output.bam} > {output.hist}
        """

# rule make_gene_counts_matrix: