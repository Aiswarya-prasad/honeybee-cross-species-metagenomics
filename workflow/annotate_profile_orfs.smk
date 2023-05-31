"""
name: annotate-orfd
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

rule cdhit_clustering:
    input:
        scaffolds_ffn = expand("results/06_metagenomicORFs/{sample}/{sample}.ffn", sample=SAMPLES_INDIA+SAMPLES_MY),
        scaffolds_faa = expand("results/06_metagenomicORFs/{sample}/{sample}.faa", sample=SAMPLES_INDIA+SAMPLES_MY),
    output:
        gene_catalog_ffn="results/08_gene_content/gene_catalog_all.ffn",
        gene_catalog_faa="results/08_gene_content/gene_catalog_all.faa",
        cdhit_genes="results/08_gene_content/gene_catalog_cdhit9590.fasta"
    params:
        mailto="a...
        runtime_s=convertToSec("0-70:00:00"),
    resources:
        mem_mb = convertToMb("512G")
    threads: 16
    log: "results/08_gene_content/cdhit_clustering.log"
    # benchmark: "results/08_gene_content/cdhit_clustering.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        cat {input.scaffolds_ffn} > {output.gene_catalog_ffn}
        cat {input.scaffolds.faa} > {output.gene_catalog_faa}
        echo "starting cd-hit-est"
        cd-hit-est -i {output.gene_catalog_ffn} -o {output.cdhit_genes} \
            -c 0.95 -T 64 -M 0 -G 0 -aS 0.9 -g 1 -r 1 -d 0
        """

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

rule index_gene_catalog:
    input:
        gene_catalog = "results/08_gene_content/gene_catalog_cdhit9590.fasta",
    output:
        bwa_index = multiext("results/08_gene_content/gene_catalog_cdhit9590.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="build_bwa_index_catalog",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/08_gene_content/bwa_index_catalog.log"
    benchmark: "results/08_gene_content/bwa_index_catalog.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa index {input.gene_catalog} &> {log}
        """

rule profile_genes:
    input:
        reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
        reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
        gene_catalog = "results/08_gene_content/gene_catalog_cdhit9590.fasta",
        bwa_index = multiext("results/08_gene_content/gene_catalog_cdhit9590.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam = "results/08_gene_content/01_profiling/{sample}_mapped.bam",
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
        runtime_s=convertToSec("0-5:10:00"),
    resources:
        mem_mb = convertToMb("20G")
    threads: 4
    log: "results/08_gene_content/01_profiling/{sample}_backmapping.log"
    benchmark: "results/08_gene_content/01_profiling/{sample}_backmapping.benchmark"
    conda: "../config/envs/mapping-env.yaml"
    shell:
        """
        bwa mem -a -t {threads} {input.scaffolds} {input.reads1} {input.reads2} \
        | samtools view -F 4 -h - |  python3 {params.filter_script} -e 5 -m 50 | samtools sort -O bam -@ {threads} > {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        """