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

rule rename_scaffolds:
    input:
        scaffolds = "results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta",
    output:
        scaffolds = "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample}_scaffolds.fasta",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="build_bwa_index",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample}_rename_scaffolds.log"
    benchmark: "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample}_rename_scaffolds.benchmark"
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python3 scripts/rename_scaffolds.py --scaffolds_in {input.scaffolds} --scaffolds_out {output.scaffolds} --sample {wildcards.sample} &>> {log}
        """

rule prodigal_get_orfs:
    input:
        scaffolds = "results/07_MAG_binng_QC/00_assembled_scaffolds/{sample}/{sample}_scaffolds.fasta"
    output:
        scaffolds_ffn = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.ffn",
        scaffolds_faa = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.faa",
        scaffolds_gff = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.gff"
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
        """
        # python scripts/filt_orfs.py --ffn_in {output.scaffolds_ffn} --ffn_out {output.orfs} --sample {wildcards.sample} --log {output.filt_log}

rule run_whokaryote:
    input:
        gff_input = "results/06_metagenomicORFs/{sample}/{sample}/prodigal_out.gff",
    output:
        whokaryote_eu = "results/05_assembly/contig_fates/whokaryote/{sample}/eukaryote_contig_headers.txt",
        whokaryote_pro = "results/05_assembly/contig_fates/whokaryote/{sample}/prokaryote_contig_headers.txt",
        whokaryote_out = "results/05_assembly/contig_fates/whokaryote/{sample}/whokaryote_predictions_S.tsv "
    params:
        outdir="results/05_assembly/contig_fates/whokaryote/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    log: "results/05_assembly/contig_fates/whokaryote/{sample}_whokaryote.log"
    benchmark: "results/05_assembly/contig_fates/whokaryote/{sample}_whokaryote.benchmark.txt"
    conda: "../config/envs/mags-env.yaml"
    shell:
        """
        whokaryote.py --outdir {params.outdir} --gff {input.gff_input} --model S &> {log}
        """

rule prodigal_filt_orfs:
    input:
        scaffolds_ffn = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.ffn",
        whokaryote_result = "results/05_assembly/contig_fates/whokaryote/{sample}/whokaryote_predictions_S.tsv",
        # replace with kraken later if better
        tax_result = "results/05_assembly/contig_fates/kaiju/nr/{sample}_fullnames.txt",
    output:
        scaffolds_ffn = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.ffn",
        filt_log = "results/06_metagenomicORFs/{sample}/orfs_filt_sumary.log",
    params:
        taxtype = "kaiju", # kaiju or kraken
        outdir = "results/06_metagenomicORFs/{sample}/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 16000
    threads: 8
    log: "results/06_metagenomicORFs/{sample}/{sample}_prodigal_filt_orfs.log"
    benchmark: "results/06_metagenomicORFs/{sample}/{sample}_prodigal_filt_orfs.benchmark"
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python scripts/filt_orfs.py --who {input.whokaryote_result} \
            --tax {input.tax_result}  --taxtype {params.taxtype} \
            --ffn_in {input.scaffolds_ffn} --ffn_out {output.scaffolds_ffn} \
            --sample {wildcards.sample} --log {output.filt_log}
        """

# # ffn to gff rule also
# rule get_filt_orfs_gff:
#     input:
#         filt_ffn = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.ffn",
#         gff = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.gff",
#     output:
#         gff_filt = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.gff",

rule get_filt_orfs_faa:
    input:
        filt_ffn = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.ffn",
        faa = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.faa",
    output:
        headers = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}_headers.txt",
        filt_faa = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.faa"
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
        scaffolds_ffn = expand("results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.ffn", sample=SAMPLES_INDIA+SAMPLES_MY),
        scaffolds_faa = expand("results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.faa", sample=SAMPLES_INDIA+SAMPLES_MY),
    output:
        gene_catalog_ffn="results/08_gene_content/20230313_gene_catalog.ffn",
        gene_catalog_faa="results/08_gene_content/20230313_gene_catalog.faa",
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
        dram_config = "config/dram_config.json",
        filt_faa = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.faa",
    output:
        dram_annotations = "results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv",
    params:
        db_location = "/reference/dram_20230610",
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
        ###
        rm -rf ${{dram_outdir}} &>> {log} # snakemake creates it but DRAM will complain
        DRAM-setup.py import_config --config_loc {input.dram_config} &>> {log}
        DRAM.py annotate_genes -i {input.filt_faa} -o ${{dram_outdir}} \
                --config_loc {input.dram_config} \
                --threads {threads} --verbose &>> {log}
        """
        # source /etc/profile.d/lmodstacks.sh &>> {log}
        # dcsrsoft use old &>> {log}
        # export PATH=/dcsrsoft/spack/external/dram/v1.2.4/bin:$PATH &>> {log}
        # module load gcc/9.3.0 python &>> {log}
        # module load hmmer mmseqs2 prodigal infernal trnascan-se barrnap &>> {log}
        # which DRAM.py &>> {log}
        # dram_annotations={output.dram_annotations} &>> {log}
        # dram_outdir=${{dram_annotations/annotations.tsv}} &>> {log}
        # DRAM-setup.py version &>> {log}

rule index_nr_gene_catalog:
    input:
        gene_catalog = "results/08_gene_content/20230313_gene_catalog.ffn",
    output:
        bwa_index = multiext("results/08_gene_content/20230313_gene_catalog.ffn", ".amb", ".ann", ".bwt", ".pac", ".sa"),
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
        gene_catalog = "results/08_gene_content/20230313_gene_catalog.ffn",
        bwa_index = multiext("results/08_gene_content/20230313_gene_catalog.ffn", ".amb", ".ann", ".bwt", ".pac", ".sa"),
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
        runtime_s=convertToSec("0-10:10:00"),
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
