"""
name: annotate-orfs
description: ORF annotation and extraction of coverage
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - rename_scaffolds
        + renaming the scaffolds to have a unique name
    - prodigal_get_orfs
        + running prodigal to get ORFs
    - run_whokaryote
        + running whokaryote to predict eukaryotic and prokaryotic contigs
    - prodigal_filt_orfs
        + filtering ORFs based on whokaryote predictions
    - get_filt_orfs_gff
        + getting the gff file for the filtered ORFs
    - get_filt_orfs_faa
        + getting the faa file for the filtered ORFs
    - cdhit_clustering
        + clustering the ORFs (not used for further analysis per se)
    - parse_clustering_file
        + parsing the clustering file to get the cluster-host affiliations
    - dram_annotate_orfs
        + annotating the ORFs using DRAM KO and Pfam information referred to later
    - index_nr_gene_catalog
        + indexing the gene catalog for mapping to the reads (not used for further analysis per se)
    - profile_genes
        + mapping the reads to the gene catalog and getting coverage information
    - bowtie2_build
        + building the bowtie2 index for the gene catalog
    - map_to_assembly
        + mapping the reads to the assembly using bowtie2 (insert size filtering was not working with bwa in instrain)
    - get_sorted_genes
        + getting the sorted genes from the gff file
    - filter_bam_by_qual
        + filtering the bam file by quality by removing reads with quality less than 20
    - gene_counts_coverage
        + getting the gene counts and coverage information
    - run_kaiju_genes
        + running kaiju on the gene catalog
    - run_kaiju_genes_taxonomy
        + running kaiju on the gene catalog with full taxonomy
    - run_kraken2_genes
        + running kraken2 on the gene catalog
    - cayman_profiling
        + running cayman on the gene catalog
scripts:
    scripts/rename_scaffolds.py
    scripts/filt_orfs.py
    scripts/get_filt_gff.py
    scripts/clean_faa_header.py
    scripts/parse_cluster_file.py
    scripts/filter_bam.py
    scripts/filter_bam.py
    scripts/gff_to_bed.py
"""

rule rename_scaffolds
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

rule prodigal_get_orfs
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

rule run_whokaryote:
    input:
        gff_input = "results/06_metagenomicORFs/{sample}/{sample}/prodigal_out/{sample}.gff",
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

rule get_filt_orfs_gff:
    input:
        filt_ffn = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.ffn",
        gff = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.gff",
    output:
        filt_gff = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.gff"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    log: "results/06_metagenomicORFs/{sample}_get_filt_orfs_gff.log"
    benchmark: "results/06_metagenomicORFs/{sample}_get_filt_orfs_gff.benchmark"
    conda: "../config/envs/genes-env.yaml"
    threads: 2
    shell:
        """
        python3 scripts/get_filt_gff.py --ffn_in {input.filt_ffn} --gff_in {input.gff} --gff_out {output.filt_gff}
        """

rule get_filt_orfs_faa:
    input:
        filt_ffn = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.ffn",
        faa = "results/06_metagenomicORFs/{sample}/prodigal_out/{sample}.faa",
    output:
        headers = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}_headers.txt",
        unfilt_faa = temp("results/06_metagenomicORFs/{sample}/filt_orfs/{sample}_unfilt.faa"),
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
        python3 scripts/clean_faa_header.py -i {input.faa} > {output.unfilt_faa}
        seqtk subseq {output.unfilt_faa} {output.headers} > {output.filt_faa}
        """

rule cdhit_clustering:
    input:
        scaffolds_ffn = expand("results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.ffn", sample=SAMPLES_INDIA+SAMPLES_MY),
        scaffolds_faa = expand("results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.faa", sample=SAMPLES_INDIA+SAMPLES_MY),
    output:
        gene_catalog_ffn="results/08_gene_content/20230313_gene_catalog.ffn",
        gene_catalog_faa="results/08_gene_content/20230313_gene_catalog.faa",
        cdhit_clustering="results/08_gene_content/00_cdhit_clustering/20230313_gene_catalog_cdhit9590.fasta.clstr",
        cdhit_genes="results/08_gene_content/00_cdhit_clustering/20230313_gene_catalog_cdhit9590.fasta"
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
        cdhit_clustering="results/08_gene_content/00_cdhit_clustering/20230313_gene_catalog_cdhit9590.fasta.clstr"
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
    conda: "../config/envs/dram-env.yaml"
    shell:
        """
        which DRAM.py &>> {log}
        dram_annotations={output.dram_annotations} &>> {log}
        dram_outdir=${{dram_annotations/annotations.tsv}} &>> {log}
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
# We do not use this anymore because it is better to get covarage
# from mapping to the assembly
# rule profile_genes:
#     input:
#         reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
#         reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
#         gene_catalog = "results/08_gene_content/20230313_gene_catalog.ffn",
#         bwa_index = multiext("results/08_gene_content/20230313_gene_catalog.ffn", ".amb", ".ann", ".bwt", ".pac", ".sa"),
#     output:
#         bam = temp("results/08_gene_content/01_profiling/{sample}_mapped.bam"),
#         flagstat = "results/08_gene_content/01_profiling/{sample}_mapped.flagstat",
#         coverage = "results/08_gene_content/01_profiling/{sample}_mapped.coverage",
#         hist = "results/08_gene_content/01_profiling/{sample}_mapped.hist",
#     params:
#         match_length = 50,
#         edit_distance = 5, # methods in microbiomics recommends 95 perc identity
#         # since reads are 150 bp long, 5 mismatches is 3.3% mismatch which is almost as instrain recommends
#         # even less chances of strains mismapping
#         filter_script = "scripts/filter_bam.py",
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         jobname="{sample}_profile_genes",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-25:10:00"),
#     resources:
#         mem_mb = convertToMb("40G")
#     threads: 4
#     log: "results/08_gene_content/01_profiling/{sample}_profile_genes.log"
#     benchmark: "results/08_gene_content/01_profiling/{sample}_profile_genes.benchmark"
#     conda: "../config/envs/mapping-env.yaml"
#     shell:
#         """
#         bwa mem -a -t {threads} {input.gene_catalog} {input.reads1} {input.reads2} \
#         | samtools view -F 4 -h - |  python3 {params.filter_script} -e 5 -m 50 | samtools sort -O bam -@ {threads} > {output.bam}
#         samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
#         samtools coverage {output.bam} > {output.coverage}
#         samtools coverage -m {output.bam} > {output.hist}
#         """
#         # awk '$6 > 50' {input_f} > {output_f}

rule bowtie2_build:
    input:
        assembly = "results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta",
    output:
        bowtie_index = multiext("results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/08_gene_content/01_profiling_bowtie2/{sample}_bowtie2_build.log"
    benchmark: "results/08_gene_content/01_profiling_bowtie2/{sample}_bowtie2_build.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2-build {input.assembly} {input.assembly}
        """

rule map_to_assembly:
    input:
        reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
        reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
        assembly = "results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta",
        bowtie_index = multiext("results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    output:
        bam = temp("results/08_gene_content/01_profiling_bowtie2/{sample}_mapped.bam"),
        flagstat = "results/08_gene_content/01_profiling_bowtie2/{sample}_mapped.flagstat",
        hist = "results/08_gene_content/01_profiling_bowtie2/{sample}_mapped.hist",
    params:
        match_length = 50,
        edit_distance = 5, # methods in microbiomics recommends 95 perc identity
        # since reads are 150 bp long, 5 mismatches is 3.3% mismatch which is almost as instrain recommends
        # even less chances of strains mismapping
        filter_script = "scripts/filter_bam.py",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-15:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/08_gene_content/01_profiling_bowtie2/{sample}_map_to_assembly.log"
    benchmark: "results/08_gene_content/01_profiling_bowtie2/{sample}_map_to_assembly.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2 -x {input.assembly} -1 {input.reads1} -2 {input.reads2} -p {threads} \
        | samtools view -F 4 -h - |  python3 {params.filter_script} -e 5 -m 50 | samtools sort -O bam -@ {threads} > {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        samtools coverage {output.bam} > {output.coverage}
        samtools coverage -m {output.bam} > {output.hist}
        """

rule get_sorted_genes:
    input:
        gff = "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.gff",
    output:
        bed = "results/08_gene_content/01_profiling_bowtie2/{sample}_filt_genes.bed",
        gff = "results/08_gene_content/01_profiling_bowtie2/{sample}_filt_genes.gff"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:30:00"),
    resources:
        mem_mb = convertToMb("50G")
    log: "results/08_gene_content/01_profiling_bowtie2/{sample}_get_sorted_genes.log"
    benchmark: "results/08_gene_content/01_profiling_bowtie2/{sample}_get_sorted_genes.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        # remove the sample name from the headers
        bedtools sort -i {input.gff} | sed -e 's/{wildcards.sample}_//' > {output.gff}
        python3 scripts/gff_to_bed.py --gff {output.gff} --bed {output.bed}
        """

rule filter_bam_by_qual:
    input:
        bam = "results/08_gene_content/01_profiling_bowtie2/{sample}_mapped.bam",
    output:
        bam = "results/08_gene_content/01_profiling_bowtie2/{sample}_mapped_Q20.bam",
        flagstat = "results/08_gene_content/01_profiling_bowtie2/{sample}_mapped_Q20.flagstat",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:30:00"),
    resources:
        mem_mb = convertToMb("50G")
    log: "results/08_gene_content/01_profiling_bowtie2/{sample}_filter_bam_by_qual.log"
    benchmark: "results/08_gene_content/01_profiling_bowtie2/{sample}_filter_bam_by_qual.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        samtools view -q 20 -b {input.bam} > {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        """

rule gene_counts_coverage:
    input:
        bam = "results/08_gene_content/01_profiling_bowtie2/{sample}_mapped_Q20.bam",
        bed = "results/08_gene_content/01_profiling_bowtie2/{sample}_filt_genes.bed"
    output:
        coverage = "results/08_gene_content/01_profiling_bowtie2/{sample}_gene_coverage.txt",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-01:00:00"),
    resources:
        mem_mb = lambda wildcards: convertToMb("350G") if wildcards.sample != "A4-4" else convertToMb("550G") # this bam was 11GB big
    log: "results/08_gene_content/01_profiling_bowtie2/{sample}_gene_counts_coverage.log"
    benchmark: "results/08_gene_content/01_profiling_bowtie2/{sample}_gene_counts_coverage.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.bam} > {output.coverage}
        """

rule run_kaiju_genes:
    input:
        ffn_input = "results/08_gene_content/20230313_gene_catalog.faa",
    output:
        kaiju_out = "results/08_gene_content/04_kaiju_on_genes/nr/20230313_gene_catalog.kaiju",
    params:
        kaiju_db_nodes = "/...<kaiju_dir_path>.../kaiju_db/nr/nodes.dmp",
        kaiju_db_fmi = "/...<kaiju_dir_path>.../kaiju_db/nr/nr/kaiju_db_nr.fmi",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    resources:
        mem_mb = convertToMb("300G")
    threads: 4
    log: "results/08_gene_content/04_kaiju_on_genes/nr/kaiju.log"
    benchmark: "results/08_gene_content/04_kaiju_on_genes/nr/kaiju.benchmark"
    conda: "../config/envs/kaiju_env.yaml"
    shell:
        """
        kaiju -X -t {params.kaiju_db_nodes} -f {params.kaiju_db_fmi} -p \
                -o {output.kaiju_out} -z {threads} \
                -i {input.ffn_input} -v -a mem &> {log}
        """

rule run_kaiju_genes_taxonomy:
    input:
        kaiju_out = "results/08_gene_content/04_kaiju_on_genes/nr/20230313_gene_catalog.kaiju",
    output:
        kaiju_names = "results/08_gene_content/04_kaiju_on_genes/nr/20230313_gene_catalog_taxa.txt",
        kaiju_names_full = "results/08_gene_content/04_kaiju_on_genes/nr/20230313_gene_catalog_taxa_full.txt",
    params:
        kaiju_db_nodes = "/...<kaiju_dir_path>.../kaiju_db/nr/nodes.dmp",
        kaiju_db_names = "/...<kaiju_dir_path>.../kaiju_db/nr/names.dmp",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-1:10:00"),
    resources:
        mem_mb = convertToMb("100G")
    threads: 1
    log: "results/08_gene_content/04_kaiju_on_genes/nr/kaiju_names.log"
    benchmark: "results/08_gene_content/04_kaiju_on_genes/nr/kaiju_names.benchmark"
    conda: "../config/envs/kaiju_env.yaml"
    shell:
        """
        kaiju-addTaxonNames -t {params.kaiju_db_nodes} \
            -n {params.kaiju_db_names} -i {input.kaiju_out} \
            -o {output.kaiju_names_full} -v &> {log}
        kaiju-addTaxonNames -p -t {params.kaiju_db_nodes} \
            -n {params.kaiju_db_names} -i {input.kaiju_out} \
            -o {output.kaiju_names} -v &> {log}
        """

rule run_kraken2_genes:
    input:
        ffn_input = "results/08_gene_content/20230313_gene_catalog.ffn",
    output:
        kraken_report = "results/08_gene_content/04_kraken2_on_genes/20230313_gene_catalog_report.txt",
    params:
        kraken_out = "results/08_gene_content/04_kraken2_on_genes/20230313_gene_catalog.kraken",
        kraken_db = "data/220131_costum_kraken2db",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:10:00"),
    resources:
        mem_mb = convertToMb("300G")
    threads: 8
    log: "results/08_gene_content/04_kraken2_on_genes/kraken2.log"
    benchmark: "results/08_gene_content/04_kraken2_on_genes/kraken2.benchmark"
    conda: "../config/envs/kraken_env.yaml"
    shell:
        """
        kraken2 --use-names --threads {threads} --db {params.kraken_db} \
                --report {output.kraken_report} --output {params.kraken_out} \
                 {input.ffn_input} &> {log}
        """

rule cayman_profiling:
    input:
        reads1 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R1_repaired.fastq.gz",
        reads2 = lambda wildcards: f"results/01_cleanreads/{wildcards.sample}_R2_repaired.fastq.gz",
        annot_db = "data/cayman_gene_db/20230313_gene_catalog_db.csv",
        gene_catalog = "results/08_gene_content/20230313_gene_catalog.ffn",
        bwa_index = multiext("results/08_gene_content/20230313_gene_catalog.ffn", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        aln_stats = "results/08_gene_content/06_cayman/{sample}.aln_stats.txt.gz",
        cazy_profile = "results/08_gene_content/06_cayman/{sample}.cazy.txt.gz",
        gene_counts = "results/08_gene_content/06_cayman/{sample}.gene_counts.txt.gz"
    params:
        prefix = "results/08_gene_content/06_cayman/{sample}",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:10:00"),
    resources:
        mem_mb = convertToMb("100G")
    threads: 8
    log: "results/08_gene_content/06_cayman/{sample}_cayman.log"
    benchmark: "results/08_gene_content/06_cayman/{sample}_cayman.benchmark"
    conda: "../config/envs/cayman-env.yaml"
    shell:
        """
        # better yet check if it is in path and if not, set it up
        # for now this fix is ok
        cayman -t {threads} --db_separator , --db_coordinates hmmer \
                                     -1 {input.reads1} \
                                     -2 {input.reads2} \
                                     --out_prefix {params.prefix} \
                                     {input.annot_db} \
                                     {input.gene_catalog} &> {log}
        """