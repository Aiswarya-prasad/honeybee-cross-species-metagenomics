"""
This pipeline is written by Aiswarya Prasad (aiswarya.prasad@unil.ch) and is
intended primarily for my own use in my PhD thesis project(s). It is being written
to work on snakemake v6.15.5 and run in a cluster using the slurm profile mentioned
 here (https://github.com/RomainFeron/snakemake-slurm). It can be run with the appropriate lines being commented out
locally in a workstation also.
"""

import os
import itertools

configfile: "config/config.yaml"
LOCAL_BACKUP = config["LocalBackup"]
# get sample lists
SAMPLES_INDIA = config["SAMPLES_INDIA"]
SAMPLES_KE_CHOSEN = config["SAMPLES_KE_CHOSEN"]
# SAMPLES = SAMPLES_INDIA + SAMPLES_KE_CHOSEN
SAMPLES = SAMPLES_INDIA
# get important paths
PROJECT_IDENTIFIER = config["ProjectIdentifier"]
BACKUP_PATH = config["BackupPath"]
DBs = config["GENOME_DBs"]
ADAPTERS = config["Adapters"]
PROJECT_PATH = config["ProjectPath"]
GROUPS = ["g__Bombilactobacillus",
          "g__Lactobacillus",
          "g__Bifidobacterium",
          "g__Gilliamella",
          "g__Frischella",
          "g__Snodgrassella",
          "g__Bartonella",
          "g__Enterobacter",
          # "g__",
          "g__Pectinatus",
          "g__Apibacter",
          "g__Dysgonomonas",
          "g__Spiroplasma",
          # "g__Zymobacter",
          "g__Entomomonas",
          "g__Saezia",
          "g__Parolsenella",
          "g__WRHT01",
          "g__Commensalibacter",
          "g__Apilactobacillus",
          "g__Bombella"]

def get_g_dict_for_defined_groups(path):
    """
    Returns dictionary of genomes and groups with each value being a list of
    genomes corresponding to a given group
    """
    g_list_dict = {}
    for group in GROUPS:
        g_list_dict[group] = []
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            group = line.split("\t")[18]
            # only include groups of interest!
            if group not in GROUPS:
                continue
            g_list_dict[group].append(genome)
    return(g_list_dict)


def get_g_dict_for_groups_from_data(path):
    """
    Returns dictionary of genomes and groups with each value being a list of
    genomes corresponding to a given group
    """
    g_list_dict = {}
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            cluster = line.split("\t")[11]
            group = line.split("\t")[18]
            # only include groups of interest!
            if group == "g__":
                group = "g__"+cluster
            if group in g_list_dict.keys():
                g_list_dict[group].append(genome)
            else:
                g_list_dict[group] = [genome]
    return(g_list_dict)

def convertToMb(string):
    """
    This function can convert text in the form
    xxG to mb
    If it does not end with G, it returns the string
    It does not handle other cases of invalid input
    """
    if string.endswith("G"):
        number = int(string.split("G")[0])
        return(number*1000)
    else:
        return(string)

def convertToSec(string):
    """
    This function can convert text in the form
    D-hh:mm:ss to seconds
    D - # days
    hh # hours
    mm # mins
    ss # secs
    """
    days = string.split("-")[0]
    hrs = string.split("-")[1].split(":")[0]
    min = string.split("-")[1].split(":")[1]
    sec = string.split("-")[1].split(":")[2]
    total = int(sec)
    total = total + 60*int(min)
    total = total + 60*60*int(hrs)
    total = total + 24*60*60*int(days)
    return(total)

def make_count_dict(filepath):
    """
    This function is used by the rule "summarize mapping"
    to parse the file containing counts from bam files
    """
    count_dict = {}
    with open(filepath) as  fh:
        for line in fh:
            sample = line.split(" ")[0]
            count = int(line.split(" ")[1])
            if sample not in count_dict.keys():
                count_dict[sample] = count
    return(count_dict)

def num_genomes_in_group(group, path):
    """
    This function finds the number of genomes in group
    so make_tree is only run for those with at least 3
    """
    if group in GROUPS:
        return len(get_g_dict_for_defined_groups(path)[group])
    else:
        return len(get_g_dict_for_groups_from_data(path)[group])


def get_MAGs_list_dict(path, full_list=False):
    """
    read list of MAGs from checkpoint output
    (06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv) and return the list of Mags
    of a given sample either as a list per sample, as a complete list including
    all samples or a dictonary
    """
    g_paths_list = []
    g_list = []
    g_list_dict = {}
    for sample in SAMPLES:
        g_list_dict[sample] = []
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            sample = "_".join(genome.split("MAG_")[1].split("_")[:-1])
            # only include samples of interest!
            if sample not in SAMPLES:
                continue
            g_list_dict[sample].append(genome)
            g_list.append(genome)
        if full_list:
            return(g_list)
        return(g_list_dict)

def get_cluster_dict(path):
    """
    return a dict with keys as magOTU clusters with list of corresponding MAGs
    as values
    """
    cluster_list_dict = {}
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            type = line.split("\t")[10]
            if type != "MAGs":
                continue
            cluster = line.split("\t")[11]
            # only include groups of interest!
            if cluster not in cluster_list_dict.keys():
                cluster_list_dict[cluster] = [genome]
            else:
                cluster_list_dict[cluster].append(genome)
    return(cluster_list_dict)

if LOCAL_BACKUP:
    localrules: backup, concat_all_mags

rule all:
    input:
        html = PROJECT_IDENTIFIER+"_Report.html",
        backup_log = "logs/backup.log",


rule raw_qc:
    input:
        reads=ancient(os.path.join("00_RawData", "{sample}_{read}.fastq.gz")),
    output:
        html="fastqc/raw/{sample}_{read}_fastqc.html",
        zip="fastqc/raw/{sample}_{read}_fastqc.zip"
    params:
        outdir="fastqc/raw/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname=lambda wildcards: wildcards.sample+"_"+wildcards.read+"_qc",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{sample}_{read}_qc.log"
    benchmark: "logs/{sample}_{read}_qc.benchmark"
    threads: 2
    conda: "envs/trim-qc-env.yaml"
    shell:
        """
        fastqc -t {threads} {input.reads} -o {params.outdir}
        """

# """
# rule make_adapters:
#     input:
#         "config/index_table.csv"
#     output:
#         "config/Adapters-PE.fa"
#     script:
#         "scripts/write_adapters.py"
# """

rule trim:
    input:
        reads1=ancient(os.path.join("00_RawData", "{sample}_R1.fastq.gz")),
        reads2=ancient(os.path.join("00_RawData", "{sample}_R2.fastq.gz")),
        adapter="config/Adapters-PE.fa"
        # adapter=rules.make_adapters.output
    output:
        reads1 = "01_Trimmed/{sample}_R1_trim.fastq.gz",
        reads2 = "01_Trimmed/{sample}_R2_trim.fastq.gz",
        reads1_unpaired = temp("01_Trimmed/{sample}_R1.unpaired.fastq.gz"),
        reads2_unpaired = temp("01_Trimmed/{sample}_R2.unpaired.fastq.gz")
    params:
        # add adapter definition to config later
        adapters = lambda wildcards: ADAPTERS["SAMPLES_INDIA"] if wildcards.sample in SAMPLES_INDIA else (ADAPTERS["SAMPLES_KE_CHOSEN"] if wildcards.sample in SAMPLES_KE_CHOSEN else ADAPTERS["default"]),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_trim",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "logs/{sample}_trim.log"
    benchmark: "logs/{sample}_trim.benchmark"
    conda: "envs/trim-qc-env.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.reads1} {input.reads2} {output.reads1} {output.reads1_unpaired} {output.reads2} {output.reads2_unpaired} ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:28 TRAILING:28  MINLEN:60
        """

rule trim_qc:
    input:
        reads="01_Trimmed/{sample}_{read}_trim.fastq.gz",
    output:
        html="fastqc/trim/{sample}_{read}_trim_fastqc.html",
        zip="fastqc/trim/{sample}_{read}_trim_fastqc.zip"
    threads: 2
    conda: "envs/trim-qc-env.yaml"
    params:
        outdir="fastqc/trim/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_{read}_trim_qc",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "logs/{sample}_{read}_trim_qc.log"
    benchmark: "logs/{sample}_{read}_trim_qc.benchmark"
    resources:
        mem_mb = 8000
    shell:
        """
        fastqc -t {threads} {input.reads} -o {params.outdir}
        """

rule index_bwa:
    input:
        db = "database/{db}",
    output:
        out = multiext("database/{db}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    conda: "envs/mapping-env.yaml"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{db}_index_bwa",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{db}_index_bwa.log"
    benchmark: "logs/{db}_index_bwa.benchmark"
    shell:
        """
        # Because we want index files to be in this directory
        cd database/
        bwa index {wildcards.db}
        """

rule index_samtools:
    input:
        db = "database/{db}",
    output:
        out = "database/{db}.fai"
    conda: "envs/mapping-env.yaml"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{db}_index_samtools",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{db}_index_samtools.log"
    benchmark: "logs/{db}_index_samtools.benchmark"
    shell:
        """
        # Because we want index files to be in this directory
        cd database/
        samtools faidx {wildcards.db}
        """

rule make_genome_list:
    input:
        db = "database/{db}"
    output:
        txt = "database/{db}_list.txt"
    log: "logs/{db}_make_list.log"
    benchmark: "logs/{db}_make_list.benchmark"
    shell:
        "cat {input.db} | grep \">\" | sed \"s/>//\" > {output.txt}"

rule run_motus:
    input:
        reads1 = "01_Trimmed/{sample}_R1_trim.fastq.gz",
        reads2 = "01_Trimmed/{sample}_R2_trim.fastq.gz",
    output:
        motus_temp = "08_motus_profile/{sample}.motus"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
        jobname="run_motus_{sample}"
    resources:
        mem_mb = convertToMb("50G")
    threads: 8
    log: "logs/{sample}_run_motus.log"
    benchmark: "logs/{sample}_run_motus.benchmark"
    conda: "envs/motus-env.yaml"
    shell:
        """
        motus downloadDB # just leave it here to be safe - it warns and continues
        motus profile -f {input.reads1} -r {input.reads2} -n {wildcards.sample} -o {output.motus_temp}  -t {threads}
        """

rule merge_motus:
    input:
        motus_temp = expand("08_motus_profile/{sample}.motus", sample=SAMPLES)
    output:
        motus_merged = "08_motus_profile/samples_merged.motus"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
        jobname="merge_motus"
    resources:
        mem_mb = convertToMb("50G")
    threads: 8
    log: "logs/merge_motus.log"
    benchmark: "logs/merge_motus.benchmark"
    conda: "envs/motus-env.yaml"
    shell:
        """
        motus merge -i $(echo \"{input.motus_temp}\" | sed -e 's/ /,/g' ) > {output.motus_merged}
        """

rule host_mapping:
    input:
        reads1=rules.trim.output.reads1,
        reads2=rules.trim.output.reads2,
        index=multiext("database/4_host_db", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam=temp("02_HostMapping/{sample}_host.bam"),
        bam_unmapped=temp("02_HostMapping/{sample}_host_unmapped.bam"),
        sam=temp("02_HostMapping/{sample}_host.sam"),
        bam_mapped=temp("02_HostMapping/{sample}_host_mapped.bam"), # follow up on what the host mapped reads contain later and make this bam temporary
        flagstat_02 = "02_HostMapping/{sample}_flagstat.tsv",
        coverage_host = "02_HostMapping/{sample}_coverage.tsv",
        coverage_host_hist = "02_HostMapping/{sample}_coverage_histogram.txt",
    params:
        host_db = os.path.join(os.getcwd(), "database/4_host_db"), # from my database directory
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_host_mapping",
        account="pengel_spirit",
        runtime_s=convertToSec("0-4:10:00"),
    resources:
        mem_mb = convertToMb("20G")
    threads: 8
    log: "logs/{sample}_host_mapping.log"
    benchmark: "logs/{sample}_host_mapping.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        # map reads to host and extract the unmapped reads as fastq
        bwa mem -t {threads} {params.host_db} {input.reads1} {input.reads2} > {output.sam}
        samtools view -bh {output.sam} | samtools sort - > {output.bam}
        samtools view -bh -f4 {output.bam} | samtools sort - > {output.bam_unmapped}
        samtools view -bh -F4 {output.bam} | samtools sort - > {output.bam_mapped}
        samtools coverage {output.bam_mapped} > {output.coverage_host}
        samtools coverage {output.bam_mapped} -m > {output.coverage_host_hist}
        samtools flagstat -O tsv {output.bam} > {output.flagstat_02}
        """

rule host_mapping_extract_host_mapped_reads:
    input:
        bam_mapped=rules.host_mapping.output.bam_mapped
    output:
        reads1 = temp("02_HostMapping/{sample}_R1_host_mapped.fastq"),
        reads2 = temp("02_HostMapping/{sample}_R2_host_mapped.fastq"),
        bai_mapped=temp("02_HostMapping/{sample}_host_mapped.bam.bai")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_host_mapping_extract_host_mapped_reads",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "logs/{sample}_host_mapping_extract_host_mapped_reads.log"
    benchmark: "logs/{sample}_host_mapping_extract_host_mapped_reads.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        samtools index {input.bam_mapped}
        picard -Xmx8g SamToFastq I={input.bam_mapped} F={output.reads1} F2={output.reads2} VALIDATION_STRINGENCY=SILENT
        """

rule host_mapping_extract_host_filtered_reads:
    input:
        bam_unmapped=rules.host_mapping.output.bam_unmapped
    output:
       reads1 = "02_HostMapping/{sample}_R1_host_unmapped.fastq",
       reads2 = "02_HostMapping/{sample}_R2_host_unmapped.fastq",
       bai_unmapped=temp("02_HostMapping/{sample}_host_unmapped.bam.bai")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_host_mapping_extract_host_filtered_reads",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "logs/{sample}_host_mapping_extract_host_filtered_reads.log"
    benchmark: "logs/{sample}_host_mapping_extract_host_filtered_reads.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        samtools index {input.bam_unmapped}
        picard -Xmx8g SamToFastq I={input.bam_unmapped} F={output.reads1} F2={output.reads2} VALIDATION_STRINGENCY=SILENT
        """

rule assemble_host_mapped:
    input:
        reads1 = rules.host_mapping_extract_host_mapped_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_mapped_reads.output.reads2,
    output:
        scaffolds_unparsed = "05_Assembly/host_mapped/{sample}_scaffolds_unparsed.fasta",
        scaffolds = "05_Assembly/host_mapped/{sample}_scaffolds.fasta",
        graph = "05_Assembly/host_mapped/{sample}_assembly_graph.fastg",
        spades_log = "05_Assembly/host_mapped/{sample}_spades.log",
    params:
        outdir = lambda wildcards: "05_Assembly/host_mapped/"+wildcards.sample,
        filt_py = "scripts/parse_spades_metagenome.py",
        length_t = 1000,
        cov_t = 1,
        memory_limit = 350,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_assemble_host_mapped",
        account="pengel_spirit",
        runtime_s=convertToSec("0-20:00:00"),
    resources:
        mem_mb = convertToMb("200G"),
    retries: 3
    threads: 4
    log: "logs/{sample}_assemble_host_mapped.log"
    benchmark: "logs/{sample}_assemble_host_mapped.benchmark"
    conda: "envs/spades-env.yaml"
    shell:
        """
        if [ -f \"{params.outdir}/scaffolds.fasta\" ]; then
          echo {params.outdir}/scaffolds.fasta\" exists copying and cleaning\"
        else
          if [ -d \"{params.outdir}/\" ]; then
            echo {params.outdir}\" exists resuming spades\"
            spades.py --continue -o {params.outdir} || true
          else
            echo \"{params.outdir} not found. Starting new spades run.\"
            spades.py -m {params.memory_limit} --meta -1 {input.reads1} -2 {input.reads2} -t {threads} -o {params.outdir} || true
          fi
        fi
        if [ -f \"{params.outdir}/scaffolds.fasta\" ]; then
          cp {params.outdir}/scaffolds.fasta {output.scaffolds_unparsed}
        else
          echo \"assembly may have failed for \"{wildcards.sample}
          echo \"touching file \"{output.scaffolds_unparsed}
          touch {output.scaffolds_unparsed}
        fi
        python3 {params.filt_py} -i {output.scaffolds_unparsed} -o {output.scaffolds} -l {params.length_t} -c {params.cov_t}
        cp {params.outdir}/assembly_graph.fastg {output.graph}
        cp {params.outdir}/spades.log {output.spades_log}
        rm -rf {params.outdir}
        """

# use metaeuk to annotate scaffolds

rule microbiomedb_mapping:
    input:
        reads1=rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2=rules.host_mapping_extract_host_filtered_reads.output.reads2,
        index=multiext("database/"+DBs["microbiome"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        genome_db = os.path.join(os.getcwd(), "database/"+DBs["microbiome"]),
    output:
        bam_mapped=temp("03_MicrobiomeMapping/{sample}_microbiome_mapped.bam"),
        bam=temp("03_MicrobiomeMapping/{sample}_microbiome.bam"),
        sam=temp("03_MicrobiomeMapping/{sample}.sam"),
        sam_mapped=temp("03_MicrobiomeMapping/{sample}_microbiome_mapped.sam"),
        flagstat_03 = "03_MicrobiomeMapping/{sample}_flagstat.tsv",
    params:
        perl_extract_mapped = "scripts/filter_sam_aln_length.pl",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_microbiomedb_mapping",
        account="pengel_spirit",
        runtime_s=convertToSec("0-4:10:00"),
    resources:
        mem_mb = convertToMb("20G")
    threads: 8
    log: "logs/{sample}_microbiomedb_mapping.log"
    benchmark: "logs/{sample}_microbiomedb_mapping.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        bwa mem -t {threads} {input.genome_db} {input.reads1} {input.reads2} > {output.sam}
        samtools view -h -F4 -@ {threads} {output.sam} | samtools view - -F 0x800 -h | perl {params.perl_extract_mapped} - > {output.sam_mapped}
        samtools view -bh {output.sam_mapped} | samtools sort - > {output.bam_mapped}
        samtools view -bh {output.sam} | samtools sort - > {output.bam}
        samtools flagstat -O tsv {output.bam} > {output.flagstat_03}
        """

rule microbiomedb_extract_reads:
    input:
        bam_mapped=rules.microbiomedb_mapping.output.bam_mapped
    output:
        reads1=temp("03_MicrobiomeMapping/{sample}_R1_microbiome_mapped.fastq"),
        reads2=temp("03_MicrobiomeMapping/{sample}_R2_microbiome_mapped.fastq"),
        bai_mapped=temp("03_MicrobiomeMapping/{sample}_microbiome_mapped.bam.bai"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_microbiomedb_extract_reads",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "logs/{sample}_microbiomedb_extract_reads.log"
    benchmark: "logs/{sample}_microbiomedb_extract_reads.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        samtools index {input.bam_mapped}
        picard -Xmx8g SamToFastq I={input.bam_mapped} F={output.reads1} F2={output.reads2} VALIDATION_STRINGENCY=SILENT
        """

rule microbiomedb_direct_mapping:
    input:
        reads1=rules.trim.output.reads1,
        reads2=rules.trim.output.reads2,
        index=multiext("database/"+DBs["microbiome"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        genome_db = os.path.join(os.getcwd(), "database/"+DBs["microbiome"]),
    output:
        bam_mapped=temp("04_MicrobiomeMappingDirect/{sample}_microbiome_mapped.bam"),
        bam=temp("04_MicrobiomeMappingDirect/{sample}_microbiome.bam"),
        sam=temp("04_MicrobiomeMappingDirect/{sample}.sam"),
        sam_mapped=temp("04_MicrobiomeMappingDirect/{sample}_microbiome_mapped.sam"),
        flagstat_04 = "04_MicrobiomeMappingDirect/{sample}_flagstat.tsv"
    params:
        perl_extract_mapped = "scripts/filter_sam_aln_length.pl",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_microbiomedb_mapping",
        account="pengel_spirit",
        runtime_s=convertToSec("0-4:10:00"),
    resources:
        mem_mb = convertToMb("20G")
    threads: 8
    log: "logs/{sample}_microbiomedb_mapping.log"
    benchmark: "logs/{sample}_microbiomedb_mapping.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        bwa mem -t {threads} {input.genome_db} {input.reads1} {input.reads2} > {output.sam}
        samtools view -h -F4 -@ {threads} {output.sam} | samtools view - -F 0x800 -h | perl {params.perl_extract_mapped} - > {output.sam_mapped}
        samtools view -bh {output.sam_mapped} | samtools sort - > {output.bam_mapped}
        samtools view -bh {output.sam} | samtools sort - > {output.bam}
        samtools flagstat -O tsv {output.bam} > {output.flagstat_04}
        """

rule microbiomedb_direct_extract_reads:
    input:
        bam_mapped=rules.microbiomedb_direct_mapping.output.bam_mapped
    output:
        reads1=temp("04_MicrobiomeMappingDirect/{sample}_R1_microbiome_mapped.fastq"),
        reads2=temp("04_MicrobiomeMappingDirect/{sample}_R2_microbiome_mapped.fastq"),
        bai_mapped=temp("04_MicrobiomeMappingDirect/{sample}_microbiome_mapped.bam.bai"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_microbiomedb_direct__extract_reads",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "logs/{sample}_microbiomedb_extract_reads.log"
    benchmark: "logs/{sample}_microbiomedb_extract_reads.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        samtools index {input.bam_mapped}
        picard -Xmx8g SamToFastq I={input.bam_mapped} F={output.reads1} F2={output.reads2} VALIDATION_STRINGENCY=SILENT
        """

###############################
###############################
###############################
# MAGs
###############################
###############################
###############################

###############################
# MAGs - own pipeline
# explore the use of based on Shini lab's pipeline
# and others such as Anvio
###############################

rule assemble_host_unmapped:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
    output:
        scaffolds_unparsed = "05_Assembly/host_unmapped/{sample}_scaffolds_unparsed.fasta",
        scaffolds = "05_Assembly/host_unmapped/{sample}_scaffolds.fasta",
        graph = "05_Assembly/host_unmapped/{sample}_assembly_graph.fastg",
        spades_log = "05_Assembly/host_unmapped/{sample}_spades.log",
    params:
        outdir = lambda wildcards: "05_Assembly/host_unmapped/"+wildcards.sample,
        filt_py = "scripts/parse_spades_metagenome.py",
        length_t = 1000,
        cov_t = 1,
        memory_limit = lambda wildcards: "550" if wildcards.sample in ["M1.4", "D2.1", "D2.4", "D2.5", "F3.5"] else "200",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_assemble_host_unmapped",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("1-20:00:00") if wildcards.sample in ["M1.4", "D2.1", "D2.4", "D2.5", "F3.5"] else convertToSec("0-20:00:00"),
    resources:
        mem_mb = lambda wildcards: convertToMb("350G") if wildcards.sample in ["M1.4", "D2.1", "D2.5", "F3.5"] else ( convertToMb("450G") if wildcards.sample == "D2.4" else convertToMb("200G") ),
    retries: 3
    threads: 4
    log: "logs/{sample}_assemble_host_unmapped.log"
    benchmark: "logs/{sample}_assemble_host_unmapped.benchmark"
    conda: "envs/spades-env.yaml"
    shell:
        """
        if [ -f \"{params.outdir}/scaffolds.fasta\" ]; then
          echo {params.outdir}/scaffolds.fasta\" exists copying and cleaning\"
        else
          if [ -d \"{params.outdir}/\" ]; then
            echo {params.outdir}\" exists resuming spades\"
            spades.py --continue -o {params.outdir} || true
          else
            echo \"{params.outdir} not found. Starting new spades run.\"
            spades.py -m {params.memory_limit} --meta -1 {input.reads1} -2 {input.reads2} -t {threads} -o {params.outdir} || true
          fi
        fi
        if [ -f \"{params.outdir}/scaffolds.fasta\" ]; then
          cp {params.outdir}/scaffolds.fasta {output.scaffolds_unparsed}
        else
          echo \"assembly may have failed for \"{wildcards.sample}
          echo \"touching file \"{output.scaffolds_unparsed}
          touch {output.scaffolds_unparsed}
        fi
        python3 {params.filt_py} -i {output.scaffolds_unparsed} -o {output.scaffolds} -l {params.length_t} -c {params.cov_t}
        cp {params.outdir}/assembly_graph.fastg {output.graph}
        cp {params.outdir}/spades.log {output.spades_log}
        rm -rf {params.outdir}
        """

rule map_to_assembly:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        scaffolds = rules.assemble_host_unmapped.output.scaffolds
    output:
        assembly_mapped = temp("05_Assembly/MapToAssembly/{sample}.bam"),
        mapped = "05_Assembly/host_unmapped/check_assembly/{sample}_assembly_mapped_counts.txt",
        sam = temp("05_Assembly/host_unmapped/check_assembly/{sample}_assembly.sam"),
        bam = temp("05_Assembly/host_unmapped/check_assembly/{sample}_assembly.bam"),
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_map_to_assembly",
        account="pengel_spirit",
        runtime_s=convertToSec("1-00:00:00"),
    resources:
        mem_mb = 8000
    threads: 5
    log: "logs/{sample}_map_to_assembly.log"
    benchmark: "logs/{sample}_map_to_assembly.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        mkdir -p 05_Assembly/MapToAssembly
        bwa index {input.scaffolds}
        bwa mem -t {threads} {input.scaffolds} {input.reads1} {input.reads2} > {output.sam}
        samtools view -bh {output.sam} | samtools sort - > {output.bam}
        cp {output.bam} {output.assembly_mapped}
        samtools view -c -f 65 -F 2308 {output.bam} | sed \"s/^/{wildcards.sample} /\" >> {output.mapped};
        rm 05_Assembly/host_unmapped/{wildcards.sample}_scaffolds.fasta.amb
        rm 05_Assembly/host_unmapped/{wildcards.sample}_scaffolds.fasta.ann
        rm 05_Assembly/host_unmapped/{wildcards.sample}_scaffolds.fasta.bwt
        rm 05_Assembly/host_unmapped/{wildcards.sample}_scaffolds.fasta.pac
        rm 05_Assembly/host_unmapped/{wildcards.sample}_scaffolds.fasta.sa
        """

rule cat_and_clean_counts_assembly:
    input:
        mapped = expand("05_Assembly/host_unmapped/check_assembly/{sample}_assembly_mapped_counts.txt", sample=SAMPLES)
    output:
        mapped = "05_Assembly/host_unmapped/check_assembly/Assembly_mapped_counts.txt"
    log:
        "logs/cat_and_clean_counts_assembly.log"
    shell:
        """
        cat {input.mapped} > {output.mapped}
        """

rule summarize_mapping_assembly:
    input:
        reads1 = expand("02_HostMapping/{sample}_R1_host_unmapped.fastq", sample=SAMPLES),
        scaffolds = expand("05_Assembly/host_unmapped/{sample}_scaffolds.fasta", sample=SAMPLES),
        scaffolds_unparsed = expand("05_Assembly/host_unmapped/{sample}_scaffolds_unparsed.fasta", sample=SAMPLES),
        counts = "05_Assembly/host_unmapped/check_assembly/Assembly_mapped_counts.txt"
    output:
        summary_assembly = "05_Assembly/host_unmapped/check_assembly/Assembly_mapping_summary.csv",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="summarize_mapping_assembly",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 5
    log: "logs/summarize_mapping_assembly.log"
    benchmark: "logs/summarize_mapping_assembly.benchmark"
    run:
        mapped_dict = make_count_dict("05_Assembly/host_unmapped/check_assembly/Assembly_mapped_counts.txt")

        Sample_list = list()
        MinContigLen_list = list()
        MaxContigLen_list = list()
        NumberOfContigs_list = list()
        NumberOfContigsUn_list = list()
        AssemblySize_list = list()
        NumReads_list = list()
        AssemblyMapped_list = list()
        ProportionMapped_list = list()

        for i in range(len(input.reads1)):
            print(i)
            Sample = input.reads1[i].split("/")[1].split("_R1")[0]
            print(Sample)
            Sample_list.append(Sample)

            MinContigLen = int(os.popen("echo $(cat "+input.scaffolds[i]+" | grep -c \'>\' | cut -d\'_\' -f4 | sort -n | head -1 )").read())
            MinContigLen_list.append(MinContigLen)

            MaxContigLen = int(os.popen("echo $(cat "+input.scaffolds[i]+" | grep -c \'>\' | cut -d\'_\' -f4 | sort -n | tail -1 )").read())
            MaxContigLen_list.append(MaxContigLen)

            NumberOfContigs = int(os.popen("echo $(cat "+input.scaffolds[i]+" | grep -c \'>\')").read())
            NumberOfContigs_list.append(NumberOfContigs)

            NumberOfContigsUn = int(os.popen("echo $(cat "+input.scaffolds_unparsed[i]+" | grep -c \'>\')").read())
            NumberOfContigsUn_list.append(NumberOfContigsUn)

            AssemblySize = int(os.popen("echo $(cat "+input.scaffolds[i]+" | grep -v \'>\' | tr -d \'\n\' | wc -m)").read())
            AssemblySize_list.append(AssemblySize)

            NumReads = int(os.popen("echo $(cat "+input.reads1[i]+" | wc -l)").read())/4
            NumReads_list.append(NumReads)

            AssemblyMapped = mapped_dict[Sample]
            AssemblyMapped_list.append(AssemblyMapped)

        ProportionMapped_list = [str(round(x/y*100, 2)) for (x,y) in zip(AssemblyMapped_list, NumReads_list)]

        Summary_path = os.path.join(os.getcwd(), output.summary_assembly)

        with open(Summary_path, "w") as file:
            file.write("Sample, Assembly size, Number of reads, " +
            "Total nummber of scaffolds, Number of filtered scaffolds," +
            "Min contig length, Max contig length," +
            "Number mapped, Percent mapped\n")
            for i in range(len(input.reads1)):
                file.write(f"{Sample_list[i]}, {AssemblySize_list[i]}, {NumReads_list[i]}, {NumberOfContigsUn_list[i]}, {NumberOfContigs_list[i]}, {MinContigLen_list[i]}, {MaxContigLen_list[i]}, {AssemblyMapped_list[i]}, {ProportionMapped_list[i]}\n")

rule backmapping:
    input:
        assembly = lambda wildcards: expand(rules.assemble_host_unmapped.output.scaffolds, sample=wildcards.sample),
        reads1 = expand("02_HostMapping/{sample}_R1_host_unmapped.fastq", sample=SAMPLES),
        reads2 = expand("02_HostMapping/{sample}_R2_host_unmapped.fastq", sample=SAMPLES),
    output:
        depth_files = expand("06_MAG_binning/backmapping/{{sample}}/{each_sample}_mapped_to_{{sample}}.depth", each_sample=SAMPLES),
    params:
        samples = SAMPLES,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="backmapping_{sample}",
        account="pengel_spirit",
        runtime_s=convertToSec("1-10:00:00"),
    resources:
        mem_mb = convertToMb("10G")
    threads: 4
    conda: "envs/mapping-env.yaml"
    log: "logs/{sample}_backmapping.log"
    benchmark: "logs/{sample}_backmapping.benchmark"
    shell:
        """
        assembly={wildcards.sample}
        bwa index {input.assembly}
        for sample in {params.samples};
        do
            bwa mem -t {threads} {input.assembly} 02_HostMapping/${{sample}}_R1_host_unmapped.fastq 02_HostMapping/${{sample}}_R1_host_unmapped.fastq > 06_MAG_binning/backmapping/${{assembly}}/${{sample}}_mapped_to_${{assembly}}.sam;
            echo "Mapping complete for ${{sample}} to ${{assembly}}"
            samtools view -bh 06_MAG_binning/backmapping/${{assembly}}/${{sample}}_mapped_to_${{assembly}}.sam | samtools sort - > 06_MAG_binning/backmapping/${{assembly}}/${{sample}}_mapped_to_${{assembly}}.bam
            echo "Deleting sam file for for ${{sample}} to ${{assembly}}"
            rm 06_MAG_binning/backmapping/${{assembly}}/${{sample}}_mapped_to_${{assembly}}.sam
            echo "Making depth file for for ${{sample}} to ${{assembly}}"
            export OMP_NUM_THREADS={threads}
            jgi_summarize_bam_contig_depths --outputDepth 06_MAG_binning/backmapping/${{assembly}}/${{sample}}_mapped_to_${{assembly}}.depth 06_MAG_binning/backmapping/${{assembly}}/${{sample}}_mapped_to_${{assembly}}.bam
            rm -f 06_MAG_binning/backmapping/${{assembly}}/${{sample}}_mapped_to_${{assembly}}.bam
        done
        """

rule merge_depths:
    input:
        depth_files = expand("06_MAG_binning/backmapping/{{sample}}/{each_sample}_mapped_to_{{sample}}.depth", each_sample=SAMPLES)
    output:
        depth_file_merged = "06_MAG_binning/backmapping/{sample}/{sample}_merged.depth"
    params:
        samples = SAMPLES,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="give_name",
        account="pengel_spirit",
        runtime_s=convertToSec("0-01:00:00"),
    resources:
        mem_mb = convertToMb("4G")
    threads: 4
    conda: "envs/mapping-env.yaml"
    log: "logs/{sample}_merge_depths.log"
    benchmark: "logs/{sample}_merge_depths.benchmark"
    shell:
        """
        scripts/merge_depths.pl {input.depth_files} > {output.depth_file_merged}
        """

rule binning:
    input:
        assembly = rules.assemble_host_unmapped.output.scaffolds,
        depth_file_merged = "06_MAG_binning/backmapping/{sample}/{sample}_merged.depth",
    output:
        bins = directory("06_MAG_binning/bins/{sample}-metabat2/"),
    params:
        min_contig_size=2500, # Metabat2 default
        min_bin_size=200000, # Metabat2 default
        max_edges=200, # Metabat2 default
        min_cv=1, # Metabat2 default
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="give_name",
        account="pengel_spirit",
        runtime_s=convertToSec("0-01:00:00"),
    resources:
        mem_mb = convertToMb("2G")
    threads: 16
    conda: "envs/mapping-env.yaml"
    log: "logs/{sample}_binning.log"
    benchmark: "logs/{sample}_binning.benchmark"
    shell:
        # metabat2 -i {input.assembly} -a {input.depth_file_merged} -o {params.prefix} --minContig {params.min_contig_size} --maxEdges {params.max_edges} -x {params.min_cv} --numThreads {threads} --minClsSize {params.min_bin_length} --saveCls
        """
        metabat2 -i {input.assembly} -a {input.depth_file_merged} -o 06_MAG_binning/bins/{wildcards.sample}-metabat2/MAG --minContig {params.min_contig_size} --maxEdges {params.max_edges} -x {params.min_cv} --numThreads {threads}
        """


###############################
###############################
# DREP
# dRep check_dependencies
# Near Essential
#
# Mash - Makes primary clusters (v1.1.1 confirmed works)
# MUMmer - Performs default ANIm comparison method (v3.23 confirmed works)
# Recommended
#
# fastANI - A fast secondary clustering algorithm
# CheckM - Determines contamination and completeness of genomes (v1.0.7 confirmed works)
# gANI (aka ANIcalculator) - Performs gANI comparison method (v1.0 confirmed works)
# Prodigal - Used be both checkM and gANI (v2.6.3 confirmed works)
# Accessory
#
# NSimScan - Only needed for goANI algorithm
# Centrifuge - Deprecated; not used anymore

###############################
###############################

rule process_metabat2:
    input:
        sample_bins = "06_MAG_binning/bins/{sample}-metabat2/",
    output:
        all_mags = directory("06_MAG_binning/bins_renamed/{sample}"),
    log: "logs/{sample}_process_metabat2.log"
    benchmark: "logs/{sample}_process_metabat2.benchmark"
    params:
        outpath = "06_MAG_binning/bins_renamed/"
    shell:
        """
        mkdir -p {output.all_mags}
        dir={input.sample_bins}
        echo "getting mags from "${{dir}}
        sample={wildcards.sample}
        mags_list=$(ls ${{dir}})
        for mag in ${{mags_list}}
        do
            mag_name=${{mag/MAG./}}
            mag_num=${{mag_name%%.fa}}
            cp ${{dir}}/${{mag}} {output.all_mags}/MAG_${{sample}}_${{mag_num}}.fa
            ln -s {output.all_mags}/MAG_${{sample}}_${{mag_num}}.fa {params.outpath}//MAG_${{sample}}_${{mag_num}}.fa
        done
        """

rule summarize_metabat2_contig_fates:
    input:
        sample_bins = "06_MAG_binning/bins_renamed/{sample}/",
        scaffolds_unparsed = "05_Assembly/host_unmapped/{sample}_scaffolds_unparsed.fasta",
        scaffolds = "05_Assembly/host_unmapped/{sample}_scaffolds.fasta",
    output:
        contig_fates = "06_MAG_binning/contig_fates/{sample}_contig_fates.csv",
    log: "logs/{sample}_summarize_metabat2_contig_fates"
    benchmark: "logs/{sample}_summarize_metabat2_contig_fabenchmark"
    run:
        with open(output.contig_fates, "w") as out_fh:
            out_fh.write(f"sample,contig_name,length,coverage,passed_filter,binned,bin_name\n")
            scaffolds_dict = {}
            with open(input.scaffolds_unparsed, "r") as unparsed_fh:
                for line in unparsed_fh:
                    if not line.startswith(">"):
                        pass
                    else:
                        line = line.strip()
                        contig = line.split(">")[1]
                        if contig not in scaffolds_dict.keys():
                            scaffolds_dict[contig] = ["F", "N", "NA"]
            with open(input.scaffolds, "r") as parsed_fh:
                for line in parsed_fh:
                    if not line.startswith(">"):
                        pass
                    else:
                        line = line.strip()
                        contig = line.split(">")[1]
                        scaffolds_dict[contig] = ["P", "N", "NA"]
            for mag_file in [x for x in os.listdir(input.sample_bins) if x.endswith(".fa")]:
                with open(os.path.join(input.sample_bins, mag_file), "r") as mag_fh:
                    for line in mag_fh:
                        if not line.startswith(">"):
                            pass
                        else:
                            line = line.strip()
                            contig = line.split(">")[1]
                            scaffolds_dict[contig] = ["P", "Y", mag_file.split(".fa")[0]]
            for contig_name in scaffolds_dict.keys():
                out_list = scaffolds_dict[contig_name]
                length = contig_name.split("_")[3]
                cov = contig_name.split("_")[5]
                out_fh.write(f"{wildcards.sample},{contig_name},{length},{cov},{out_list[0]},{out_list[1]},{out_list[2]}\n")

rule summarize_metabat2_contig_coverages:
    input:
        contig_fates = "06_MAG_binning/contig_fates/{sample}_contig_fates.csv",
        depth_file_merged = "06_MAG_binning/backmapping/{sample}/{sample}_merged.depth",
        all_mags = "06_MAG_binning/bins_renamed/{sample}",
    output:
        sample_contig_coverages = "06_MAG_binning/contig_fates/backmapping_coverages/{sample}_contig_coverages.csv",
    params:
        contig_coverages_per_bin_dir="06_MAG_binning/contig_fates/backmapping_coverages/contig_coverage_by_bin",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="give_name",
        account="pengel_spirit",
        runtime_s=convertToSec("0-01:00:00"),
    resources:
        mem_mb = convertToMb("4G")
    threads: 8
    log: "logs/{sample}_summarize_metabat2_contig_coverages.log"
    benchmark: "logs/{sample}_summarize_metabat2_contig_coverages.benchmark"
    run:
        sample = wildcards.sample
        contig_bin_map = {}
        with open(input.contig_fates, "r") as in_fh:
            is_header = True
            for line in in_fh:
                line = line.strip()
                if is_header:
                    is_header = False
                else:
                    line = line.strip()
                    # sample = line.split(",")[0]
                    contig = line.split(",")[1]
                    # length = line.split(",")[2]
                    bin = line.split(",")[6]
                    contig_bin_map[contig] = bin
        """
        list of all scaffolds and coverages
        """
        with open(output.sample_contig_coverages, "w") as out_fh:
            out_fh.write(f"bin_name,contig,avg_coverage,variance,sample\n")
            # for sample in samples_list:
            with open("06_MAG_binning/backmapping/"+sample+"/"+sample+"_merged.depth", "r") as depths_fh:
                is_header = True
                for line in depths_fh:
                    line = line.strip()
                    line_split = line.split("\t")
                    if is_header:
                        header = line_split
                        is_header = False
                    else:
                        contig = line_split[0]
                        coverages_dict = {key.split("_mapped_to_")[0]: val for key, val in zip(header, line_split) if key.endswith(".bam")}
                        vars_dict = {key.split("_mapped_to_")[0]: val for key, val in zip(header, line_split) if key.endswith(".bam-var")}
                        for sample, coverage in coverages_dict.items():
                            out_fh.write(f"{contig_bin_map[contig]},{contig},{coverage},{vars_dict[sample]},{sample}\n")

rule checkm_evaluation:
    input:
        bin = "06_MAG_binning/bins_renamed/{sample}"
    output:
        checkm_summary = "06_MAG_binning/evaluate_bins/{sample}_checkm.summary",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="give_name",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:00:00"),
    resources:
        mem_mb = convertToMb("50G") # needs at least roughly 40G
    threads: 8
    conda: "envs/mags-env.yaml"
    log: "logs/{sample}_checkm_evaluation.log"
    benchmark: "logs/{sample}_checkm_evaluation.benchmark"
    shell:
        """
        out_file={output.checkm_summary}
        checkm lineage_wf -x fa {input.bin} ${{out_file/_checkm.summary/}} --threads {threads} -f {output.checkm_summary} --tab_table
        """

rule prepare_info_for_drep:
    input:
        checkm_summary = expand("06_MAG_binning/evaluate_bins/{sample}_checkm.summary", sample=SAMPLES),
    output:
        compiled_checkm_summary = "06_MAG_binning/evaluate_bins/checkm_drep_summary.txt",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="give_name",
        account="pengel_spirit",
        runtime_s=convertToSec("0-04:00:00"),
    resources:
        mem_mb = convertToMb("2G")
    threads: 4
    log: "logs/prepare_info_for_drep.log"
    benchmark: "logs/prepare_info_for_drep.benchmark"
    run:
        with open(output.compiled_checkm_summary, "w") as out_fh:
            out_fh.write("genome,completeness,contamination\n")
            for file in input.checkm_summary:
                with open(file, "r") as in_fh:
                    for line in in_fh:
                        if line.startswith("MAG_"):
                            genome = line.split()[0]
                            completeness = line.split()[12]
                            contamination = line.split()[13]
                            out_fh.write(f"{genome}.fa,{completeness},{contamination}\n")
                        else:
                            continue

rule drep:
    input:
        checkm_summary = expand("06_MAG_binning/evaluate_bins/{sample}_checkm.summary", sample=SAMPLES),
        bins = expand("06_MAG_binning/bins_renamed/{sample}", sample=SAMPLES),
        compiled_checkm_summary = "06_MAG_binning/evaluate_bins/checkm_drep_summary.txt",
    output:
        drep_S = "06_MAG_binning/drep_results/data_tables/Sdb.csv",
        drep_N = "06_MAG_binning/drep_results/data_tables/Ndb.csv",
        drep_Wi = "06_MAG_binning/drep_results/data_tables/Widb.csv",
        drep_gI = "06_MAG_binning/drep_results/data_tables/genomeInformation.csv",
    params:
        bins = expand("06_MAG_binning/bins_renamed/{sample}/*.fa", sample=SAMPLES),
        overlap = 0.2, # ask Lucas why
        ani = 0.95,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="give_name",
        account="pengel_spirit",
        runtime_s=convertToSec("0-04:00:00"),
    resources:
        mem_mb = convertToMb("5G")
    threads: 24
    conda: "envs/mags-env.yaml"
    log: "logs/drep.log"
    benchmark: "logs/drep.benchmark"
    shell:
        """
        out_file={output.drep_S}
        dRep dereplicate ${{out_file/\/data_tables\/Sdb.csv/}} -g {params.bins} --genomeInfo {input.compiled_checkm_summary} -comp 0 -con 1000 -sa {params.ani} -nc {params.overlap} -p {threads}
        """

rule extract_mag_lists:
    input:
        checkm_summary = expand("06_MAG_binning/evaluate_bins/{sample}_checkm.summary", sample=SAMPLES),
        drep_C = "06_MAG_binning/drep_results/data_tables/Cdb.csv",
        drep_W = "06_MAG_binning/drep_results/data_tables/Wdb.csv",
        drep_gI = "06_MAG_binning/drep_results/data_tables/genomeInformation.csv"
    output:
        all_list = "06_MAG_binning/all_genomes.csv",
        high_quality = "06_MAG_binning/all_hq_genomes.csv",
        med_quality = "06_MAG_binning/all_medq_genomes.csv",
    log: "logs/extract_mag_list.log"
    benchmark: "logs/extract_mag_list.benchmark"
    run:
        write_header = False
        with open(output.all_list, "w+") as all_fh:
            for summary in input.checkm_summary:
                with open(summary, "r") as summary_fh:
                    if write_header == False:
                        print("reading summary_fh")
                    for line in summary_fh:
                        if write_header == False:
                            all_fh.write(f"{line}")
                            write_header = True
                        if line.startswith("MAG_"):
                            all_fh.write(f"{line}")
        with open(output.all_list, "r") as all_fh:
            print("reading all_fh")
            lengths = {}
            clusters = {}
            with open(input.drep_W, "r") as cluster_info_fh:
                print("reading from "+input.drep_W)
                for line in cluster_info_fh:
                    if line.startswith("MAG_"):
                        genome = line.split(",")[0].split(".fa")[0]
                        cluster = line.split(",")[1]
                        clusters[genome] = cluster
                with open(input.drep_gI, "r") as length_info_fh:
                    print("reading from "+input.drep_gI)
                    for line in length_info_fh:
                        line = line.strip()
                        if line.startswith("MAG_"):
                            genome = line.split(",")[0].split(".fa")[0]
                            length = line.split(",")[3]
                            N50 = line.split(",")[4]
                            lengths[genome] = [length, N50]
                    with open(output.high_quality, "w+") as hq_fh:
                        hq_fh.write(f"genome, lineage, completeness, contamination, length, N50, cluster\n")
                        with open(output.med_quality, "w+") as mq_fh:
                            mq_fh.write(f"genome, lineage, completeness, contamination, length, N50, cluster\n")
                            for line in all_fh:
                                if line.startswith("MAG_"):
                                    print(line)
                                    genome = line.split()[0]
                                    lineage = line.split()[1]
                                    completeness = line.split()[12]
                                    contamination = line.split()[13]
                                    g_length = lengths[genome][0]
                                    N50 = lengths[genome][1]
                                    print(f"{clusters}")
                                    if genome in clusters.keys():
                                        cluster = clusters[genome]
                                    else:
                                        clusters[genome] = "NA"
                                    if float(completeness) >= 95 and float(contamination) <= 5:
                                        hq_fh.write(f"{genome}, {lineage}, {completeness}, {contamination}, {g_length}, {N50}, {cluster}\n")
                                    if float(completeness) >= 50 and float(contamination) <= 10:
                                        mq_fh.write(f"{genome}, {lineage}, {completeness}, {contamination}, {g_length}, {N50}, {cluster}\n")

rule gtdb_annotate:
    input:
        genomes = expand("06_MAG_binning/bins_renamed/{sample}", sample=SAMPLES),
        all_list = "06_MAG_binning/all_genomes.csv",
    output:
        tax_info = "06_MAG_binning/gtdbtk_out_dir/classify/gtdbtk.bac120.summary.tsv",
    params:
        path_to_db="/work/FAC/FBM/DMF/pengel/spirit/aprasad/gtdb/release207_v2/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="give_name",
        account="pengel_spirit",
        runtime_s=convertToSec("2-00:00:00"),
    resources:
        mem_mb = convertToMb("150G") # needs at least roughly 40G
    threads: 8
    conda: "envs/gtdb-env.yaml"
    log: "logs/gtdb_annotate.log"
    benchmark: "logs/gtdb_annotate.benchmark"
    shell:
        """
        # set up path to database
        export GTDBTK_DATA_PATH={params.path_to_db}
        # find a more snakemakey way to do this
        ln -sf {PROJECT_PATH}/06_MAG_binning/bins_renamed/*/*.fa 06_MAG_binning/bins_renamed/
        genomes_dir="06_MAG_binning/bins_renamed/"
        out_file={output.tax_info}
        gtdbtk classify_wf --genome_dir ${{genomes_dir}} --out_dir ${{out_file/\/classify\/gtdbtk.bac120.summary.tsv/}} --extension ".fa" --write_single_copy_genes --keep_intermediates
        """

checkpoint make_phylo_table:
    input:
        drep_res = "06_MAG_binning/drep_results/data_tables/Cdb.csv",
        isolates = "config/IsolateGenomeInfo.csv",
        MAG_info = "06_MAG_binning/all_genomes.csv",
        tax_info = "06_MAG_binning/gtdbtk_out_dir/classify/gtdbtk.bac120.summary.tsv"
        # move the family - phylotype dictonary to config or make it a separate input file later.
    output:
        out_all = "06_MAG_binning/all_GenomeInfo_auto.tsv",
        out_tree = "06_MAG_binning/ForTree_GenomeInfo_auto.tsv",
        out_mags = "06_MAG_binning/MAGs_GenomeInfo_auto.tsv",
        out_mags_filt = "06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv",
        out_t_all = temp("06_MAG_binning/all_GenomeInfo_auto.csv"),
        out_t_tree = temp("06_MAG_binning/ForTree_GenomeInfo_auto.csv"),
        out_t_mags = temp("06_MAG_binning/MAGs_GenomeInfo_auto.csv"),
        out_t_mags_filt = temp("06_MAG_binning/MAGs_filt_GenomeInfo_auto.csv"),
    log: "logs/make_phylo_table.log"
    benchmark: "logs/make_phylo_table.benchmark"
    conda: "envs/rmd-env.yaml"
    threads: 2
    shell:
        """
        which R
        ./scripts/make_phylo_table.R {PROJECT_PATH}
        ./scripts/csv_to_tsv.py {output.out_t_all}
        ./scripts/csv_to_tsv.py {output.out_t_tree}
        ./scripts/csv_to_tsv.py {output.out_t_mags}
        ./scripts/csv_to_tsv.py {output.out_t_mags_filt}
        """

rule prepare_genomes:
    input:
        mags = "06_MAG_binning/bins_renamed"
    output:
        genome = "07_AnnotationAndPhylogenies/00_genomes/{genome}.fa",
    threads: 4
    params:
        info = "config/IsolateGenomeInfo.csv",
        sample_name = lambda wildcards: "_".join(wildcards.genome.split("MAG_")[1].split("_")[:-1]),
        mags_dir = "06_MAG_binning/bins_renamed/",
        ftp_summary = "https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt",
        assembly_summary_genbank = "assembly_summary_genbank.txt",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
    log: "logs/{genome}_prepare_genomes.log"
    benchmark: "logs/{genome}_prepare_genomes.benchmark"
    shell:
        """
        if [[ \"{wildcards.genome}\" == *\"MAG_\"* ]]; then
            echo \"{wildcards.genome} is a MAG.\"
            echo \"Copying {wildcards.genome}.fa from {params.mags_dir}/{params.sample_name}\"
            cp {params.mags_dir}/{params.sample_name}/{wildcards.genome}.fa 07_AnnotationAndPhylogenies/00_genomes/{wildcards.genome}.fa
        else
            if [ -f {output.genome} ]; then
                # if snakemake is working, this part should never run!
                echo \"{output.genome} exists.\" 2>&1 | tee -a {log}
            else
                echo \"cannot find {output.genome}. Going to try to download!\" 2>&1 | tee -a {log}
                wget -N \"{params.ftp_summary}\" 2>&1 | tee -a {log}
                strain=$(cat {params.info} | grep \"{wildcards.genome}\" | cut -f4 -d\",\" )
                number_strain_matches=$(cat {params.assembly_summary_genbank} | grep \"${{strain}}\" | wc -l)
                echo \"The strain is ${{strain}} --- 1\" 2>&1 | tee -a {log}
                if [[ ${{number_strain_matches}} != 1 ]]; then
                    echo "Strain name ${{strain}} matches multiple entries"
                    genbank_acc=$(cat {params.info} | grep \"{wildcards.genome},\" | cut -f2 -d\",\" )
                    if [[ ${{genbank_acc}} == \"NA\" ]]; then
                        echo "Genbank path for ${{strain}} is not available ${{genbank_acc}}"
                        genbank_path=$(cat {params.assembly_summary_genbank} | grep \"${{strain}}\" | cut -f15,20 | sort -k1 -r | head -1 | cut -f2)
                    else
                        genbank_path=$(cat {params.assembly_summary_genbank} | grep \"${{strain}}\" | cut -f1,15,20 | grep \"${{genbank_acc}}\" | cut -f3)
                        echo \"${{genbank_path}}_here\" 2>&1 | tee -a {log}
                    fi
                else
                    echo \"The strain is ${{strain}} --- 2\" 2>&1 | tee -a {log}
                    genbank_path=$(cat {params.assembly_summary_genbank} | grep \"${{strain}}\" | cut -f15,20 | sort -k1 -r | head -1 | cut -f2)
                fi
                    echo \"The strain is ${{strain}}\" 2>&1 | tee -a {log}
                    echo \"Downloading from ${{genbank_path}}/$(basename ${{genbank_path}})_genomic.fna.gz\" 2>&1 | tee -a {log}
                    wget -O {output.genome}.gz \"${{genbank_path}}/$(basename ${{genbank_path}})_genomic.fna.gz\" 2>&1 | tee -a {log}
                    ls {output.genome}.gz
                    echo \"Downloaded! unzipping {output.genome}.gz\" 2>&1 | tee -a {log}
                    gunzip < {output.genome}.gz > {output.genome} 2>&1 | tee -a {log}
                    echo \"Unzipped! Editing header to contain only ID of the genome\" 2>&1 | tee -a {log}
                    echo \">{wildcards.genome}\" > {output.genome}.temp
                    cat {output.genome} | grep -v \">\" >> {output.genome}.temp
                    mv {output.genome}.temp {output.genome}
            fi
        fi
        fz_file=\"07_AnnotationAndPhylogenies/00_genomes/{wildcards.genome}.fa.gz\"
        if [ -f ${{fz_file}} ]; then
            rm ${{fz_file}}
        fi
        """

rule annotate:
    input:
        genome = "07_AnnotationAndPhylogenies/00_genomes/{genome}.fa"
    output:
        faa = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa",
        ffn = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.ffn",
        gff = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.gff",
        fna = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.fna",
    params:
        outdir = "07_AnnotationAndPhylogenies/01_prokka/{genome}/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-4:10:00"),
        # jobname="annotate_{genome}"
    resources:
        mem_mb = 8000
    threads: 8
    log: "logs/{genome}_annotate.log"
    benchmark: "logs/{genome}_annotate.benchmark"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        prokka --compliant --force \
            --outdir {params.outdir} \
            --locustag {wildcards.genome} \
            --prefix {wildcards.genome} \
            --evalue 0.001 \
            {input.genome} 2>&1 | tee -a {log}
        """

rule prepare_faa:
    input:
        info = lambda wildcards: checkpoints.make_phylo_table.get().output.out_tree,
        prokka_faa_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa"
    output:
        faa_file = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/{genome}.faa"
    params:
        prokka_path="07_AnnotationAndPhylogenies/01_prokka",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit"
    log: "logs/prepare_faa_{genome}_{group}.log"
    shell:
        """
        cp {input.prokka_faa_file} {output.faa_file}
        """

rule run_orthofinder_phylo:
    input:
        faa_files = lambda wildcards: expand("07_AnnotationAndPhylogenies/02_orthofinder/{{group}}/{genome}.faa", genome=get_g_dict_for_defined_groups(checkpoints.make_phylo_table.get().output.out_tree)[wildcards.group]),
    output:
        orthogroups = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.txt",
        orthogroup_counts = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.GeneCount.tsv",
        orthogroup_stats = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/Results_{group}/Comparative_Genomics_Statistics/Statistics_Overall.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("1-2:10:00") if wildcards.group == "g__Lactobacillus" or wildcards.group == "g__Bifidobacterium" else convertToSec("0-5:10:00"),
        # jobname="run_orthofinder_phylo_{group}"
    resources:
         mem_mb = convertToMb("20G")
    threads: 8
    log: "logs/{group}_run_orthofinder_phylo.log"
    benchmark: "logs/{group}_run_orthofinder_phylo.benchmark"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        rm -rf $(dirname {input.faa_files[0]})/OrthoFinder/Results_{wildcards.group}
        orthofinder -og -t {threads} -n {wildcards.group} -f $(dirname {input.faa_files[0]})  2>&1 | tee -a {log}
        """

rule summarise_orthogroups:
    input:
        ortho_file = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.txt",
        genomes_list = lambda wildcards: checkpoints.make_phylo_table.get().output.out_tree,
    output:
        summary_orthogroups = "07_AnnotationAndPhylogenies/02_orthofinder/{group}_Orthogroups_summary.csv"
    params:
        group = lambda wildcards: wildcards.group,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
    log: "logs/{group}_summarise_orthogroups.log"
    benchmark: "logs/{group}_summarise_orthogroups.benckmark"
    script:
        "scripts/summarise_orthogroups.py"

rule get_single_ortho_phylo:
    input:
        ortho_file = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.txt",
        genomes_list = lambda wildcards: checkpoints.make_phylo_table.get().output.out_tree,
    output:
        single_ortho = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/{group}_single_ortho.txt",
        single_ortho_MAGs = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/{group}_single_ortho_MAGs.txt",
    params:
        group = lambda wildcards: wildcards.group,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="{group}_extract_orthologs",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{group}_get_single_ortho_phylo.log"
    benchmark: "logs/{group}_get_single_ortho_phylo.benchmark"
    script:
        "scripts/get_single_ortho_phylo.py"


# eventually use this instead of single ortho (if present in more than half the MAGs cut-off) for core genes
# rule motupan:
#     input:
#         fnas = lambda wildcards: expand("07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.fna", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags, full_list = True)),
#         checkm_files = expand("06_MAG_binning/evaluate_bins/{sample}_checkm.summary", sample=SAMPLES),
#     output:
#         outfile = "06_MAG_binning/mOTUlizer/mOTUlizer_output.tsv",
#         checkm_concat = "06_MAG_binning/mOTUlizer/checkm_concat.tsv",
#         ani_parsed = "06_MAG_binning/mOTUlizer/ani_parsed.tsv"
#     params:
#         annotations_dir = "07_AnnotationAndPhylogenies/01_prokka/",
#         seed_completeness = 50,
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-4:10:00"),
#         # jobname="annotate_{genome}"
#     resources:
#         mem_mb = 8000
#     threads: 8
#     log: "logs/motupan.log"
#     benchmark: "logs/motupan.benchmark"
#     conda: "envs/mags-env.yaml"
#     shell:
#         """
#         echo -e \"query\tsubject\tani\" > {output.ani_parsed}
#         awk -F\',\' \'{{ print($2,"\t",$1,"\t",$3) }}\' {input.ani_file} | sed -e \'s/ \t /\t/g\' | grep -v \"ani\" >> {output.ani_parsed}
#         head -1 {input.checkm_files} > {output.checkm_concat}
#         cat {input.checkm_files} | grep -v \"Bin Id\" >> {output.checkm_concat}
#         mOTUlize.py --fnas {params.annotations_dir}/*/*.fna -o {output.outfile} --checkm {output.checkm_concat} --similarities {output.ani_parsed} --MAG-completeness {params.seed_completeness} --prefix \"magOTU_\"
#         """

rule extract_orthologs_phylo:
    input:
        ortho_single = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/{group}_single_ortho.txt",
        faa_files = lambda wildcards: expand("07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa", genome=get_g_dict_for_defined_groups(checkpoints.make_phylo_table.get().output.out_tree)[wildcards.group]),
        ffn_files = lambda wildcards: expand("07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.ffn", genome=get_g_dict_for_defined_groups(checkpoints.make_phylo_table.get().output.out_tree)[wildcards.group]),
    output:
        ortho_seq_dir = directory("07_AnnotationAndPhylogenies/02_orthofinder/{group}/single_ortholog_sequences/"),
        done = touch("07_AnnotationAndPhylogenies/02_orthofinder/{group}/single_ortholog_sequences.done")
    params:
        # to the prefix, script adds, "genome/genome.xxx"
        faaffndir = "07_AnnotationAndPhylogenies/01_prokka/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname=lambda wildcards: wildcards.group+"_extract_orthologs",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{group}_extract_orthologs_phylo.log"
    benchmark: "logs/{group}_extract_orthologs_phylo.benchmark"
    conda: "envs/phylogenies-env.yaml"
    threads: 5
    shell:
        """
        python3 scripts/extract_orthologs_phylo.py --orthofile {input.ortho_single} --outdir {output.ortho_seq_dir} --extracted_ffndir {params.faaffndir}
        """

rule align_orthologs:
    input:
        orthogroups_sequences = "07_AnnotationAndPhylogenies/02_orthofinder/{group}/single_ortholog_sequences/",
    output:
        out_dir = directory("07_AnnotationAndPhylogenies/03_aligned_orthogroups/{group}/"),
        done = touch("07_AnnotationAndPhylogenies/03_aligned_orthogroups/{group}/mafft.done"),
    conda: "envs/phylogenies-env.yaml"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-18:10:00"),
        jobname=lambda wildcards:"align_orthologs_"+wildcards.group
    resources:
        mem_mb = convertToMb("20G")
    threads: 8
    log: "logs/{group}_align_orthologs.log"
    benchmark: "logs/{group}_align_orthologs.benchmark"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        mkdir -p {output.out_dir}
        for OG in $(ls {input.orthogroups_sequences})
        do
            echo "starting alignment for ${{OG}}" 2>&1 | tee -a {log}
            mafft --amino --inputorder --localpair --maxiterate 1000 {input.orthogroups_sequences}/${{OG}} > {output.out_dir}/${{OG/.fa/_aligned.fa}}
        done
        touch {output.done}
        """
# no pruning, just filling
rule prune_and_concat:
    input:
        aligned_dir = "07_AnnotationAndPhylogenies/03_aligned_orthogroups/{group}/",
        done = "07_AnnotationAndPhylogenies/03_aligned_orthogroups/{group}/mafft.done",
        genomes_list = lambda wildcards: checkpoints.make_phylo_table.get().output.out_tree,
    output:
        pruned_dir = directory("07_AnnotationAndPhylogenies/04_pruned_and_concat_alignments/{group}/"),
        pruned_cat = "07_AnnotationAndPhylogenies/04_pruned_and_concat_alignments/{group}/CoreGeneAlignment.fasta",
    params:
        pipe_names = False,
        group = lambda wildcards: wildcards.group,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
        jobname=lambda wildcards: "prune_and_concat_"+wildcards.group
    resources:
        mem_mb = 8000
    threads: 4
    log: "logs/{group}_prune_and_concat.log"
    benchmark: "logs/{group}_prune_and_concat.benchmark"
    conda: "envs/scripts-env.yaml"
    shell:
        """
        python3 scripts/prune_and_concat_alns.py --aligned_dir {input.aligned_dir} --pruned_dir {output.pruned_dir} --pipe_names {params.pipe_names} --genomes_list {input.genomes_list} --group {params.group}
        """

rule make_tree:
    input:
        pruned_cat = "07_AnnotationAndPhylogenies/04_pruned_and_concat_alignments/{group}/CoreGeneAlignment.fasta"
    output:
        treefile = "07_AnnotationAndPhylogenies/05_IQTree/{group}/{group}_Phylogeny.treefile",
        contree = "07_AnnotationAndPhylogenies/05_IQTree/{group}/{group}_Phylogeny.contree",
        iqlog = "07_AnnotationAndPhylogenies/05_IQTree/{group}/{group}_Phylogeny.log"
    params:
        outdir = "07_AnnotationAndPhylogenies/05_IQTree/{group}/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("1-2:10:00") if wildcards.group == "g__Lactobacillus" or wildcards.group == "g__Bifidobacterium" else convertToSec("0-10:10:00"),
        jobname=lambda wildcards: "make_tree_"+wildcards.group
    resources:
        mem_mb = convertToMb("50G")
    threads: 16
    log: "logs/{group}_make_tree.log"
    benchmark: "logs/{group}_make_tree.benchmark"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        if [ -f {output.iqlog} ]; then
            iqtree -s {input.pruned_cat} \
                    -st AA -nt {threads} \
                    -bb 1000 -seed 12345 -m TEST --undo \
                    -pre {params.outdir}{wildcards.group}_Phylogeny
        else
            mkdir -p 07_AnnotationAndPhylogenies/05_IQTree
            mkdir -p {params.outdir}
            iqtree -s {input.pruned_cat} \
                    -st AA -nt {threads} \
                    -bb 1000 -seed 12345 -m TEST \
                    -pre {params.outdir}{wildcards.group}_Phylogeny
        fi
        """

# rule make_contig_tracker:
#     input:
#         faa_files = lambda wildcards: "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa",
#     output:
#         outfile = "06_MAG_binning/contig_tracker_after_prokka/{genome}_contig_tracker.tsv"
#     shell:
#         """
#         # prokka renames contigs
#         # write the names of the contig according to spades and according to the prokka fna file here
#         """

rule concat_all_mags:
    input:
        all_mags = lambda wildcards: expand("07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.fna", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        ffn_files = lambda wildcards: expand("database/ffn_files/{genome}.ffn", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        gff_files = lambda wildcards: expand("database/gff_files/{genome}.gff", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        faa_files = lambda wildcards: expand("database/faa_files/{genome}.faa", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
    output:
        mag_database = "database/MAGs_database",
        bwa_index = multiext("database/MAGs_database", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        samtools_index = "database/MAGs_database.fai",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    conda: "envs/mapping-env.yaml"
    log: "logs/concat_all_mags.log"
    benchmark: "logs/concat_all_mags.benchmark"
    shell:
        """
        cat {input.all_mags} > {output.mag_database}
        bwa index {output.mag_database}
        samtools faidx {output.mag_database}
        """

rule select_cluster_ref_genomes_mag_database:
    input:
        genomes_info = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
        genome_scores = lambda wildcards: rules.drep.output.drep_S
    output:
        mapfile = "database/MAGs_database_ref_info.txt"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/select_cluster_ref_genomes_mag_database.benchmark"
    log: "logs/select_cluster_ref_genomes_mag_database.log"
    run:
        genome_scores = {}
        with open(input.genome_scores, "r") as win_fh:
            for line in win_fh:
                line = line.strip()
                if line.startswith("genome"):
                    continue
                genome_scores[line.split(",")[0].split(".fa")[0]] = float(line.split(",")[1])
        with open(input.genomes_info, "r") as info_fh:
            # make rep genomes set
            group_max_score = {}
            rep_genome_dict = {}
            for line in info_fh:
                if line.startswith("ID"):
                    continue
                genome_id = line.split("\t")[0]
                if "MAG_" not in genome_id:
                    continue
                cluster = line.split("\t")[11]
                # genus = line.split("\t")[13]
                # genus = genus+cluster if genus == "g__" else genus
                if cluster in group_max_score.keys():
                    if int(genome_scores[genome_id]) > group_max_score[cluster]:
                        group_max_score[cluster] = genome_scores[genome_id]
                        rep_genome_dict[cluster] = genome_id
                else:
                    group_max_score[cluster] = genome_scores[genome_id]
                    rep_genome_dict[cluster] = genome_id
        with open(input.genomes_info, "r") as info_fh:
            with open(output.mapfile, "w") as out_fh:
                out_fh.write("genome_id\tcluster\trep_genome\tgroup\n")
                for line in info_fh:
                    if line.startswith("ID"):
                        continue
                    genome_id = line.split("\t")[0]
                    if "MAG_" not in genome_id:
                        continue
                    cluster = line.split("\t")[11]
                    genus = line.split("\t")[13]
                    if genus == "NA":
                        print(f"genus for {genome_id} is {genus}. This group may not make sense.")
                        continue
                    species_id = cluster
                    group_name = genus+cluster if genus == "g__" else genus
                    rep_genome = 1 if genome_id in rep_genome_dict.values() else 0
                    out_fh.write(f"{genome_id}\t{species_id}\t{rep_genome}\t{group_name}\n")

rule copy_mag_database_annotation:
    input:
        ffn_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.ffn",
        gff_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.gff",
        faa_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa",
    output:
        ffn_file = "database/ffn_files/{genome}.ffn",
        gff_file = "database/gff_files/{genome}.gff",
        faa_file = "database/faa_files/{genome}.faa",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "logs/copy_mag_database_annotation_{genome}.log"
    shell:
        """
        rsync -av {input.ffn_file} {output.ffn_file}
        rsync -av {input.gff_file} {output.gff_file}
        rsync -av {input.faa_file} {output.faa_file}
        """

rule prepare_faa_mag_database:
    input:
        info = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
        prokka_faa_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa"
    output:
        faa_file = "database/MAGs_database_Orthofinder/{group}/{genome}.faa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit"
    log: "logs/prepare_faa_mag_database_{genome}_{group}.log"
    shell:
        """
        cp {input.prokka_faa_file} {output.faa_file}
        """

rule run_orthofinder_mag_database:
    input:
        faa_files = lambda wildcards: expand("database/MAGs_database_Orthofinder/{{group}}/{genome}.faa", genome=get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt)[wildcards.group]),
    output:
        orthogroups = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.txt",
        orthogroup_counts = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.GeneCount.tsv",
        orthogroup_stats = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/Results_{group}/Comparative_Genomics_Statistics/Statistics_Overall.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("1-2:10:00") if wildcards.group == "g__Lactobacillus" or wildcards.group == "g__Bifidobacterium" else convertToSec("0-5:10:00"),
        # jobname="run_orthofinder_mag_database_{group}"
    resources:
         mem_mb = convertToMb("20G")
    threads: 8
    log: "logs/{group}_run_orthofinder_mag_database.log"
    benchmark: "logs/{group}_run_orthofinder_mag_database.benchmark"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        rm -rf $(dirname {input.faa_files[0]})/OrthoFinder/Results_{wildcards.group}
        orthofinder -og -t {threads} -n {wildcards.group} -f $(dirname {input.faa_files[0]})  2>&1 | tee -a {log}
        """

# add a file which summarises which MAGs may not be represented at all in the
# half core orthogroups
rule summarise_orthogroups_mag_database:
    input:
        ortho_file = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.txt",
        genomes_list = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
    output:
        summary_orthogroups = "database/MAGs_database_Orthofinder/{group}_Orthogroups_summary.csv"
    params:
        group = lambda wildcards: wildcards.group,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
    log: "logs/{group}_summarise_orthogroups_mag_database.log"
    benchmark: "logs/{group}_summarise_orthogroups_mag_database.benckmark"
    script:
        "scripts/summarise_orthogroups.py"

rule get_single_ortho_mag_database:
    input:
        ortho_file = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.txt",
        genomes_list = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
    output:
        single_ortho = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho.txt",
    params:
        group = lambda wildcards: wildcards.group,
        mailto = "aiswarya.prasad@unil.ch",
        mailtype = "BEGIN,END,FAIL,TIME_LIMIT_80",
        # jobname="{group}_extract_orthologs",
        account = "pengel_spirit",
        runtime_s = convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{group}_get_single_ortho_mag_database.log"
    benchmark: "logs/{group}_get_single_ortho_mag_database.benchmark"
    script:
        "scripts/get_single_ortho_phylo.py"

rule extract_orthologs_mag_database:
    input:
        ortho_single = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho.txt",
        ffn_files = lambda wildcards: expand("database/ffn_files/{genome}.ffn", genome=get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt)[wildcards.group]),
        faa_files = lambda wildcards: expand("database/faa_files/{genome}.faa", genome=get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt)[wildcards.group]),
    output:
        done = touch("database/MAGs_database_Orthofinder/{group}/single_ortholog_sequences.done")
    params:
        # to the prefix, script adds, "genome/genome.xxx"
        ortho_seq_dir = lambda wildcards: "database/MAGs_database_Orthofinder/"+wildcards.group+"/single_ortholog_sequences/",
        ffndir = "database/ffn_files/",
        faadir = "database/faa_files/",
        # faaffndir = "07_AnnotationAndPhylogenies/01_prokka/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname=lambda wildcards: wildcards.group+"_extract_orthologs",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{group}_extract_orthologs_mag_database.log"
    benchmark: "logs/{group}_extract_orthologs_mag_database.benchmark"
    conda: "envs/phylogenies-env.yaml"
    threads: 5
    shell:
        """
        python3 scripts/extract_orthologs_phylo.py --orthofile {input.ortho_single} --outdir {params.ortho_seq_dir} --faadir {params.faadir} --ffndir {params.ffndir}
        """

rule calc_perc_id_mag_database:
    input:
        single_ortho_seq_done = "database/MAGs_database_Orthofinder/{group}/single_ortholog_sequences.done",
        genome_info_path = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
    output:
        perc_id = "database/MAGs_database_Orthofinder/{group}/{group}_perc_id.txt",
        done = touch("database/MAGs_database_Orthofinder/{group}/{group}_perc_id_all.done")
    params:
        ortho_seq_dir = lambda wildcards: "database/MAGs_database_Orthofinder/"+wildcards.group+"/single_ortholog_sequences/",
        pwd_prefix = os.path.join(os.getcwd()),
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:10:00"),
        jobname=lambda wildcards: wildcards.group+"_calc_perc_id",
    resources:
        mem_mb = 8000
    conda: "envs/phylogenies-env.yaml"
    log: "logs/{group}_calc_perc_id_mag_database.log"
    benchmark: "logs/{group}_calc_perc_id_mag_database.benchmark"
    threads: 5
    shell:
        """
        scripts_dir={params.pwd_prefix}/scripts
        genome_db_meta={params.pwd_prefix}/{input.genome_info_path}
        cd {params.pwd_prefix}/{params.ortho_seq_dir}
        outfile={params.pwd_prefix}/{output.perc_id}
        rm -rf $outfile
        for j in $( ls *.faa ); do
            OG=${{j%.\"faa\"}}
            echo \"Processing: $OG\"
            aln_faa=$OG\"_aln.fasta\"
            #Aligning amino-acid sequences
            if [ ! -f $aln_faa ];
            then
                mafft --auto --quiet $j > $aln_faa
            fi
            #Back-translating alignment (codon-aligned nucleotide alignment)
            ffn_file=$OG\".ffn\"
            if [ ! -f $ffn_file ];
            then
                python3 \"${{scripts_dir}}/aln_aa_to_dna.py\" \"$aln_faa\" \"$ffn_file\"
            fi
            #Trimming alignment for gaps
            aln_nuc=$OG\"_aln_nuc.fasta\"
            trim_file=$OG\"_aln_trim.fasta\"
            if [ ! -f $trim_file ];
            then
                python3 \"${{scripts_dir}}/trim_aln.py\"  \"$aln_nuc\" \"$trim_file\"
            fi
            #Calculating inter-SDP alignment stats
            python3 \"${{scripts_dir}}/calc_perc_id_orthologs.py\" --meta \"$genome_db_meta\" --trim_file \"$trim_file\" --outfile \"$outfile\"
        done
        """

rule filter_orthologs_mag_database:
    input:
        single_ortho = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho.txt",
        perc_id = "database/MAGs_database_Orthofinder/{group}/{group}_perc_id.txt",
        done = "database/MAGs_database_Orthofinder/{group}/{group}_perc_id_all.done",
        single_ortho_seq_done = "database/MAGs_database_Orthofinder/{group}/single_ortholog_sequences.done",
    output:
        ortho_filt = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho_filt.txt",
    params:
        ortho_seq_dir = lambda wildcards: "database/MAGs_database_Orthofinder/"+wildcards.group+"/single_ortholog_sequences/",
        mailto = "aiswarya.prasad@unil.ch",
        mailtype = "BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname = lambda wildcards: wildcards.group+"_filter_orthologs_mag_database",
        account = "pengel_spirit",
        runtime_s = convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{group}_filter_orthologs_mag_database.log"
    benchmark: "logs/{group}_filter_orthologs_mag_database.benchmark"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        python3 scripts/filter_orthologs_phylo.py --single_ortho {input.single_ortho} --perc_id {input.perc_id} --extracted_ffndir {params.ortho_seq_dir} --ortho_filt {output.ortho_filt}
        """

rule summarise_orthogroups_filtered_mag_database:
    input:
        ortho_file = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho.txt",
        ortho_file_filt = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho_filt.txt",
        perc_id = "database/MAGs_database_Orthofinder/{group}/{group}_perc_id.txt",
        done = "database/MAGs_database_Orthofinder/{group}/{group}_perc_id_all.done",
        genomes_list = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
        ref_info = "database/MAGs_database_ref_info.txt"
    output:
        summary_orthogroups_filt = "database/MAGs_database_Orthofinder/{group}_Orthogroups_filtered_summary.csv"
    params:
        group = lambda wildcards: wildcards.group,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
    log: "logs/{group}_summarise_orthogroups_filtered_mag_database.log"
    benchmark: "logs/{group}_summarise_orthogroups_filtered_mag_database.benckmark"
    script:
        "scripts/summarise_orthogroups_filtered.py"


rule make_bed_files_mag_database:
    input:
        gff_file = "database/gff_files/{genome}.gff",
        faa_file = "database/faa_files/{genome}.faa",
    output:
        outfile = "database/bed_files/{genome}.bed",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
    log: "logs/{genome}_make_bed_files_mag_database.log"
    benchmark: "logs/{genome}_make_bed_files_mag_database.benckmark"
    # needs biopython
    # convert to shell accorsingly
    #  re-run summarise after filter - consider re-running filter
    # with the additional info
    script:
        "scripts/parse_gff_to_bed.py"
#
#
# edit core_cov to work one bam file at a time
# and from shell within a conda env
# rule core_cov:
#     input:
#         genome_db_meta = "database/"+DBs["microbiome"]+"_metafile.txt",
#         bed_files = expand("database/bed_files/{genome}.bed", genome=get_g_list(DBs["microbiome"])),
#         bam_files = expand("03_MicrobiomeMapping/{sample}_microbiome_mapped.bam", sample=config["SAMPLES"]),
#         ortho_file = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_single_ortho_filt.txt",
#         core_cov = "scripts/core_cov.py",
#         ref_info = "database/MAGs_database_ref_info.txt"
#     output:
#         txt = "04_CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov.txt"
#     params:
#         out_dir = "04_CoreCov_"+ProjectIdentifier,
#         bed_files = "database/bed_files/",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     log: "logs/{phylotype}_core_cov.log"
#     conda: "envs/core-cov-env.yaml"
#     threads: 5
#     shell:
#         """
#         {input.core_cov} -d {input.genome_db_meta} -L {input.bam_files} -g {input.ortho_file} -b {params.bed_files} -o {params.out_dir} -n {wildcards.phylotype};
#         """
#
# rule core_cov_plots:
#     input:
#         txt = "04_CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov.txt",
#         core_covR = "scripts/core_cov.R",
#     output:
#         txt = "04_CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov_coord.txt",
#         pdf = "04_CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov_coord.pdf"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-1:10:00"),
#     resources:
#         mem_mb = 8000
#     conda: "envs/core-cov-env.yaml"
#     log: "logs/{phylotype}_core_cov_plots.log"
#     threads: 1
#     # at some point do:
#     # script:
#     #     "scripts/core_cov.R"
#     shell:
#         """
#         {input.core_covR} {input.txt};
#         """


rule map_to_MAGs:
    input:
        reads1 = rules.trim.output.reads1,
        reads2 = rules.trim.output.reads2,
        reads1_hostfiltered = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2_hostfiltered = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        index = multiext("database/MAGs_database", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        genome_db = "database/MAGs_database"
    output:
        bam_mapped = "09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_mapped.bam",
        bam = temp("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}.bam"),
        sam = temp("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}.sam"),
        sam_mapped = temp("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}mapped.sam"),
        mag_mapping_flagstat = "09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_flagstat.tsv",
        bam_mapped_hostfiltered = temp("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_host-filtered_mapped.bam"),
        bam_hostfiltered = temp("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_host-filtered.bam"),
        sam_hostfiltered = temp("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_host-filtered.sam"),
        sam_mapped_hostfiltered = temp("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_host-filtered_mapped.sam"),
        mag_mapping_hostfiltered_flagstat = "09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_host-filtered_flagstat.tsv",
    params:
        perl_extract_mapped = "scripts/filter_sam_aln_length.pl",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_microbiomedb_mapping",
        account="pengel_spirit",
        runtime_s=convertToSec("0-4:10:00"),
    resources:
        mem_mb = convertToMb("20G")
    threads: 8
    log: "logs/{sample}_map_to_mags.log"
    benchmark: "logs/{sample}_map_to_mags.benchmark"
    conda: "envs/mapping-env.yaml"
    shell:
        """
        bwa mem -t {threads} {input.genome_db} {input.reads1} {input.reads2} > {output.sam}
        samtools view -h -F4 -@ {threads} {output.sam} | samtools view - -F 0x800 -h | perl {params.perl_extract_mapped} - > {output.sam_mapped}
        samtools view -bh {output.sam_mapped} | samtools sort - > {output.bam_mapped}
        samtools view -bh {output.sam} | samtools sort - > {output.bam}
        samtools flagstat -O tsv {output.bam} > {output.mag_mapping_flagstat}

        bwa mem -t {threads} {input.genome_db} {input.reads1_hostfiltered} {input.reads2_hostfiltered} > {output.sam_hostfiltered}
        samtools view -h -F4 -@ {threads} {output.sam_hostfiltered} | samtools view - -F 0x800 -h | perl {params.perl_extract_mapped} - > {output.sam_mapped_hostfiltered}
        samtools view -bh {output.sam_mapped_hostfiltered} | samtools sort - > {output.bam_mapped_hostfiltered}
        samtools view -bh {output.sam_hostfiltered} | samtools sort - > {output.bam_hostfiltered}
        samtools flagstat -O tsv {output.bam_hostfiltered} > {output.mag_mapping_hostfiltered_flagstat}
        """

rule setup_midas_parse_gff:
    input:
        gff_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.gff",
    output:
        outfile = "11_MIDAS/MAGs_database/{genome}/{genome}.genes",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/setup_midas_parse_gff_{genome}.benchmark"
    log: "logs/setup_midas_parse_gff_{genome}.log"
    shell:
        """
        python3 scripts/parse_gff_for_midas.py --outfile {output.outfile} --gff {input.gff_file}
        """

rule setup_midas_for_custom_db_prepare_files:
    input:
        fna_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.fna",
        faa_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa",
        ffn_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.ffn",
    output:
        fna_file = "11_MIDAS/MAGs_database/{genome}/{genome}.fna",
        faa_file = "11_MIDAS/MAGs_database/{genome}/{genome}.faa",
        ffn_file = "11_MIDAS/MAGs_database/{genome}/{genome}.ffn"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "logs/setup_midas_for_custom_db_prepare_files_{genome}.log"
    shell:
        """
        mkdir -p 11_MIDAS/MAGs_database/
        mkdir -p 11_MIDAS/MAGs_database/{wildcards.genome}
        rsync -av {input.fna_file} {output.fna_file}
        rsync -av {input.faa_file} {output.faa_file}
        rsync -av {input.ffn_file} {output.ffn_file}
        """

rule setup_midas_for_custom_db_write_mapfile:
    input:
        genomes_info = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
        genome_scores = lambda wildcards: rules.drep.output.drep_S
    output:
        mapfile = "11_MIDAS/midas_mapfile.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/setup_midas_for_custom_db_write_mapfile.benchmark"
    log: "logs/setup_midas_for_custom_db_write_mapfile.log"
    run:
        genome_scores = {}
        with open(input.genome_scores, "r") as win_fh:
            for line in win_fh:
                line = line.strip()
                if line.startswith("genome"):
                    continue
                genome_scores[line.split(",")[0].split(".fa")[0]] = float(line.split(",")[1])
        with open(input.genomes_info, "r") as info_fh:
            # make rep genomes set
            group_max_score = {}
            rep_genome_dict = {}
            for line in info_fh:
                if line.startswith("ID"):
                    continue
                genome_id = line.split("\t")[0]
                if "MAG_" not in genome_id:
                    continue
                cluster = line.split("\t")[11]
                # genus = line.split("\t")[13]
                # genus = genus+cluster if genus == "g__" else genus
                if cluster in group_max_score.keys():
                    if int(genome_scores[genome_id]) > group_max_score[cluster]:
                        group_max_score[cluster] = genome_scores[genome_id]
                        rep_genome_dict[cluster] = genome_id
                else:
                    group_max_score[cluster] = genome_scores[genome_id]
                    rep_genome_dict[cluster] = genome_id
        with open(input.genomes_info, "r") as info_fh:
            with open(output.mapfile, "w") as out_fh:
                out_fh.write("genome_id\tspecies_id\trep_genome\n")
                for line in info_fh:
                    if line.startswith("ID"):
                        continue
                    genome_id = line.split("\t")[0]
                    if "MAG_" not in genome_id:
                        continue
                    cluster = line.split("\t")[11]
                    # genus = line.split("\t")[13]
                    # if genus == "NA":
                    #     print(f"genus for {genome_id} is {genus}. This group may not make sense.")
                    #     continue
                    species_id = cluster
                    # species_id = genus+cluster if genus == "g__" else genus
                    rep_genome = 1 if genome_id in rep_genome_dict.values() else 0
                    out_fh.write(f"{genome_id}\t{species_id}\t{rep_genome}\n")

rule create_midas_for_custom_db:
    input:
        fna_files = lambda wildcards: expand("11_MIDAS/MAGs_database/{genome}/{genome}.fna", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        faa_files = lambda wildcards: expand("11_MIDAS/MAGs_database/{genome}/{genome}.faa", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        ffn_files = lambda wildcards: expand("11_MIDAS/MAGs_database/{genome}/{genome}.ffn", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        genes_file = lambda wildcards: expand("11_MIDAS/MAGs_database/{genome}/{genome}.genes", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        mapfile = "11_MIDAS/midas_mapfile.tsv"
    output:
        done = touch("11_MIDAS/Midas_DB.done")
    params:
        outdir = "11_MIDAS/MAGs_database/Midas_DB",
        # consider putting these paths in config
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        midas_db = "11_MIDAS/MAGs_database/Midas_DB",
        # midas_db_default = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/midas_db_v1.2/",
        indir = "11_MIDAS/MAGs_database/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:10:00"),
    benchmark: "logs/create_midas_for_custom_db.benchmark"
    log: "logs/create_midas_for_custom_db.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        build_midas_db.py {params.indir} {input.mapfile} {params.outdir} --threads {threads} --resume
        if [ ! -f {params.midas_db}/marker_genes/phyeco.fa.bwt ];
        then
            hs-blastn index {params.midas_db}/marker_genes/phyeco.fa
        fi
        """

rule midas_species_abundance:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        done = "11_MIDAS/Midas_DB.done"
    output:
        output = "11_MIDAS/MidasProfiling/{sample}/species/species_profile.txt",
    params:
        outdir = "11_MIDAS/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_species_abundance_{sample}.benchmark"
    log: "logs/midas_species_abundance_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        run_midas.py species {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_species_abundance_merge:
    input:
        per_sample_output = expand("11_MIDAS/MidasProfiling/{sample}/species/species_profile.txt", sample=SAMPLES)
    output:
        relative_abundance = "11_MIDAS/MidasProfiling/relative_abundance.txt",
        count_reads = "11_MIDAS/MidasProfiling/count_reads.txt",
        coverage = "11_MIDAS/MidasProfiling/coverage.txt",
        species_prevalence = "11_MIDAS/MidasProfiling/species_prevalence.txt",
    params:
        outdir = "11_MIDAS/MidasProfiling/",
        midas_db = "11_MIDAS/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_species_abundance_merged.benchmark"
    log: "logs/midas_species_abundance_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py species {params.outdir} -i {params.outdir} -t dir
        """

rule midas_gene_content:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        output = "11_MIDAS/MidasProfiling/{sample}/species/species_profile.txt",
        done = "11_MIDAS/Midas_DB.done",
    output:
        output = "11_MIDAS/MidasProfiling/{sample}/genes/summary.txt",
    params:
        outdir = "11_MIDAS/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        species_cov = 10, # default is 3 use later if needed
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_gene_content_{sample}.benchmark"
    log: "logs/midas_gene_content_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        run_midas.py genes {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_gene_content_merge:
    input:
        per_sample_output = expand("11_MIDAS/MidasProfiling/{sample}/genes/summary.txt", sample=SAMPLES)
    output:
        # genes_coverage = "11_MIDAS/MidasProfiling/{cluster}/genes_coverage.txt",
        # genes_presence = "11_MIDAS/MidasProfiling/{cluster}/genes_presence_absence.txt",
        # genes_summary = "11_MIDAS/MidasProfiling/{cluster}/genes_summary.txt",
        # genes_copy_number = "11_MIDAS/MidasProfiling/{cluster}/genes_copy_number.txt"
        outdone = touch("11_MIDAS/midas_genes_summary.done")
    params:
        outdir = "11_MIDAS/MidasProfiling/",
        midas_db = "11_MIDAS/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_gene_content_merged.benchmark"
    log: "logs/midas_gene_content_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py genes {params.outdir} -i {params.outdir} -t dir
        for readme in 11_MIDAS/MidasProfiling/*/readme.txt;
        do
            cp ${{readme%%readme.txt}}/readme.txt ${{readme%%readme.txt}}/genes_readme.txt
        done
        """

rule midas_profile_snps:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        output = "11_MIDAS/MidasProfiling/{sample}/species/species_profile.txt",
        done = "11_MIDAS/Midas_DB.done"
    output:
        summary = "11_MIDAS/MidasProfiling/{sample}/snps/summary.txt",
    params:
        outdir = "11_MIDAS/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        species_cov = 1, # default is 3 use later if needed
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_profile_snps_{sample}.benchmark"
    log: "logs/midas_profile_snps_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        run_midas.py snps {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_profile_snps_merge:
    input:
        per_sample_output = expand("11_MIDAS/MidasProfiling/{sample}/snps/summary.txt", sample=SAMPLES)
    output:
        # snps_freq = "11_MIDAS/MidasProfiling/{cluster}/snps_freq.txt",
        # snps_depth = "11_MIDAS/MidasProfiling/{cluster}/snps_depth.txt",
        # snps_summary = "11_MIDAS/MidasProfiling/{cluster}/snps_summary.txt",
        # snps_info = "11_MIDAS/MidasProfiling/{cluster}/snps_info.txt"
        outdone = touch("11_MIDAS/snps_summary.done")
    params:
        outdir = "11_MIDAS/MidasProfiling/",
        midas_db = "11_MIDAS/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_profile_snps_merged.benchmark"
    log: "logs/midas_profile_snps_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py snps {params.outdir} -i {params.outdir} -t dir
        for readme in 11_MIDAS/MidasProfiling/*/readme.txt;
        do
            cp ${{readme%%readme.txt}}/readme.txt ${{readme%%readme.txt}}/snps_readme.txt
        done
        """

rule setup_midas_genus_parse_gff:
    input:
        gff_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.gff",
    output:
        outfile = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.genes",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/setup_midas_genus_parse_gff_{genome}.benchmark"
    log: "logs/setup_midas_genus_parse_gff_{genome}.log"
    shell:
        """
        python3 scripts/parse_gff_for_midas.py --outfile {output.outfile} --gff {input.gff_file}
        """

rule setup_midas_genus_for_custom_db_prepare_files:
    input:
        fna_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.fna",
        faa_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa",
        ffn_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.ffn",
    output:
        fna_file = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.fna",
        faa_file = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.faa",
        ffn_file = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.ffn"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "logs/setup_midas_genus_for_custom_db_prepare_files_{genome}.log"
    shell:
        """
        mkdir -p 11_MIDAS_genus_level/MAGs_database/
        mkdir -p 11_MIDAS_genus_level/MAGs_database/{wildcards.genome}
        rsync -av {input.fna_file} {output.fna_file}
        rsync -av {input.faa_file} {output.faa_file}
        rsync -av {input.ffn_file} {output.ffn_file}
        """

rule setup_midas_genus_for_custom_db_write_mapfile:
    input:
        genomes_info = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
        genome_scores = lambda wildcards: rules.drep.output.drep_S,
    output:
        mapfile = "11_MIDAS_genus_level/midas_mapfile.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/setup_midas_genus_for_custom_db_write_mapfile.benchmark"
    log: "logs/setup_midas_genus_for_custom_db_write_mapfile.log"
    run:
        genome_scores = {}
        with open(input.genome_scores, "r") as win_fh:
            for line in win_fh:
                line = line.strip()
                if line.startswith("genome"):
                    continue
                genome_scores[line.split(",")[0].split(".fa")[0]] = float(line.split(",")[1])
        with open(input.genomes_info, "r") as info_fh:
            # make rep genomes set
            group_max_score = {}
            rep_genome_dict = {}
            for line in info_fh:
                if line.startswith("ID"):
                    continue
                genome_id = line.split("\t")[0]
                if "MAG_" not in genome_id:
                    continue
                cluster = line.split("\t")[11]
                genus = line.split("\t")[13]
                genus = genus+cluster if genus == "g__" else genus
                if genus in group_max_score.keys():
                    if int(genome_scores[genome_id]) > group_max_score[genus]:
                        group_max_score[genus] = genome_scores[genome_id]
                        rep_genome_dict[genus] = genome_id
                else:
                    group_max_score[genus] = genome_scores[genome_id]
                    rep_genome_dict[genus] = genome_id
        with open(input.genomes_info, "r") as info_fh:
            with open(output.mapfile, "w") as out_fh:
                out_fh.write("genome_id\tspecies_id\trep_genome\n")
                for line in info_fh:
                    if line.startswith("ID"):
                        continue
                    genome_id = line.split("\t")[0]
                    if "MAG_" not in genome_id:
                        continue
                    cluster = line.split("\t")[11]
                    genus = line.split("\t")[13]
                    # if genus == "NA":
                    #     print(f"genus for {genome_id} is {genus}.")
                    #     continue
                    species_id = genus+cluster if genus == "g__" else genus
                    rep_genome = 1 if genome_id in rep_genome_dict.values() else 0
                    out_fh.write(f"{genome_id}\t{species_id}\t{rep_genome}\n")

rule create_midas_genus_for_custom_db:
    input:
        fna_files = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.fna", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        faa_files = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.faa", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        ffn_files = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.ffn", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        genes_file = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.genes", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        mapfile = "11_MIDAS_genus_level/midas_mapfile.tsv"
    output:
        done = touch("11_MIDAS_genus_level/Midas_DB.done")
    params:
        outdir = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        # consider putting these paths in config
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        # midas_db_default = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/midas_db_v1.2/",
        indir = "11_MIDAS_genus_level/MAGs_database/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:10:00"),
    benchmark: "logs/create_midas_genus_for_custom_db.benchmark"
    log: "logs/create_midas_genus_for_custom_db.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        build_midas_db.py {params.indir} {input.mapfile} {params.outdir} --threads {threads} --resume
        if [ ! -f {params.midas_db}/marker_genes/phyeco.fa.bwt ];
        then
            hs-blastn index {params.midas_db}/marker_genes/phyeco.fa
        fi
        """

rule midas_genus_species_abundance:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        done = "11_MIDAS_genus_level/Midas_DB.done"
    output:
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_species_abundance_{sample}.benchmark"
    log: "logs/midas_genus_species_abundance_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        run_midas.py species {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_genus_species_abundance_merge:
    input:
        per_sample_output = expand("11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt", sample=SAMPLES)
    output:
        relative_abundance = "11_MIDAS_genus_level/MidasProfiling/relative_abundance.txt",
        count_reads = "11_MIDAS_genus_level/MidasProfiling/count_reads.txt",
        coverage = "11_MIDAS_genus_level/MidasProfiling/coverage.txt",
        species_prevalence = "11_MIDAS_genus_level/MidasProfiling/species_prevalence.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_species_abundance_merged.benchmark"
    log: "logs/midas_genus_species_abundance_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py species {params.outdir} -i {params.outdir} -t dir
        """

rule midas_genus_gene_content:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt",
        done = "11_MIDAS_genus_level/Midas_DB.done",
    output:
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/genes/summary.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        species_cov = 10, # default is 3 use later if needed
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_gene_content_{sample}.benchmark"
    log: "logs/midas_genus_gene_content_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        run_midas.py genes {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_genus_gene_content_merge:
    input:
        per_sample_output = expand("11_MIDAS_genus_level/MidasProfiling/{sample}/genes/summary.txt", sample=SAMPLES)
    output:
        # genes_coverage = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_coverage.txt",
        # genes_presence = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_presence_absence.txt",
        # genes_summary = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_summary.txt",
        # genes_copy_number = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_copy_number.txt"
        outdone = touch("11_MIDAS_genus_level/midas_genes_summary.done")
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_gene_content_merged.benchmark"
    log: "logs/midas_genus_gene_content_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py genes {params.outdir} -i {params.outdir} -t dir
        for readme in 11_MIDAS_genus_level/MidasProfiling/*/readme.txt;
        do
            cp ${{readme%%readme.txt}}/readme.txt ${{readme%%readme.txt}}/genes_readme.txt
        done
        """

rule midas_genus_profile_snps:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt",
        done = "11_MIDAS_genus_level/Midas_DB.done"
    output:
        summary = "11_MIDAS_genus_level/MidasProfiling/{sample}/snps/summary.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        species_cov = 1, # default is 3 use later if needed
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_profile_snps_{sample}.benchmark"
    log: "logs/midas_genus_profile_snps_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        run_midas.py snps {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_genus_profile_snps_merge:
    input:
        per_sample_output = expand("11_MIDAS_genus_level/MidasProfiling/{sample}/snps/summary.txt", sample=SAMPLES)
    output:
        # snps_freq = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_freq.txt",
        # snps_depth = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_depth.txt",
        # snps_summary = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_summary.txt",
        # snps_info = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_info.txt"
        outdone = touch("11_MIDAS_genus_level/snps_summary.done")
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_profile_snps_merged.benchmark"
    log: "logs/midas_genus_profile_snps_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py snps {params.outdir} -i {params.outdir} -t dir
        for readme in 11_MIDAS_genus_level/MidasProfiling/*/readme.txt;
        do
            cp ${{readme%%readme.txt}}/readme.txt ${{readme%%readme.txt}}/snps_readme.txt
        done
        """

rule midas_genus_snp_diversity_intra:
    input:
        done = rules.midas_profile_snps_merge.output.outdone
    output:
        outfile = "11_MIDAS_genus_level/snp_diversity_intra/snp_diversity_{cluster}.info"
    params:
        indir = "11_MIDAS_genus_level/MidasProfiling/{cluster}/",
        maf = 0.01, #minimum minor allele frequency to be considered polymorphic
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_snp_diversity_intra_{cluster}.benchmark"
    log: "logs/midas_snp_diversity_intra_{cluster}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        snp_diversity.py {params.indir} --genomic_type genome-wide --sample_type per-sample --out {output.outfile} --snp_maf {params.maf}
        """

rule midas_genus_snp_diversity_inter:
    input:
        done = rules.midas_profile_snps_merge.output.outdone
    output:
        outfile = "11_MIDAS_genus_level/snp_diversity_inter/snp_diversity_{cluster}.info"
    params:
        indir = "11_MIDAS_genus_level/MidasProfiling/{cluster}/",
        maf = 0.01, #minimum minor allele frequency to be considered polymorphic
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_snp_diversity_inter_{cluster}.benchmark"
    log: "logs/midas_snp_diversity_inter_{cluster}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        snp_diversity.py {params.indir} --genomic_type genome-wide --sample_type pooled-samples --out {output.outfile} --snp_maf {params.maf}
        """

###############################
###############################
###############################
# Make report and backup
###############################
###############################
###############################


rule compile_report:
    input:
        rmd = PROJECT_IDENTIFIER+"_Report.Rmd",
        isolates = "config/IsolateGenomeInfo.csv",
        coverage_host = expand("02_HostMapping/{sample}_coverage.tsv", sample=SAMPLES),
        coverage_host_hist = expand("02_HostMapping/{sample}_coverage_histogram.txt", sample=SAMPLES),
        flagstat_02 = expand("02_HostMapping/{sample}_flagstat.tsv", sample=SAMPLES),
        flagstat_03 = expand("03_MicrobiomeMapping/{sample}_flagstat.tsv", sample=SAMPLES),
        flagstat_04 = expand("04_MicrobiomeMappingDirect/{sample}_flagstat.tsv", sample=SAMPLES),
        qc_raw = expand("fastqc/raw/{sample}_{read}_fastqc.html", sample=SAMPLES, read=config["READS"]),
        qc_trimmed = expand("fastqc/trim/{sample}_{read}_trim_fastqc.html", sample=SAMPLES, read=config["READS"]),
        motus_merged = "08_motus_profile/samples_merged.motus",
        summary_assembly = "05_Assembly/host_unmapped/check_assembly/Assembly_mapping_summary.csv",
        contig_fates = expand("06_MAG_binning/contig_fates/{sample}_contig_fates.csv", sample=SAMPLES),
        sample_contig_coverages = expand("06_MAG_binning/contig_fates/backmapping_coverages/{sample}_contig_coverages.csv", sample=SAMPLES),
        out_all = "06_MAG_binning/all_GenomeInfo_auto.tsv",
        out_tree = "06_MAG_binning/ForTree_GenomeInfo_auto.tsv",
        checkpoint = lambda wildcards: checkpoints.make_phylo_table.get().output.out_tree,
        trees = lambda wildcards: ["07_AnnotationAndPhylogenies/05_IQTree/"+group+"/"+group+"_Phylogeny.contree" for group in [x for x in GROUPS if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_tree) > 4]],
        single_ortho_MAGs = expand("07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/{group}_single_ortho_MAGs.txt", group=GROUPS),
        ortho_summary = expand("07_AnnotationAndPhylogenies/02_orthofinder/{group}_Orthogroups_summary.csv", group=GROUPS),
        midas_species_prevalence = "11_MIDAS/MidasProfiling/species_prevalence.txt",
        midas_genes_summary = "11_MIDAS/midas_genes_summary.done",
        snps_summary = "11_MIDAS/snps_summary.done",
        midas_genus_species_prevalence = "11_MIDAS_genus_level/MidasProfiling/species_prevalence.txt",
        midas_genus_genes_summary = "11_MIDAS_genus_level/midas_genes_summary.done",
        snps_summary_genus = "11_MIDAS_genus_level/snps_summary.done",
        # snp_intra = lambda wildcards: expand("11_MIDAS/snp_diversity_intra/snp_diversity_{cluster}.info", cluster = get_cluster_dict(checkpoints.make_phylo_table.get().output.out_all).keys()),
        # snp_inter = lambda wildcards: expand("11_MIDAS/snp_diversity_inter/snp_diversity_{cluster}.info", cluster = get_cluster_dict(checkpoints.make_phylo_table.get().output.out_all).keys())
        # mag_mapping_flagstat = expand("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_flagstat.tsv", sample=SAMPLES),
        # mag_mapping_hostfiltered_flagstat = expand("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_host-filtered_flagstat.tsv", sample=SAMPLES),
        summarise_db_ortho = lambda wildcards: ["database/MAGs_database_Orthofinder/"+group+"_Orthogroups_summary.csv" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
        summarise_db_ortho_filt = lambda wildcards: ["database/MAGs_database_Orthofinder/"+group+"_Orthogroups_filtered_summary.csv" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
        bam_mapped_mag_database = expand("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_mapped.bam", sample=SAMPLES),
        single_ortho_filt_mag_database = lambda wildcards: ["database/MAGs_database_Orthofinder/"+group+"/OrthoFinder/"+group+"_single_ortho_filt.txt" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
    output:
        html = PROJECT_IDENTIFIER+"_Report.html",
        # fig_00a = "Figures/00a-Total_reads_trimming.pdf",
        # fig_00b = "Figures/00b-Total_reads_species_after_trimming.pdf",
        # fig_00c = "Figures/00c-Mapped_Unmapped_reads_prop.pdf",
        # fig_00d = "Figures/00d-Mapping_non_specific.pdf",
        # fig_04b = "Figures/04b-motus_filt.pdf",
        # fig_04a = "Figures/04a-motus.pdf",
        # fig_03a = "Figures/03a-NumberOfFilteredcontigs.pdf",
        # fig_03b = "Figures/03b-Assembly_summary.pdf",
        # fig_05c_am = "Figures/05c-contig_length_histogram_am.pdf",
        # fig_05d_am = "Figures/05d-contig_length_histogram_am_passed.pdf",
        # fig_05c_ac = "Figures/05c-contig_length_histogram_ac.pdf",
        # fig_05d_ac = "Figures/05d-contig_length_histogram_ac_passed.pdf",
        # fig_05c_ad = "Figures/05c-contig_length_histogram_ad.pdf",
        # fig_05d_ad = "Figures/05d-contig_length_histogram_ad_passed.pdf",
        # fig_05c_af = "Figures/05c-contig_length_histogram_af.pdf",
        # fig_05d_af = "Figures/05d-contig_length_histogram_af_passed.pdf",
        # fig_05g_am = "Figures/05g-length_vs_depth_contigs_all_am.pdf",
        # fig_05h_am = "Figures/05h-length_vs_depth_contigs_binned_am.pdf",
        # fig_05g_ac = "Figures/05g-length_vs_depth_contigs_all_ac.pdf",
        # fig_05h_ac = "Figures/05h-length_vs_depth_contigs_binned_ac.pdf",
        # fig_05g_ad = "Figures/05g-length_vs_depth_contigs_all_ad.pdf",
        # fig_05h_ad = "Figures/05h-length_vs_depth_contigs_binned_ad.pdf",
        # fig_05g_af = "Figures/05g-length_vs_depth_contigs_all_af.pdf",
        # fig_05h_af = "Figures/05h-length_vs_depth_contigs_binned_af.pdf",
        # fig_05a = "Figures/05a-data-contigs_passed_failed.pdf",
        # fig_05b = "Figures/05b-contigs_binned_unbinned.pdf",
        # fig_05e_1 = "Figures/05e-contigs_binned_unbinned_by_genus.pdf",
        # fig_05e_2 = "Figures/05e-contigs_binned_unbinned_by_genus_scaled.pdf",
        # fig_05f_1 = "Figures/05f-contigs_by_genus_binned.pdf",
        # fig_05f_2 = "Figures/05f-contigs_by_genus_binned_scaled.pdf",
        # fig_07a_1 = "Figures/07a-prev_vs_coverage_all_MAGs_genus_by_host.pdf",
        # fig_07a_2 = "Figures/07a-prev_vs_coverage_filtered_MAGs_genus_by_host.pdf",
        # fig_07a_3 = "Figures/07a-prev_vs_coverage_all_MAGs_genus_by_host_no_size.pdf",
        # fig_07a_4 = "Figures/07a-prev_vs_coverage_filtered_MAGs_genus_by_host_no_size.pdf",
        # fig_06a = "Figures/06a-QC_MAG_histogram.pdf",
        # fig_06b = "Figures/06b-QC_MAG_per_sample.pdf",
        # fig_07a_5 = "Figures/07a-prev_vs_coverage_MAGs_genus.pdf",
        # fig_07a_6 = "Figures/07a-prev_overall_vs_coverage_MAGs_genus.pdf",
        # fig_07_1 = "Figures/07-QC_per_Genus_per_host_all.pdf",
        # fig_07_2 = "Figures/07-QC_per_Genus_per_host_passed.pdf",
        # fig_08a = "Figures/08a-magOTUs_per_sample.pdf",
        # fig_08b = "Figures/08b-magOTU_per_sample_genus.pdf",
        # fig_08c_1 = "Figures/08c-magOTU_by_host_genus.pdf",
        # fig_08c_2 = "Figures/08c-magOTU_by_host_MAGs_coverage.pdf",
        # fig_08c_3 = "Figures/08c-magOTU_by_host_MAGs_prevalence.pdf",
        # fig_08d_1 = "Figures/08d-magOTU_shared_per_sample_genus.pdf",
        # fig_08d_2 = "Figures/08d-magOTUs_shared_per_sample_prev_abund.pdf",
        # fig_08e = "Figures/08e-magOTUs_by_host_completeness.pdf",
        # fig_08f_1 = "Figures/08f-mags_prev_vs_abd_by_genus_completeness.pdf",
        # fig_08f_2 = "Figures/08f-mags_prev_vs_abd_by_genus_host.pdf",
        # fig_09a_1 = "Figures/09a-magOTU_number_per_sample_all.pdf",
        # fig_09a_2 = "Figures/09a-magOTU_number_per_sample_passed.pdf",
        # fig_09c = "Figures/09c-magOTUs_accumulation_curve.pdf",
        # fig_09d_1 = "Figures/09d-magOTUs_pcoa.pdf",
        # fig_09d_2 = "Figures/09d-magOTUs_pcoa_location.pdf",
        # fig_09d_3 = "Figures/09d-magOTUs_pcoa_colony.pdf",
        # fig_09d_4 = "Figures/09d-magOTUs_pcoa_country.pdf",
        # fig_09d_5 = "Figures/09d-magOTUs_pcoa_run_id.pdf",
        # fig_10a = "Figures/10a-Number_of_genomes_per_group.pdf",
        # fig_10b = "Figures/10b-Number_of_orthogroups_per_group.pdf",
        # fig_10c = "Figures/10c-Orthogroups_status_per_group.pdf",
        # fig_11a = "Figures/11a-mag_drep_scores_histogram.pdf",
        # fig_11b = "Figures/11b-cluster_size_histogram.pdf",
    conda: "envs/rmd-env.yaml"
    threads: 2
    params:
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = convertToMb("15G")
    log: "logs/compile_report.log"
    benchmark: "logs/compile_report.benchmark"
    shell:
        """
        rm -rf PROJECT_IDENTIFIER_Report_cache
        R -e \"rmarkdown::render('{input.rmd}')\"
        """

rule backup:
    input:
        html = PROJECT_IDENTIFIER+"_Report.html",
        rmd = PROJECT_IDENTIFIER+"_Report.Rmd",
        isolates = "config/IsolateGenomeInfo.csv",
        coverage_host = expand("02_HostMapping/{sample}_coverage.tsv", sample=SAMPLES),
        coverage_host_hist = expand("02_HostMapping/{sample}_coverage_histogram.txt", sample=SAMPLES),
        flagstat_02 = expand("02_HostMapping/{sample}_flagstat.tsv", sample=SAMPLES),
        flagstat_03 = expand("03_MicrobiomeMapping/{sample}_flagstat.tsv", sample=SAMPLES),
        flagstat_04 = expand("04_MicrobiomeMappingDirect/{sample}_flagstat.tsv", sample=SAMPLES),
        qc_raw = expand("fastqc/raw/{sample}_{read}_fastqc.html", sample=SAMPLES, read=config["READS"]),
        qc_trimmed = expand("fastqc/trim/{sample}_{read}_trim_fastqc.html", sample=SAMPLES, read=config["READS"]),
        motus_merged = "08_motus_profile/samples_merged.motus",
        summary_assembly = "05_Assembly/host_unmapped/check_assembly/Assembly_mapping_summary.csv",
        contig_fates = expand("06_MAG_binning/contig_fates/{sample}_contig_fates.csv", sample=SAMPLES),
        sample_contig_coverages = expand("06_MAG_binning/contig_fates/backmapping_coverages/{sample}_contig_coverages.csv", sample=SAMPLES),
        out_all = "06_MAG_binning/all_GenomeInfo_auto.tsv",
        out_tree = "06_MAG_binning/ForTree_GenomeInfo_auto.tsv",
        assembly_mapped = expand("05_Assembly/MapToAssembly/{sample}.bam", sample=SAMPLES),
        checkpoint_tree = lambda wildcards: checkpoints.make_phylo_table.get().output.out_tree,
        checkpoint_all = lambda wildcards: checkpoints.make_phylo_table.get().output.out_all,
        trees = lambda wildcards: ["07_AnnotationAndPhylogenies/05_IQTree/"+group+"/"+group+"_Phylogeny.contree" for group in [x for x in GROUPS if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_tree) > 4]],
        single_ortho_MAGs = expand("07_AnnotationAndPhylogenies/02_orthofinder/{group}/OrthoFinder/{group}_single_ortho_MAGs.txt", group=GROUPS),
        ortho_summary = expand("07_AnnotationAndPhylogenies/02_orthofinder/{group}_Orthogroups_summary.csv", group=GROUPS),
        midas_species_prevalence = "11_MIDAS/MidasProfiling/species_prevalence.txt",
        midas_genes_summary = "11_MIDAS/midas_genes_summary.done",
        snps_summary = "11_MIDAS/snps_summary.done",
        midas_genus_species_prevalence = "11_MIDAS_genus_level/MidasProfiling/species_prevalence.txt",
        midas_genus_genes_summary = "11_MIDAS_genus_level/midas_genes_summary.done",
        snps_summary_genus = "11_MIDAS_genus_level/snps_summary.done",
        # snp_intra = lambda wildcards: expand("11_MIDAS/snp_diversity_intra/snp_diversity_{cluster}.info", cluster = get_cluster_dict(checkpoints.make_phylo_table.get().output.out_all).keys()),
        # snp_inter = lambda wildcards: expand("11_MIDAS/snp_diversity_inter/snp_diversity_{cluster}.info", cluster = get_cluster_dict(checkpoints.make_phylo_table.get().output.out_all).keys())
        # mOTUlize = "06_MAG_binning/mOTUlizer/mOTUlizer_output.tsv",
        # mag_mapping_flagstat = expand("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_flagstat.tsv", sample=SAMPLES),
        # mag_mapping_hostfiltered_flagstat = expand("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_host-filtered_flagstat.tsv", sample=SAMPLES),
        summarise_db_ortho = lambda wildcards: ["database/MAGs_database_Orthofinder/"+group+"_Orthogroups_summary.csv" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
        summarise_db_ortho_filt = lambda wildcards: ["database/MAGs_database_Orthofinder/"+group+"_Orthogroups_filtered_summary.csv" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
        bam_mapped_mag_database = expand("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_mapped.bam", sample=SAMPLES),
        single_ortho_filt_mag_database = lambda wildcards: ["database/MAGs_database_Orthofinder/"+group+"/OrthoFinder/"+group+"_single_ortho_filt.txt" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
    output:
        outfile = touch("logs/backup.done")
    threads: 2
    log: "logs/backup.log"
    benchmark: "logs/backup.benchmark"
    params:
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = convertToMb("4G")
    run:
        if LOCAL_BACKUP:
            shell("echo \' ensure that "+BACKUP_PATH+" exists \'")
        shell("scripts/backup.sh "+PROJECT_PATH+" "+os.path.join(BACKUP_PATH, PROJECT_IDENTIFIER)+" logs/backup.log")


# ###############################
# ###############################
# ###############################
# # SNV calling
# ###############################
# ###############################
# ###############################
#
# rule mapping_red_db:
#     input:
#         db_red = os.path.join(os.getcwd(), "database/"+DBs["microbiome_red"]),
#         reads1 = ancient(rules.trim.output.reads1),
#         reads2 = ancient(rules.trim.output.reads2),
#         index =  multiext("database/"+DBs["microbiome_red"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
#     output:
#         bam_red = temp("07_SNVProfiling/Mapping_red/{sample}_mapped_red.bam"),
#         bam_red1 = temp("07_SNVProfiling/Mapping_red/{sample}_mapped_red1.bam"),
#         sam_red1 = temp("07_SNVProfiling/Mapping_red/{sample}_mapped_red1.sam"),
#         sam_red2 = temp("07_SNVProfiling/Mapping_red/{sample}_mapped_red2.sam"),
#     params:
#         perl_extract_mapped = "scripts/filter_sam_aln_length.pl",
#         picard = "picard",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-6:10:00"),
#     resources:
#         mem_mb = 8000
#     threads:
#         20
#     log: "logs/{sample}_mapping_red_db.log"
#     conda: "envs/mapping-env.yaml"
#     shell:
#         """
#         bwa mem -t {threads} {input.db_red} {input.reads1} {input.reads2} > {output.sam_red1}
#         samtools view -h -F4 -@ {threads} {output.sam_red1} | samtools view - -F 0x800 -h | grep -E "NM:i:[0-4]|^\@" | perl {params.perl_extract_mapped} - > {output.sam_red2}
#         samtools view -bh {output.sam_red2} | samtools sort - > {output.bam_red1}
#         {params.picard} -Xmx8g AddOrReplaceReadGroups I={output.bam_red1} O={output.bam_red} RGID={wildcards.sample} RGLB=lib1 RGPL=illumina RGPU=none RGSM={wildcards.sample}
#         samtools index {output.bam_red}
#         """
#
# rule de_duplicate:
#     input:
#         bam_red = rules.mapping_red_db.output.bam_red,
#     output:
#         bam_red = temp("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam"),
#         txt = "07_SNVProfiling/Mapping_red/{sample}_dedup_metrics.txt"
#     params:
#         db_red = "database/"+DBs["microbiome_red"],
#         perl_extract_mapped = "scripts/filter_sam_aln_length.pl",
#         picard = "picard",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-6:10:00"),
#     resources:
#         mem_mb = 8000
#     threads:
#         20
#     conda: "envs/mapping-env.yaml"
#     log: "logs/{sample}_de_duplicate.log"
#     shell:
#         """
#         picard -Xmx8g MarkDuplicates INPUT={input.bam_red} OUTPUT={output.bam_red} METRICS_FILE={output.txt};
#         """
#
# rule de_duplicate_index:
#     input:
#         bam_red = rules.de_duplicate.output.bam_red,
#     output:
#         bam_red = "07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam.bai",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-6:10:00"),
#     resources:
#         mem_mb = 8000
#     threads:
#         20
#     conda: "envs/mapping-env.yaml"
#     log: "logs/{sample}_de_duplicate_index.log"
#     shell:
#         """
#         samtools index {input.bam_red}
#         """
#
# rule subset_ortho_and_db:
#     input:
#         db = "database/"+DBs["microbiome"],
#         metafile = "database/"+DBs["microbiome"]+"_metafile.txt",
#         orthofiles = expand("database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_single_ortho_filt.txt", phylotype=Phylos),
#     output:
#         db = "database/"+DBs["microbiome_red"],
#         orthofiles = expand("database/OrthoFiles_"+DBs["microbiome_red"]+"/{phylotype}_single_ortho_filt.txt", phylotype=Phylos),
#         orthodir = directory("database/OrthoFiles_"+DBs["microbiome_red"]+"/")
#     log: "logs/subset_ortho_and_meta.log"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 4000
#     script:
#         "scripts/subset_metagenome_db.py"
#
# rule core_cov_red:
#     input:
#         ortho_file = "database/OrthoFiles_"+DBs["microbiome_red"]+"/{phylotype}_single_ortho_filt.txt",
#         bam_files = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam", sample=config["SAMPLES"]),
#         bed_files = expand("database/bed_files/{genome}.bed", genome=get_g_list(DBs["microbiome_red"])),
#         genome_db_meta = "database/"+DBs["microbiome"]+"_metafile.txt",
#         core_cov = "scripts/core_cov.py",
#     output:
#         corecovred_txt = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov.txt"
#     params:
#         bed_files = "database/bed_files",
#         out_dir = "07_SNVProfiling/CoreCov_"+ProjectIdentifier,
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-7:10:00")
#     resources:
#         mem_mb = 8000
#     log: "logs/{phylotype}_core_cov.log"
#     conda: "envs/core-cov-env.yaml"
#     threads: 5
#     shell:
#         """
#         {input.core_cov} -d {input.genome_db_meta} -L {input.bam_files} -g {input.ortho_file} -b {params.bed_files} -o {params.out_dir} -n {wildcards.phylotype};
#         """
#
# rule corecov_split_by_sdp:
#     input:
#         coords = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov_coord.txt", phylotype=Phylos),
#         corecovs = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov.txt", phylotype=Phylos),
#     output:
#         sdp_coords = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_corecov_coord.txt", sdp=get_g_sdp_dict(DBs["microbiome"]).keys()),
#         sdp_corecovs = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_corecov.txt", sdp=get_g_sdp_dict(DBs["microbiome"]).keys()),
#     log: "logs/corecov_split_by_sdp.log"
#     threads: 1
#     script:
#         "scripts/corecov_split_by_sdp.py"
#
# rule core_cov_red_plots:
#     input:
#         corecovred_txt = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov.txt",
#         core_covR = "scripts/core_cov.R",
#     output:
#         corecovred_pdf = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov_coord.pdf",
#         corecovred_txt = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov_coord.txt"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-1:10:00"),
#     resources:
#         mem_mb = 4000
#     conda: "envs/core-cov-env.yaml"
#     log: "logs/{phylotype}_core_cov_red_plots.log"
#     threads: 1
#     shell:
#         """
#         {input.core_covR} {input.corecovred_txt};
#         """
#
# CHROMS = get_chroms("database/"+DBs["microbiome"])
#
# rule generate_freebayes_regions:
#     input:
#         ref_idx = "database/"+DBs["microbiome_red"],
#         index = "database/"+DBs["microbiome_red"]+".fai",
#         bams = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam", sample=config["SAMPLES"]),
#         bam_files_index = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam.bai", sample=config["SAMPLES"]),
#         script = "scripts/fasta_generate_regions.py",
#     output:
#         regions = expand("database/freebayes_regions/"+DBs["microbiome_red"]+".{chrom}.region.{i}.bed", chrom=CHROMS, i = CHUNKS)
#     log:
#         "logs/generate_freebayes_regions.log"
#     params:
#         chunks = CHUNKS,
#         nchunks = N,
#         chroms = CHROMS,
#         outdir = "database/freebayes_regions/"+DBs["microbiome_red"],
#     shell:
#         """
#         python3 {input.script} {input.index} {params.nchunks} --chunks --bed {params.outdir} --chromosome {params.chroms}
#         """
#
# rule freebayes_profiling:
#     input:
#         bam_files = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam", sample=config["SAMPLES"]),
#         bam_files_index = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam.bai", sample=config["SAMPLES"]),
#         db_red_ind = "database/"+DBs["microbiome_red"]+".fai",
#         db_red = "database/"+DBs["microbiome_red"],
#         regions = "database/freebayes_regions/"+DBs["microbiome_red"]+".{chrom}.region.{i}.bed"
#     output:
#         vcf = temp("07_SNVProfiling/VCFs/{chrom}/variants.{i}.vcf")
#     params:
#         help_script = os.path.join(os.getcwd(), "scripts/fasta_generate_regions.py"), # if not in path
#         freebayes = "freebayes", # path
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("1-00:00:00"),
#     resources:
#         mem_mb = 10000
#     threads: 15
#     conda: "envs/snv-env.yaml"
#     log: "logs/freebayes_profiling_{chrom}_{i}.log"
#     shell:
#         """
#         freebayes --pooled-continuous --min-alternate-count 5 --min-alternate-fraction 0.1 --min-coverage 10 -f {input.db_red} -t {input.regions} {input.bam_files} > {output.vcf}
#         """
#
# rule concat_vcfs:
#     input:
#         calls = expand("07_SNVProfiling/VCFs/{chrom}/variants.{i}.vcf", chrom=CHROMS, i=CHUNKS)
#     output:
#         vcf = "07_SNVProfiling/freebayes_raw.vcf"
#     log:
#         "logs/concat_vcfs.log"
#     conda:
#         "envs/snv-env.yaml"
#     threads: 15
#     shell:
#         "bcftools concat {input.calls} | vcfuniq > {output.vcf}"
#
# rule vcf_filter_indels:
#     input:
#         vcf_raw = rules.concat_vcfs.output.vcf
#     output:
#         vcf_filt = "07_SNVProfiling/All_CandSNVs.vcf"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     conda: "envs/snv-env.yaml"
#     log: "logs/vcf_filter_indels.log"
#     shell:
#         """
#         bcftools view --exclude-types indels {input.vcf_raw} > {output.vcf_filt}
#         """
#
# rule filt_bedfiles:
#     input:
#         metafile = "database/"+DBs["microbiome"]+"_metafile.txt",
#         sdp_coords = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_corecov_coord.txt", sdp=get_g_sdp_dict(DBs["microbiome"]).keys()),
#         sdp_corecovs = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_corecov.txt", sdp=get_g_sdp_dict(DBs["microbiome"]).keys()),
#         bedfile = "database/bed_files/{genome}.bed",
#         ortho_dir = "database/OrthoFiles_"+DBs["microbiome_red"]
#     output:
#         core_red_bed = "07_SNVProfiling/bed_files/{genome}_core_filt.bed"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 2000
#     log: "logs/{genome}_filt_bedfiles.log"
#     script:
#         "scripts/filt_core_bed.py"
#
# rule calc_core_lengths:
#     input:
#         core_red_beds = expand("07_SNVProfiling/bed_files/{genome}_core_filt.bed", genome=get_g_list(DBs["microbiome_red"])),
#         db = rules.subset_ortho_and_db.output.db, # so that it is repeated if it changes - not used by script
#     output:
#         lengths = "07_SNVProfiling/table_core_length.txt"
#     params:
#         genomes_sdp_dict=get_g_sdp_dict(DBs["microbiome"]),
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 2000
#     log: "logs/calc_core_length.log"
#     script:
#         "scripts/calc_filt_core_length.py"
#
# rule cat_bedfiles:
#     input:
#         bedfiles = expand("07_SNVProfiling/bed_files/{genome}_core_filt.bed", genome=get_g_list(DBs["microbiome_red"]))
#     output:
#         bed_all = "07_SNVProfiling/bed_files/all_filt.bed"
#     log: "logs/cat_bedfiles.log"
#     threads: 1
#     shell:
#         """
#         cat {input.bedfiles} > {output.bed_all}
#         """
#
# rule vcf_filter_intersect:
#     input:
#         vcf_filt = rules.vcf_filter_indels.output.vcf_filt,
#         bed_all = rules.cat_bedfiles.output.bed_all
#     output:
#         vcf_intersect = "07_SNVProfiling/Intersect_CandSNVs.vcf",
#     params:
#         vcfintersect = "vcfintersect",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     threads:
#         20
#     conda: "envs/snv-env.yaml"
#     log: "logs/vcf_filter_intersect.log"
#     shell:
#         """
#         {params.vcfintersect} --bed {input.bed_all} {input.vcf_filt} | vcfbreakmulti > {output.vcf_intersect}
#         """
#
# rule vcf_split_by_sdp:
#     input:
#         vcf = rules.vcf_filter_intersect.output.vcf_intersect,
#         metafile = "database/"+DBs["microbiome"]+"_metafile.txt",
#     output:
#         vcfs = expand("07_SNVProfiling/{sdp}_all_samples.vcf", sdp=get_g_sdp_dict(DBs["microbiome"]).keys())
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:30:00"),
#     threads:
#         4
#     log: "logs/vcf_split_by_sdp.log"
#     script:
#         "scripts/vcf_split_by_sdp.py"
#
# rule vcf_filter_samples:
#     input:
#         vcf = "07_SNVProfiling/{sdp}_all_samples.vcf",
#         coord = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_corecov_coord.txt",
#     output:
#         vcf = "07_SNVProfiling/{sdp}.vcf",
#     params:
#         script = "scripts/filt_vcf_samples.py",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     threads:
#         20
#     conda: "envs/snv-env.yaml"
#     log: "logs/{sdp}_vcf_filter_samples.log"
#     shell:
#         """
#         SAMPLES=$( python3 {params.script} {input.vcf} {input.coord} )
#         echo $SAMPLES
#         if [ -z "$SAMPLES" ]; then
#             echo "No samples in list!"
#             touch {output.vcf}
#         else
#             vcfkeepsamples {input.vcf} $SAMPLES > {output.vcf}
#         fi
#         """
#
# rule vcf_filter_missing:
#     input:
#         vcf = "07_SNVProfiling/{sdp}.vcf"
#     output:
#         freq = "07_SNVProfiling/{sdp}_filt.freq",
#     params:
#         summary = "07_SNVProfiling/summary_filtering.txt",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     threads:
#         4
#     log: "logs/{sdp}_vcf_filter_missing.log"
#     script:
#         "scripts/filter_snvs.py"
#
# rule summarise_snps:
#     input:
#         metadata = "Metadata_"+ProjectIdentifier+".csv",
#         lengths = "07_SNVProfiling/table_core_length.txt",
#         freq = rules.vcf_filter_missing.output.freq,
#         coord = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_corecov_coord.txt",
#     output:
#         txt_sample = "07_SNVProfiling/{sdp}_sample_var_host.txt",
#         txt_host = "07_SNVProfiling/{sdp}_tot_var.txt",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:20:00"),
#     threads:
#         4
#     conda: "envs/mapping-env.yaml" # for perl
#     log: "logs/{sdp}_summarise_snps.log"
#     script:
#         "scripts/summarise_snps.py"
#
# rule summarise_shared_fractions:
#     input:
#         freq = rules.vcf_filter_missing.output.freq
#     output:
#         shared_frac = "07_SNVProfiling/{sdp}_filt_shared_fraction.txt",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-3:10:00"),
#     threads:
#         4
#     log: "logs/{sdp}_summarise_shared_fractions.log"
#     conda: "envs/mapping-env.yaml" # for perl
#     shell:
#         """
#         perl scripts/calc_jaccard.pl {input.freq}
#         """
#
# rule make_cumulative_curves:
#     input:
#         metadata = "Metadata_"+ProjectIdentifier+".csv",
#         freq = rules.vcf_filter_missing.output.freq,
#         lengths = rules.calc_core_lengths.output.lengths,
#     output:
#         cum_curve = "07_SNVProfiling/{sdp}_cum_curve.txt",
#         cum_curve_by_group = "07_SNVProfiling/{sdp}_cum_curve_by_group.txt",
#         cum_curve_by_colony = "07_SNVProfiling/{sdp}_cum_curve_by_colony.txt",
#     params:
#         nb_curves = 10,
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:20:00"),
#     threads:
#         4
#     conda: "envs/mapping-env.yaml" # for perl
#     log: "logs/{sdp}_make_cumulative_curves.log"
#     script:
#         "scripts/cum_curve_SNVs_host.py"
#
# rule make_dist_matrix:
#     input:
#         file_frac = rules.summarise_shared_fractions.output.shared_frac
#     output:
#         dist_matrix = "07_SNVProfiling/{sdp}_dist_matrix.txt",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:20:00"),
#     threads:
#         4
#     conda: "envs/rmd-env.yaml"
#     log: "logs/{sdp}_make_dist_matrix.log"
#     script:
#         "scripts/distance_matrix.R"
