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
    conda: "../config/envs/snv-env.yaml"
    shell:
        """
        if [ ! -f {output.scaffolds_ffn} ]; then
            prodigal -i {input.scaffolds} -o {output.scaffolds_gff} -f gff -a {output.scaffolds_faa} -d {output.scaffolds_ffn} -p meta &> {log}
        fi
        python scripts/filt_orfs.py --ffn_in {output.scaffolds_ffn} --ffn_out {output.orfs} --sample {wildcards.sample} --log {output.filt_log}
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

# to do next
# cluster genes...?
# rule to get coverage matrix of genes in all orfs