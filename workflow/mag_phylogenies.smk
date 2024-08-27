#!/usr/bin/env python

"""
name: mag_phylogenies
description: xxx
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - (checkpoint) make_phylo_metadata
    - (checkpoint) collect_prodigal_from_checkm
    - download_assembly_summary
        + download assembly_summary frm NCBI
    - rename_prodigal_checkm:
        + rename prodigal genes from checkm output for MAGs
    - prepare_genomes_faa:
        + prepare genomes for orthofinder by running prodigal on them if needed and renaming the genes
    - rename_faa_and_ffn:
        + rename the genes in the faa and ffn files for orthofinder with simple headers
    - collect_faa_orthofinder:
        + collect the faa files for orthofinder from the renamed faa files for each MAG
    - run_orthofinder:
        + run orthofinder on the faa files for each genus to get orthogroups
    - run_orthofinder_iqtree:
        + run orthofinder with iqtree on the orthogroups to get a species tree
    - get_OG_nuc_sequences:
        + get the nucleotide sequences for each orthogroup
    - extract_bac120_nucleotide:
        + extract the nucleotide sequences for the bac120 markers using the script (not a rule)
    - align_bac120_nucleotide_macse:
        + align the bac120 nucleotide sequences using macse
    - replace_unknown_characters:
        + replace unknown characters in the aligned sequences with N for iqtree
    - make_bac120_nucleotide_tree:
        + make a tree from the aligned bac120 nucleotide sequences using iqtree
    - dram_annotate_mags:
        + annotate the MAGs using DRAM
    - dram_rename_annotations:
        + rename the annotations from DRAM to be more informative and easier to use in downstream analysis
    - dram_annotate_concat_mags:
        + concatenate the DRAM annotations for all MAGs into a single file for running distill
    - dram_distill_mags:
        + distill the DRAM annotations to get a summary of the metabolism of the MAGs
scripts:
    scripts/make_phylo_metadata.py
    scripts/collect_prodigal_from_checkm.sh
    scripts/gff_to_bed.py
    scripts/extract_marker_nucleotide_sequences.py
"""
# orthofinder and motupan

checkpoint make_phylo_metadata:
    input:
        mag_metadata = lambda wildcards: checkpoints.mag_metadata_summary.get().output.metadata,
        Isolate_metadata = "config/Outgroup_isolate_genomes.tsv",
        assembly_summary = "results/11_phylogenies/assembly_summary.txt"
    output:
        phylo_metadata = "results/11_phylogenies/phylo_genomes_metadata.tsv",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_phylo_metadata",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 2
    log: "results/11_phylogenies/phylo_genomes_metadata.log"
    # benchmark:
    # conda: 
    shell:
        """
        python3 scripts/make_phylo_metadata.py \
            --mag_metadata {input.mag_metadata} \
            --Isolate_metadata {input.Isolate_metadata} \
            --phylo_metadata {output.phylo_metadata}
        """

rule download_assembly_summary:
    output:
        assembly_summary = "results/11_phylogenies/assembly_summary.txt"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="download_assembly_summary",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 2
    log: "results/11_phylogenies/download_assembly_summary.log"
    shell:
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -O {output.assembly_summary}
        """

# this is a checkpoint becuase it does the collection for all the MAGs for which checkm already ran prodigal
# rules using the results of prodigal for mags should only be run after this rule
# puts the checkm outputs at "outdir"
checkpoint collect_prodigal_from_checkm:
    input:
        # ensure checkm is run for all the samples from which MAGs were made
        checkm_merged = "results/09_MAGs_collection/checkm_merged.tsv"
    output:
        collected = "results/09_MAGs_collection/prodigal_output/collect_from_checkm.done"
    params:
        # pattern to use in script to collect prodigal genes from checkm
        # checkm_prodigal_genes = "results/07_MAG_binng_QC/03_checkm_results/*/bins/*",
        outdir = "results/09_MAGs_collection/prodigal_output/from_checkm",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/prodigal_output/collect_prodigal_from_checkm.log"
    benchmark: "results/09_MAGs_collection/prodigal_output/collect_prodigal_from_checkm.benchmark"
    shell:
        """
        bash scripts/collect_prodigal_from_checkm.sh \
            {params.outdir} \
            {output.collected}
        """

rule rename_prodigal_checkm:
    input:
        # prodigal_checkm_faa = lambda wildcards: f"results/09_MAGs_collection/prodigal_output/from_checkm/{wildcards.mag}.faa",
        # prodigal_checkm_gff = lambda wildcards: f"results/09_MAGs_collection/prodigal_output/from_checkm/{wildcards.mag}.gff",
        mag_fa = "results/09_MAGs_collection/MAGs/{mag}.fa",
        collected = lambda wildcards: checkpoints.collect_prodigal_from_checkm.get().output.collected
        # results/09_MAGs_collection/prodigal_output/collect_from_checkm.done
    output:
        renamed_ffn = "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{mag}/{mag}.ffn",
        renamed_faa = "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{mag}/{mag}.faa",
        renamed_gff = "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{mag}/{mag}.gff",
        renamed_bed = "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{mag}/{mag}.bed"
    params:
        # sample_name = lambda wildcards: wildcards.mag.split("_")[0],
        outdir = "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{mag}_rename_prodigal.log"
    benchmark: "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{mag}_rename_prodigal.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        collected={input.collected}
        collected_dir=${{collected/collect_from_checkm.done/from_checkm}}
        cat ${{collected_dir}}/{wildcards.mag}.faa | sed -e 's/ID=/ID={wildcards.mag}_/g' > {output.renamed_faa}
        cat ${{collected_dir}}/{wildcards.mag}.gff | sed -e 's/ID=/ID={wildcards.mag}_/g' > {output.renamed_gff}
        python3 scripts/gff_to_bed.py --gff {output.renamed_gff} --bed {output.renamed_bed}
        bedtools getfasta -fi {input.mag_fa} -bed {output.renamed_bed} -fo {output.renamed_ffn}
        """
# needs them to be at "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/.."
# as expected by get_genome_path()
# if GCA... fails, edit the config/Outgroup_isolate_genomes.tsv file as needed and continue from there ->
# often is it the last didgit after the "." that is different
# phylometadata will be remade
rule prepare_genomes_faa:
    input:
        phylo_genomes_metadata = "results/11_phylogenies/phylo_genomes_metadata.tsv",
        assembly_summary = "results/11_phylogenies/assembly_summary.txt",
        genome_file_faa = lambda wildcards: get_genome_path(wildcards.genome, checkpoints.mag_metadata_summary.get().output.metadata, "faa")["path"],
        genome_file_ffn = lambda wildcards: get_genome_path(wildcards.genome, checkpoints.mag_metadata_summary.get().output.metadata, "ffn")["path"]
    output:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.faa",
        genome_ffn = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.ffn",
    params:
        input_type = lambda wildcards: get_genome_path(wildcards.genome, checkpoints.mag_metadata_summary.get().output.metadata, "faa")["type"],
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{genome}_prepare_genomes_faa",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 2
    log: "results/11_phylogenies/00_genomes/{genome}/{genome}_prepare_genomes.log"
    benchmark: "results/11_phylogenies/00_genomes/{genome}/{genome}_prepare_genomes.benchmark.json"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        genome_ffn={output.genome_ffn}
        genome_gff=${{genome_ffn/.ffn/.gff}}
        if [ {params.input_type} == "fna" ]
        then
            echo "running prodigal as fna was returned"
            # remember genoe_file_faa is actually fna here!
            prodigal -i {input.genome_file_faa} -a {output.genome_faa} -o ${{genome_gff}} -f gff
            bedtools getfasta -fi {input.genome_file_faa} -bed ${{genome_gff}} -fo ${{genome_ffn}}
        else
            cp {input.genome_file_faa} {output.genome_faa}
            echo "copied faa file to final destination"
            cp {input.genome_file_ffn} {output.genome_ffn}
            echo "copied ffn file to final destination"
        fi
        """

rule rename_faa_and_ffn:
    input:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.faa",
        genome_ffn = "results/11_phylogenies/00_genomes/{genome}/{genome}_original.ffn",
    output:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}.faa",
        genome_ffn = "results/11_phylogenies/00_genomes/{genome}/{genome}.ffn",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_phylo_metadata",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 4
    log: "results/11_phylogenies/00_genomes/{genome}/{genome}_rename_faa_and_ffn.log"
    run:
        i = 1
        with open(output.genome_faa, "w") as fh,  open(f"{input.genome_faa}") as f:
            for line in f:
                if line.startswith(">"):
                    fh.write(f">{wildcards.genome}_{i}\n")
                    i += 1
                else:
                    fh.write(line)
        j = 1
        with open(output.genome_ffn, "w") as fh,  open(f"{input.genome_ffn}") as f:
            for line in f:
                if line.startswith(">"):
                    fh.write(f">{wildcards.genome}_{j}\n")
                    j += 1
                else:
                    fh.write(line)
        if i != j:
            # consider adding additional checks to make sure the faa and fna record names are correctly matched
            sys.exit("Number of genes in faa and ffn files are not equal for {wildcards.genome}")


rule collect_faa_orthofinder:
    input:
        genome_faa = "results/11_phylogenies/00_genomes/{genome}/{genome}.faa"
    output:
        genome_faa = "results/11_phylogenies/01_orthofinder_input/{genus}/{genome}.faa",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("16G")
    threads: 2
    log: "results/11_phylogenies/01_orthofinder_input/{genus}_collect_input_logs/{genome}_collect.log"
    benchmark: "results/11_phylogenies/01_orthofinder_input/{genus}_collect_input_logs/{genome}_collect.benchmark"
    run:
        shell("cp {input.genome_faa} {output.genome_faa}")


rule run_orthofinder:
    input:
        faa_files = lambda wildcards: expand("results/11_phylogenies/01_orthofinder_input/{{genus}}/{genome}.faa", genome=get_mags_for_genus_phylogeny(wildcards.genus, checkpoints.make_phylo_metadata.get().output.phylo_metadata)),
    output:
        trees = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Species_Tree/SpeciesTree_rooted_node_labels.txt",
        orthofile = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Orthogroups/Orthogroups.txt"
    params:
        outdir = "results/11_phylogenies/02_orthofinder_results/{genus}",
        prev_outdir = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00"),
        # runtime_s=lambda wildcards: convertToSec("1-2:10:00") if wildcards.genus == "g__Lactobacillus" or wildcards.genus in ["g__Bifidobacterium", "g__Lactobacillus"] else convertToSec("0-15:10:00"),
    resources:
        mem_mb = convertToMb("200G")
    threads: 16
    log: "results/11_phylogenies/02_orthofinder_results/{genus}_run_orthofinder.log"
    benchmark: "results/11_phylogenies/02_orthofinder_results/{genus}_run_orthofinder.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        # https://github.com/davidemms/OrthoFinder#running-orthofinder
        orthofinder -ot -t {threads} -n {wildcards.genus} \
        -ft {params.prev_outdir} -ot \
        # -f $(dirname {input.faa_files[0]}) -o {params.outdir} \
        -M msa -T fasttree \
        2>&1 | tee -a {log}
        """
        # since we are only resuming now, I have edited this but do this in a more robust way later
        # ie decide if it is fresh or continued based on data in the output dir

rule run_orthofinder_iqtree:
    input:
        trees = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Species_Tree/SpeciesTree_rooted_node_labels.txt",
        orthofile = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Orthogroups/Orthogroups.txt"
    output:
        trees = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}_iqtree/Species_Tree/SpeciesTree_rooted_node_labels.txt"
    params:
        prev_outdir = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}",
        new_outdir = "results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}_iqtree",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=lambda wildcards: convertToSec("3-00:00:00"),
    resources:
        mem_mb = convertToMb("200G")
    threads: 16
    log: "results/11_phylogenies/02_orthofinder_results/{genus}_iqtree_run_orthofinder.log"
    benchmark: "results/11_phylogenies/02_orthofinder_results/{genus}_iqtree_run_orthofinder.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        rm -rf {params.new_outdir}
        orthofinder -ft {params.prev_outdir} \
            -M msa -T iqtree -ot -t {threads} \
            -n {wildcards.genus}_iqtree 2>&1 | tee {log}
        """

# identify markers in MAGs for each genus
def get_mag_and_id(gene):
    dict_return = {}
    dict_return["id"] = gene.split("_")[-1]
    dict_return["mag"] = "_".join(gene.split("_")[:-1])
    return dict_return

# # For each species that could have potentially codiversified, the pairwise bacdiv for the MAGs 
# # of that species should be lower for pairs of MAGs from the same species than for pairs of MAGs
# # from different species and the slope of the regression line should be such that this increases
# # with host divergence time and the slope should be a value that is realistic for #substitutions/site/year

# mags_to_use = {"g__Lactobacillus" : ["A4-2_9"]} # only an example

# rule get_OG_nuc_sequences:
#     input:
#         nuc_sequences
#     output:
#         nuc_sequences = "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/{marker}.fa",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         account="pengel_spirit",
#         runtime_s=convertToSec("3-00:00:00") # 25h worked for ~118/120 markers
#     resources:
#         mem_mb = convertToMb("100G")
#     threads: 8
#     log: "results/11_phylogenies/nucleotide_trees/{genus}/{marker}_get_OG_nuc_sequences.log"
#     benchmark: "results/11_phylogenies/nucleotide_trees/{genus}/{marker}_get_OG_nuc_sequences.benchmark"
#     conda: "../config/envs/phylogenies-env.yaml"
#     run:
#         og = wildcards.marker
#         og_file = input.orthofile
#         gene_id = get_gene_id_num(og)
#         mag = get_mag_and_id(gene_id)["mag"]
#         mag_faa = get_genome_path(mag, checkpoints.mag_metadata_summary.get().output.metadata, "faa")["path"]
#         mag_ffn = get_genome_path(mag, checkpoints.mag_metadata_summary.get().output.metadata, "ffn")["path"]
#         # resume from here



# rule align_bac120_nucleotide_macse:
#     input:
#         input_seq = "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/{marker}.fa"
#     output:
#         out_nuc = "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/nuc/{marker}.fa",
#         out_aa = "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/aa/{marker}.fa"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         account="pengel_spirit",
#         runtime_s=convertToSec("3-00:00:00") # 25h worked for ~118/120 markers
#     resources:
#         mem_mb = convertToMb("100G")
#     threads: 8
#     log: "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/{marker}_aligned.log"
#     benchmark: "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/{marker}_aligned.benchmark"
#     conda: "../config/envs/macse-env.yaml"
#     shell:
#         """
#         macse -prog alignSequences -seq {input.input_seq} \
#               -out_NT {output.out_nuc} -out_AA {output.out_aa} \
#               -gc_def 11 
#         """

# rule replace_unknown_characters:
#     input:
#         input_seq = "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/nuc/{marker}.fa"
#     output:
#         output_seq = "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/nuc/clean/{marker}_clean.fa"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-00:02:00")
#     resources:
#         mem_mb = convertToMb("4G")
#     threads: 8
#     log: "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/logs/{marker}_clean.log"
#     benchmark: "results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/logs/{marker}_clean.benchmark"
#     conda: "../config/envs/phylogenies-env.yaml"
#     shell:
#         """
#         # since iqtree treats all the unknown characters the same way and complains about !
#         cat {input.input_seq} | sed -e 's/!/N/g' > {output.output_seq}
#         """

def get_genus_markers(genus):
    return [os.path.basename(x).split(".fa")[0] for x in glob.glob(f"results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/nuc/clean/*.fa")]

rule make_bac120_nucleotide_tree:
    input:
        input_seqs = lambda wildcards: [f"results/11_phylogenies/nucleotide_trees/{genus}/sequences_aligned/nuc/clean/{marker}_clean.fa" for marker in get_genus_markers(wildcards.genus)]
    output:
        marker_file = "results/11_phylogenies/nucleotide_trees/iqtree/{genus}/{genus}.done"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00")
    resources:
        mem_mb = convertToMb("150G")
    threads: 8
    log: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/MAGs_bac120_nuc/MAGs_bac120_nuc.log"
    benchmark: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/MAGs_bac120_nuc/MAGs_bac120_nuc.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        outdir=$(dirname {output.out_tree})
        indir=$(dirname {input.input_seqs[0]})
        iqtree -s ${{indir}} \
            -nt {threads} \
            -bb 10000 \
            -seed 1234 \
            -m MFP \
            -pre ${{outdir}}/MAGs_bac120_nuc
        """


# rule extract_bac120_nucleotide:
    # this is done by the script
    # scripts/extract_marker_nucleotide_sequences.py
    # make this into a rule later maybe

bac120_markers = [os.path.basename(x).split(".fa")[0] for x in glob.glob("results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences/*.fa")]

rule align_bac120_nucleotide_macse:
    input:
        input_seq = "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences/{marker}.fa"
    output:
        out_nuc = "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/nuc/{marker}.fa",
        out_aa = "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/aa/{marker}.fa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00") # 25h worked for ~118/120 markers
    resources:
        mem_mb = convertToMb("100G")
    threads: 8
    log: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/{marker}_aligned.log"
    benchmark: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/{marker}_aligned.benchmark"
    conda: "../config/envs/macse-env.yaml"
    shell:
        """
        macse -prog alignSequences -seq {input.input_seq} \
              -out_NT {output.out_nuc} -out_AA {output.out_aa} \
              -gc_def 11 
        """

rule replace_unknown_characters:
    input:
        input_seq = "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/nuc/{marker}.fa"
    output:
        output_seq = "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/nuc/clean/{marker}_clean.fa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:02:00")
    resources:
        mem_mb = convertToMb("4G")
    threads: 8
    log: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/logs/{marker}_clean.log"
    benchmark: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/logs/{marker}_clean.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        # since iqtree treats all the unknown characters the same way and complains about !
        cat {input.input_seq} | sed -e 's/!/N/g' > {output.output_seq}
        """

rule make_bac120_nucleotide_tree:
    input:
        input_seqs = expand("results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/nuc/clean/{marker}_clean.fa", marker=bac120_markers)
    output:
        out_tree = "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/MAGs_bac120_nuc/MAGs_bac120_nuc.treefile"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00")
    resources:
        mem_mb = convertToMb("150G")
    threads: 8
    log: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/MAGs_bac120_nuc/MAGs_bac120_nuc.log"
    benchmark: "results/11_phylogenies/05_MAG_bac120_nucleotide_trees/MAGs_bac120_nuc/MAGs_bac120_nuc.benchmark"
    conda: "../config/envs/phylogenies-env.yaml"
    shell:
        """
        outdir=$(dirname {output.out_tree})
        indir=$(dirname {input.input_seqs[0]})
        iqtree -s ${{indir}} \
            -nt {threads} \
            -bb 10000 \
            -seed 1234 \
            -m MFP \
            -pre ${{outdir}}/MAGs_bac120_nuc
        """
    

rule dram_annotate_mags:
    input:
        dram_config = "config/dram_config.json",
        filt_faa = "results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{mag}/{mag}.faa",
    output:
        dram_annotations = "results/09_MAGs_collection/dram_output/{mag}/annotations.tsv",
    params:
        dram_outdir = lambda wildcards: os.path.join("results/09_MAGs_collection/dram_output/", wildcards.mag),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:00:00")
    resources:
        mem_mb = convertToMb("100G")
    threads: 8
    log: "results/09_MAGs_collection/dram_output/logs/{mag}_dram_annotate.log"
    benchmark: "results/09_MAGs_collection/dram_output/logs/{mag}_dram_annotate.benchmark"
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

rule dram_rename_annotations:
    input:
        dram_annotations = "results/09_MAGs_collection/dram_output/{mag}/annotations.tsv",
    output:
        dram_annotations = "results/09_MAGs_collection/dram_output/{mag}/annotations_renamed.tsv",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:15:00")
    log: "results/09_MAGs_collection/dram_output/logs/{mag}_dram_rename.log"
    benchmark: "results/09_MAGs_collection/dram_output/logs/{mag}_dram_rename.benchmark"
    run:
        with open(input.dram_annotations) as f:
            header = f.readline()
            with open(output.dram_annotations, "w") as fh:
                fh.write(header)
                for line in f:
                    # first column has no header and it contains the gene name
                    # second column has fasta and contains the mag name this is to be replaced
                    # with the handmade species name of that mag
                    # all other columns are the same
                    mag_name = line.split("\t")[1]
                    species_name = handmade_species_name(get_mag_info(mag_name, checkpoints.mag_metadata_summary.get().output.metadata))
                    fh.write(line.replace(mag_name, str(species_name)+"_"+str(mag_name)))


rule dram_annotate_concat_mags:
    input:
        dram_annotations = lambda wildcards: expand("results/09_MAGs_collection/dram_output/{mag}/annotations_renamed.tsv", mag = get_medium_mags(checkpoints.mag_metadata_summary.get().output.metadata)),
    output:
        dram_annotation_merged = "results/09_MAGs_collection/dram_distill/annotations.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:15:00")
    log: "results/09_MAGs_collection/dram_output/logs/dram_annotate_concat_mags.log"
    benchmark: "results/09_MAGs_collection/dram_output/logs/dram_annotate_concat_mags.benchmark"
    run:
        with open(input.dram_annotations[0]) as f:
            header = f.readline()
        with open(output.dram_annotation_merged, "w") as fh:
            fh.write(header)
            for file in input.dram_annotations:
                with open(file, "r") as f:
                    next(f)
                    for line in f:
                        fh.write(line)


rule dram_distill_mags:
    input:
        dram_config = "config/dram_config.json",
        dram_annotation_merged = "results/09_MAGs_collection/dram_distill/annotations.tsv",
    output:
        dram_distilled = "results/09_MAGs_collection/dram_distill/output/metabolism_summary.xlsx",
    params:
        dram_outdir = lambda wildcards: os.path.join("results/09_MAGs_collection/dram_distill/output/"),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:00:00")
    resources:
        mem_mb = convertToMb("100G")
    threads: 8
    log: "results/09_MAGs_collection/dram_output/logs/dram_distill_mags.log"
    benchmark: "results/09_MAGs_collection/dram_output/logs/dram_distill_mags.benchmark"
    conda: "../config/envs/dram-env.yaml"
    shell:
        """
        rm -rf {params.dram_outdir}
        DRAM.py distill -i {input.dram_annotation_merged} -o {params.dram_outdir} --config_loc {input.dram_config}
        """

"""
mkdir results/09_MAGs_collection/functions_list/
cat results/09_MAGs_collection/dram_distill/annotations.tsv | cut -f2,4 | grep -P "\tK" > results/09_MAGs_collection/functions_list/all_kos.txt
"""

# for each mag, make a list of KOs from non-empty lines of annotation
# and name the file with the mag name containing handmade species name
# this will be used for ko_mapper.py and minpath

"""
did minpath and ko_mapper outside of snakemake
"""
