#!/usr/bin/env python

"""
name: backmapping-binning
description: Takes binning results from metabat2 and summarises it then runs checkm, gtdbtk, drep on mags and creates a filtered mag database for instrain (non-redundant) and redundant for mapping to
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - checkm_evaluation
scripts:
    - *.py
targets:
    - ...
"""

checkpoint collect_mags:
    input:
        bins_dir = expand("results/07_MAG_binng_QC/02_bins/{sample}/", sample=SAMPLES_INDIA+SAMPLES_MY)
    output:
        collect_mags_marker = "results/09_MAGs_collection/MAGs/collect_mags.done"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/MAGs/collect_mags.log"
    benchmark: "results/09_MAGs_collection/MAGs/collect_mags.benchmark"
    shell:
        """
        mags_path={output.collect_mags_marker}
        mags_path=${{mags_path/collect_mags.done/}}
        for dir in {input.bins_dir}; do
            cp ${{dir}}/*.fa ${{mags_path}}
        done
        echo "renaming MAG files - replacing . with _"
        for mag in ${{mags_path}}/*.fa; do
            mv ${{mag}} ${{mag/\./_}}
        done
        rm ${{mags_path}}/*unbinned.fa
        touch {output.collect_mags_marker}
        """

rule checkm_evaluate:
    input:
        bins = "results/07_MAG_binng_QC/02_bins/{sample}/"
    output:
        checkm_summary = "results/07_MAG_binng_QC/03_checkm_results/{sample}_checkm.summary",
        lineage_ms = "results/07_MAG_binng_QC/03_checkm_results/{sample}/lineage.ms",
        plots_marker = touch("results/07_MAG_binng_QC/03_checkm_results/{sample}/plots.done"),
    params:
        extension = "fa",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="checkm_evaluate_{sample}",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:30:00"),
    resources:
        mem_mb = convertToMb("200G")
    threads: 4
    log: "results/07_MAG_binng_QC/03_checkm_results/{sample}_checkm_evaluate.log"
    benchmark: "results/07_MAG_binng_QC/03_checkm_results/{sample}_checkm_evaluate.benchmark"
    conda: "../config/envs/mags-env.yaml"
    shell:
        """
        out_file={output.checkm_summary}
        out_dir=${{out_file/_checkm.summary/}}
        checkm lineage_wf -t {threads} -x {params.extension} {input.bins} ${{out_dir}} &>> {log}
        markers_file={output.lineage_ms}
        plts_dir=${{out_dir}}/plots/
        checkm qa -t {threads} -o 2 --tab_table ${{markers_file}} ${{out_dir}} -f {output.checkm_summary} &>> {log}
        checkm gc_plot -x {params.extension} {input.bins} ${{plts_dir}} 95 &>> {log}
        checkm coding_plot -x {params.extension} ${{out_dir}} {input.bins} ${{plts_dir}} 95 &>> {log}
        checkm nx_plot -x {params.extension} {input.bins} ${{plts_dir}} &>> {log}
        checkm marker_plot -x {params.extension} ${{out_dir}} {input.bins} ${{plts_dir}} || true &>> {log}
        """

rule merge_checkm_output:
    input:
        checkm_summary = expand("results/07_MAG_binng_QC/03_checkm_results/{sample}_checkm.summary", sample=SAMPLES_sub)
    output:
        checkm_merged = "results/09_MAGs_collection/checkm_merged.tsv"
    params:
        all_mags_path="results/09_MAGs_collection",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/merge_checkm_output.log"
    benchmark: "results/09_MAGs_collection/merge_checkm_output.benchmark"
    run:
        header_written = False
        for file in input.checkm_summary:
            with open(file, "r") as in_fh:
                with open(output.checkm_merged, "a") as out_fh:
                    header = in_fh.readline()
                    if not header_written:
                        out_fh.write(header)
                        header_written = True
                    for line in in_fh:
                        line_rest = line.split("\t")[1:]
                        genome = line.split("\t")[0]
                        sample = "".join(genome.split(".")[:-1])
                        mag_num = genome.split(".")[-1]
                        genome_renamed = f"{sample}_{mag_num}"
                        line_final = "\t".join([genome_renamed] + line_rest)
                        out_fh.write(line_final)

rule gtdbtk_batchfile:
    input:
        all_mags_marker = lambda wildcards: checkpoints.collect_mags.get().output.collect_mags_marker
    output:
        batchfile = "results/09_MAGs_collection/gtdb_input_batchfile.tsv"
    params:
        genomes_dir="results/09_MAGs_collection/mags",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("4G")
    threads: 2
    log: "results/09_MAGs_collection/gtdb_input_batchfile.log"
    run:
        with open(output.batchfile, "w") as out_fh:
            for (dirpath, dirnames, filenames) in os.walk(params.genomes_dir):
                for filename in filenames:
                    if not filename.endswith(".fa"):
                        continue
                    fasta_path = os.path.join(os.getcwd(), dirpath, filename)
                    genome_id = filename.split(".fa")[0]
                    out_fh.write(f"{fasta_path}\t{genome_id}\n")

rule gtdb_annotate:
    input:
        batchfile = "results/09_MAGs_collection/gtdb_input_batchfile.tsv"
    output:
        tax_info = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/gtdb_output/classify/All_mags_sub.bac120.summary.tsv",
        tax_info_ar = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/gtdb_output/classify/All_mags_sub.ar53.summary.tsv",
    params:
        path_to_db="/work/FAC/FBM/DMF/pengel/spirit/aprasad/gtdb/release214/", # new db
        prefix = "All_mags_sub",
        mash_path = "results/09_MAGs_collection/gtdb_output/mash_sketch/cli/mash_db.msh",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="gtdb_annotate",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00"),
    resources:
        mem_mb = convertToMb("512G")
    threads: 8
    conda: "../config/envs/mags-env.yaml"
    log: "results/09_MAGs_collection/gtdb_annotate.log"
    benchmark: "results/09_MAGs_collection/gtdb_annotate.benchmark"
    shell:
        """
        # set up path to database
        export GTDBTK_DATA_PATH={params.path_to_db}
        out_file={output.tax_info}
        gtdbtk classify_wf --batchfile {input.batchfile} --out_dir ${{out_file/\/classify\/gtdbtk.bac120.summary.tsv/}} \
            --extension ".fa" --write_single_copy_genes --keep_intermediates \
            --prefix {params.prefix} --mash_db {params.mash_path} \
            --cpus {threads}
        """

rule make_drep_genome_info:
    input:
        checkm_merged = "results/09_MAGs_collection/checkm_merged.tsv",
        collect_mags_marker = lambda wildcards: checkpoints.collect_mags.get().output.collect_mags_marker
    output:
        drep_genomeinfo = "results/09_MAGs_collection/drep_genome_info.tsv",
        mags_collected = touch("results/09_MAGs_collection/high_medium_mags/collect_mags.done")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_drep_genome_info",
        account="pengel_spirit",
        runtime_s=convertToSec("0-04:00:00"),
    resources:
        mem_mb = convertToMb("2G")
    threads: 4
    log: "results/09_MAGs_collection/drep_genome_info.log"
    benchmark: "results/09_MAGs_collection/drep_genome_info.benchmark"
    run:
        shell("mkdir -p results/09_MAGs_collection/high_medium_mags")
        with open(output.drep_genomeinfo, "w") as out_fh:
            out_fh.write("genome,completeness,contamination\n")
            with open(input.checkm_merged, "r") as in_fh:
                header = in_fh.readline()
                header = header.strip().split("\t")
                ind_genome = header.index("Bin Id")
                ind_comp = header.index("Completeness")
                ind_cont = header.index("Contamination")
                for line in in_fh:
                    if "unbinned" in line:
                        continue
                    genome = line.split("\t")[ind_genome]
                    completeness = line.split("\t")[ind_comp]
                    contamination = line.split("\t")[ind_cont]
                    if float(completeness) > 50 and float(contamination) < 10:
                        shell(f"cp results/09_MAGs_collection/MAGs/{genome}.fa results/09_MAGs_collection/high_medium_mags/")
                        out_fh.write(f"{genome}.fa,{completeness},{contamination}\n")

# only dereplicate high and medium quality MAGs
# as the step is computationally expensive

rule drep_dereplicate:
    input:
        drep_genomeinfo = "results/09_MAGs_collection/drep_genome_info.tsv",
        collect_mags_marker = "results/09_MAGs_collection/high_medium_mags/collect_mags.done",
    output:
        drep_S = "results/09_MAGs_collection/drep_output/data_tables/Sdb.csv",
        drep_N = "results/09_MAGs_collection/drep_output/data_tables/Ndb.csv",
        drep_M = "results/09_MAGs_collection/drep_output/data_tables/Mdb.csv",
        drep_C = "results/09_MAGs_collection/drep_output/data_tables/Cdb.csv",
        drep_B = "results/09_MAGs_collection/drep_output/data_tables/Bdb.csv",
        drep_W = "results/09_MAGs_collection/drep_output/data_tables/Wdb.csv",
        drep_Wi = "results/09_MAGs_collection/drep_output/data_tables/Widb.csv",
        drep_gI = "results/09_MAGs_collection/drep_output/data_tables/genomeInformation.csv",
    params:
        outdir = "results/09_MAGs_collection/drep_output",
        project_dir=os.getcwd(),
        overlap = 0.2, # ask Lucas why
        ani = 0.95,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="drep_dereplicate",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00"),
    resources:
        mem_mb = convertToMb("250G")
    threads: 8
    conda: "../config/envs/mags-env.yaml"
    log: "results/09_MAGs_collection/drep_dereplicate.log"
    benchmark: "results/09_MAGs_collection/drep_dereplicate.benchmark"
    shell:
        """
        marker={input.collect_mags_marker}
        bins=${{marker/collect_mags.done/}}
        out_file={output.drep_S}
        dRep dereplicate {params.outdir} -g ${{bins}}/*.fa \
            -comp 0 -con 1000 --clusterAlg average \
            --genomeInfo {input.drep_genomeinfo} \
            -sa {params.ani} -nc {params.overlap} -p {threads} --debug
        """

checkpoint mag_metadata_summary:
    input:
        gtdb = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/gtdb_output/classify/All_mags_sub.bac120.summary.tsv",
        checkm = "results/09_MAGs_collection/checkm_merged.tsv",
        drep_gI = "results/09_MAGs_collection/drep_output/data_tables/genomeInformation.csv",
        drep_Wi = "results/09_MAGs_collection/drep_output/data_tables/Widb.csv",
        drep_W = "results/09_MAGs_collection/drep_output/data_tables/Wdb.csv",
        drep_S = "results/09_MAGs_collection/drep_output/data_tables/Sdb.csv",
        drep_C = "results/09_MAGs_collection/drep_output/data_tables/Cdb.csv",
        drep_N = "results/09_MAGs_collection/drep_output/data_tables/Ndb.csv",
        drep_M = "results/09_MAGs_collection/drep_output/data_tables/Mdb.csv",
    output:
        metadata = "results/09_MAGs_collection/MAGs_metadata_summary.tsv"
    params:
        outdir = "results/09_MAGs_collection/drep_output",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_mag_metadata_summary",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("200G")
    threads: 16
    conda: "../config/envs/scripts-env.yaml"
    log: "results/09_MAGs_collection/mag_metadata_summary.log"
    benchmark: "results/09_MAGs_collection/mag_metadata_summary.benchmark"
    shell:
        """
        python scripts/make_mag_metadata_summary.py \
            --gtdb {input.gtdb} \
            --checkm {input.checkm} \
            --drep_gI {input.drep_gI} \
            --drep_Wi {input.drep_Wi} \
            --drep_S {input.drep_S} \
            --outfile {output.metadata}
        """

rule collect_prodigal_from_checkm:
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
    shell:
        """
        bash scripts/collect_prodigal_from_checkm.sh \
            {params.outdir} \
            {output.collected}
        """

rule summarize_contig_fates:
    input:
        scaffolds = "results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta",
        whokaryote_out = "results/05_assembly/contig_fates/whokaryote/{sample}/whokaryote_predictions_S.tsv",
        kaiju_out = "results/05_assembly/contig_fates/kaiju/nr/{sample}.kaiju",
        kaiju_names = "results/05_assembly/contig_fates/kaiju/nr/{sample}_names.txt",
        kaiju_names_full = "results/05_assembly/contig_fates/kaiju/nr/{sample}_fullnames.txt",
        bins_directory = "results/07_MAG_binng_QC/02_bins/{sample}",
        gtdb_bac = "results/09_MAGs_collection/gtdb_output/classify/All_mags_sub.bac120.summary.tsv",
        gtdb_ar = "results/09_MAGs_collection/gtdb_output/classify/All_mags_sub.ar53.summary.tsv",
        # kraken too
    output:
        contig_fates = "results/05_assembly/contig_fates/{sample}_contig_fates.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    shell:
        """
        python scripts/summarize_contig_fates.py \
            --scaffolds {input.scaffolds} \
            --whokaryote_out {input.whokaryote_out} \
            --kaiju_out {input.kaiju_out} \
            --kaiju_names {input.kaiju_names} \
            --kaiju_names_full {input.kaiju_names_full} \
            --bins_directory {input.bins_directory} \
            --gtdb_bac {input.gtdb_bac} \
            --gtdb_ar {input.gtdb_ar} \
            --sample {wildcards.sample} \
            --outfile {output.contig_fates}
        """



# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution

# rule checkm_all_bins:
#     input:
#         mags_marker = checkpoint.collect_mags.get().output.collect_mags_marker,
#     output:
#         checkm_comparison = touch("results/09_MAGs_collection/checkm_output/checkm_comparison.done"),
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         jobname="checkm_all_bins",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-10:00:00"),
#     resources:
#         mem_mb = convertToMb("100G")
#     threads: 20
#     log: "results/09_MAGs_collection/checkm_output/checkm_evaluation.log"
#     benchmark: "results/09_MAGs_collection/checkm_output/checkm_evaluation.benchmark"
#     shell:
#         """
#         mags_marker={input.mags_marker}
#         bins=${{mags_marker/collect_mags.done/}}
#         outdir=results/09_MAGs_collection/checkm_output
#         checkm unique ${{bins}}
#         for mag in mags;
#         do
#             tetra=${{outdir}}/${{mag}}_tetra.tsv
#             checkm tetra ${{mag}} ${{tetra}}
#             checkm outliers ${{outdir}}/${{mag}}_outliers ${{bins}} ${{tetra}} ${{outdir}}/outliers.tsv
#         done
        
#         # checkm merge bacteria.ms ./bins ./output
#         # checkm bin_compare seqs.fna ./bins1 ./bins2 bin_comparison.tsv
#         """

# make non-redundant mag database for instrain
# ! compare corecov and instrain results from the previous analysis !
