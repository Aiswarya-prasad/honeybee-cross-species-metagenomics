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

rule run_whokaryote:
    input:
        gff_input = "results/06_metagenomicORFs/{sample}/{sample}.gff",
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

checkpoint collect_mags:
    input:
        bins_dir = expand("results/07_MAG_binng_QC/02_bins/{sample}/", sample=SAMPLES_sub)
    output:
        collect_mags_marker = "results/09_MAGs_collection/All_mags_sub/MAGs/collect_mags.done"
    params:
        all_mags_path="results/09_MAGs_collection/All_mags_sub",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/All_mags_sub/MAGs/collect_mags.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/MAGs/collect_mags.benchmark"
    shell:
        """
        for dir in {input.bins_dir}; do
            cp ${{dir}}/*.fa {params.all_mags_path}
        done
        echo "renaming MAG files - replacing . with _"
        for mag in {params.all_mags_path}/*.fa; do
            mv ${{mag}} ${{mag/\./_}}
        done
        rm {params.all_mags_path}/*unbinned.fa
        touch {output.collect_mags_marker}
        """

rule merge_checkm_output:
    input:
        checkm_summary = expand("results/07_MAG_binng_QC/03_checkm_results/{sample}_checkm.summary", sample=SAMPLES_sub)
    output:
        checkm_merged = "results/09_MAGs_collection/All_mags_sub/checkm_merged.tsv"
    params:
        all_mags_path="results/09_MAGs_collection/All_mags_sub",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/All_mags_sub/merge_checkm_output.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/merge_checkm_output.benchmark"
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
        # all_mags_marker = "results/09_MAGs_collection/All_mags_sub/MAGs/collect_mags.done",
    output:
        batchfile = "results/09_MAGs_collection/All_mags_sub/gtdb_input_batchfile.tsv"
    params:
        genomes_dir="results/09_MAGs_collection/All_mags_sub/mags",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("4G")
    threads: 2
    log: "results/09_MAGs_collection/All_mags_sub/gtdb_input_batchfile.log"
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
        batchfile = "results/09_MAGs_collection/All_mags_sub/gtdb_input_batchfile.tsv"
    output:
        tax_info = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/gtdb_output/classify/All_mags_sub.bac120.summary.tsv",
        # tax_info_ar = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/gtdb_output/classify/All_mags_sub.ar53.summary.tsv",
    params:
        path_to_db="/work/FAC/FBM/DMF/pengel/spirit/aprasad/gtdb/release214/", # new db
        # path_to_db="/work/FAC/FBM/DMF/pengel/spirit/aprasad/gtdb/release207_v2/",
        prefix = "All_mags_sub",
        mash_path = "results/09_MAGs_collection/All_mags_sub/gtdb_output/mash_sketch/cli/mash_db.msh",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="gtdb_annotate",
        account="pengel_spirit",
        runtime_s=convertToSec("0-15:00:00"),
    resources:
        mem_mb = convertToMb("512G")
    threads: 16
    conda: "../config/envs/mags-env.yaml"
    log: "results/09_MAGs_collection/All_mags_sub/gtdb_annotate.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/gtdb_annotate.benchmark"
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
        checkm_merged = "results/09_MAGs_collection/All_mags_sub/checkm_merged.tsv",
        collect_mags_marker = lambda wildcards: checkpoints.collect_mags.get().output.collect_mags_marker
    output:
        drep_genomeinfo = "results/09_MAGs_collection/All_mags_sub/drep_genome_info.tsv",
        mags_collected = touch("results/09_MAGs_collection/All_mags_sub/high_medium_mags/collect_mags.done")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_drep_genome_info",
        account="pengel_spirit",
        runtime_s=convertToSec("0-04:00:00"),
    resources:
        mem_mb = convertToMb("2G")
    threads: 4
    log: "results/09_MAGs_collection/All_mags_sub/drep_genome_info.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/drep_genome_info.benchmark"
    run:
        shell("mkdir -p results/09_MAGs_collection/All_mags_sub/high_medium_mags")
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
                        shell(f"cp results/09_MAGs_collection/All_mags_sub/MAGs/{genome}.fa results/09_MAGs_collection/All_mags_sub/high_medium_mags/")
                        out_fh.write(f"{genome}.fa,{completeness},{contamination}\n")

rule drep_dereplicate:
    input:
        drep_genomeinfo = "results/09_MAGs_collection/All_mags_sub/drep_genome_info.tsv",
        collect_mags_marker = "results/09_MAGs_collection/All_mags_sub/high_medium_mags/collect_mags.done",
    output:
        drep_S = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Sdb.csv",
        drep_N = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Ndb.csv",
        drep_M = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Mdb.csv",
        drep_C = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Cdb.csv",
        drep_B = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Bdb.csv",
        drep_W = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Wdb.csv",
        drep_Wi = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Widb.csv",
        drep_gI = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/genomeInformation.csv",
    params:
        outdir = "results/09_MAGs_collection/All_mags_sub/drep_output",
        project_dir=os.getcwd(),
        overlap = 0.2, # ask Lucas why
        ani = 0.95,
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="drep_dereplicate",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:00:00"),
    resources:
        mem_mb = convertToMb("250G")
    threads: 8
    conda: "../config/envs/mags-env.yaml"
    log: "results/09_MAGs_collection/All_mags_sub/drep_dereplicate.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/drep_dereplicate.benchmark"
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
        gtdb = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/gtdb_output/All_mags_sub.bac120.summary.tsv",
        checkm = "results/09_MAGs_collection/All_mags_sub/checkm_merged.tsv",
        drep_gI = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/genomeInformation.csv",
        drep_Wi = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Widb.csv",
        drep_W = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Wdb.csv",
        drep_S = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Sdb.csv",
        drep_C = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Cdb.csv",
        drep_N = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Ndb.csv",
        drep_M = "results/09_MAGs_collection/All_mags_sub/drep_output/data_tables/Mdb.csv",
    output:
        metadata = "results/09_MAGs_collection/All_mags_sub/All_mags_sub_metadata_summary.tsv"
    params:
        outdir = "results/09_MAGs_collection/All_mags_sub/drep_output",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_mag_metadata_summary",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = convertToMb("200G")
    threads: 16
    conda: "../config/envs/scripts-env.yaml"
    log: "results/09_MAGs_collection/All_mags_sub/mag_metadata_summary.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/mag_metadata_summary.benchmark"
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
        checkm_merged = "results/09_MAGs_collection/All_mags_sub/checkm_merged.tsv"
    output:
        collected = "results/09_MAGs_collection/All_mags_sub/prodigal_output/collect_from_checkm.done"
    params:
        # pattern to use in script to collect prodigal genes from checkm
        # checkm_prodigal_genes = "results/07_MAG_binng_QC/03_checkm_results/*/bins/*",
        outdir = "results/09_MAGs_collection/All_mags_sub/prodigal_output/from_checkm",
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

# make table with metadata for all mags 

# add reference genomes from (gtdb?..IMB?...)
# move on to instrain then phylo and then come back here

# for an alluvial would be nice to have
# colored by Prokaryote/Eukaryote:
# Vertical variables: contig, bins/unbinned, Genus, magOTU
# could also make this per genus and per sample
# ggplot(as.data.frame(...),
#        aes(y = Freq, axis1 = contig, axis2 = bin, axis3 = Genus / Sample, axis4 = magOTU)) +
#   geom_alluvium(aes(fill = whokaryote), width = 1/12) +
#   geom_stratum(width = 1/12, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Contig", "Bin", "Genus", "magOTU"), expand = c(.025, 0.025, 0.025, 0.025)) +
#   scale_fill_brewer(type = "qual", palette = "Set1")
# more examples: https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html

# rule contig_fates:
#     input:
#     output:
#         contig_fates = "results/05_assembly/contig_fates/"

# write rule to summarize contig fates and stats for each bin and assembly (check if this is the right way to do it) 
# rule summarize_contig_fates:
#     input:
#         bins = directory("results/07_MAG_binng_QC/02_bins/{sample}/") 
#     output:
#         contig_fates = "results/07_MAG_binng_QC/02_bins/{sample}_contig_fates/{sample}_contig_fates.txt",
#         contig_stats = "results/07_MAG_binng_QC/02_bins/{sample}_contig_stats/{sample}_contig_stats.txt"
#     params:
#         mailto="
# rules for checkm, gtdbtk, and drep
# "/work/FAC/FBM/DMF/pengel/spirit/aprasad/gtdb/release214"

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
# make redundant magdatabase for core coverage estimation
