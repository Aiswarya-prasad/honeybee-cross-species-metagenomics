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
    benchmark: "results/05_assembly/contig_fates/whokaryote/{sample}_whokaryote.log"
    conda: "../config/envs/mags-env.yaml"
    shell:
        """
        whokaryote.py --outdir {params.outdir} --gff {input.gff_input} --model S &> {log}
        """

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
# checkpoint collect_mags:
#     input:
        

# rule checkm_evaluation:
#     input:
#         all_mags = ("06_MAG_binning/bins/{sample}"),
#     output:
#     # coding_density_plots = "06_MAG_binning/evaluate_bins/{sample}/plots/{MAG(s)}.coding_density_plots.png",
#         # gc_plots = "06_MAG_binning/evaluate_bins/{sample}/plots/{MAG(s)}.gc_plots.png",
#         # marker_pos_plots = "06_MAG_binning/evaluate_bins/{sample}/plots/{MAG(s)}.marker_pos_plot.png",
#         # nx_plots = "06_MAG_binning/evaluate_bins/{sample}/plots/{MAG(s)}.nx_plots.png",
#         plots_marker = touch("06_MAG_binning/evaluate_bins/{sample}/plots.done"),
#         summary_extended = "06_MAG_binning/evaluate_bins/{sample}_checkm.summary_extended",
#         checkm_summary = "06_MAG_binning/evaluate_bins/{sample}_checkm.summary",
#         lineage_ms = "06_MAG_binning/evaluate_bins/{sample}/lineage.ms",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         # jobname="give_name",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-02:00:00"),
#     resources:
#         mem_mb = convertToMb("50G") # needs at least roughly 40G
#     threads: 8
#     conda: "envs/mags-env.yaml"
#     log: "logs/{sample}_checkm_evaluation.log"
#     benchmark: "logs/{sample}_checkm_evaluation.benchmark"
#     shell:
#         """
#         all_mags_marker={input.all_mags_marker}
#         bins=${{all_mags_marker/\.done*/}}
#         out_file={output.checkm_summary}
#         checkm lineage_wf -x fa ${{bins}} ${{out_file/_checkm.summary/}} --threads {threads} -f {output.checkm_summary} --tab_table
#         # checkm lineage_wf -x fa ${{bins}} ${{out_file/_checkm.summary/}} --threads {threads} -f {output.checkm_summary} --tab_table --aai_strain 0.95
#         all_mags_marker={input.all_mags_marker}
#         bins=${{all_mags_marker/\.done*/}}
#         out_file={input.checkm_summary}
#         out_dir=${{out_file/_checkm.summary/}}
#         plts_dir=${{out_dir}}/plots/{wildcards.sample}/
#         markers_file={input.lineage_ms}
#         analyze_dir=${{markers_file/lineage.ms/}}
#         checkm qa -t 24 -o 2 --tab_table ${{markers_file}} ${{out_dir}} -f {output.summary_extended}
#         checkm gc_plot -x {params.extension} ${{bins}} ${{plts_dir}} 95
#         checkm coding_plot -x {params.extension} ${{out_dir}} ${{bins}} ${{plts_dir}} 95
#         checkm nx_plot -x {params.extension} ${{bins}} ${{plts_dir}}
#         checkm marker_plot -x {params.extension} ${{out_dir}} ${{bins}} ${{plts_dir}}
#         """

# summarize magOTUs and rename headers as needed

# annotate mags

# orthofinder and motupan

# make non-redundant mag database for instrain
# ! compare corecov and instrain results from the previous analysis !
# make redundant magdatabase for core coverage estimation
