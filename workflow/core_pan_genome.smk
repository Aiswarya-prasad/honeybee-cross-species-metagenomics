#!/usr/bin/env python

"""
name: core_pan_genome
description: xxx
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - xxx
scripts:
    - x*.py
targets:
    - ...
"""

rule rename_prodigal_checkm:
    input:
        prodigal_checkm_faa = "results/09_MAGs_collection/All_mags_sub/prodigal_output/from_checkm/{mag}.faa",
        prodigal_checkm_gff = "results/09_MAGs_collection/All_mags_sub/prodigal_output/from_checkm/{mag}.gff",
        mag_fa = "results/09_MAGs_collection/All_mags_sub/MAGs/{mag}.fa",
        collected = "results/09_MAGs_collection/All_mags_sub/prodigal_output/collect_from_checkm.done"
    output:
        renamed_ffn = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed_for_pangenome/{mag}/{mag}.ffn",
        renamed_faa = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed_for_pangenome/{mag}/{mag}.faa",
        renamed_gff = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed_for_pangenome/{mag}/{mag}.gff",
        renamed_bed = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed_for_pangenome/{mag}/{mag}.bed"
    params:
        sample_name = lambda wildcards: wildcards.mag.split("_")[0],
        outdir = "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed_for_pangenome",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed_for_pangenome/{mag}_rename_prodigal.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/prodigal_output/renamed_for_pangenome/{mag}_rename_prodigal.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        cat {input.prodigal_checkm_faa} | sed -e 's/ID=/ID={wildcards.mag}_/g' > {output.renamed_faa}
        cat {input.prodigal_checkm_gff} | sed -e 's/ID=/ID={wildcards.mag}_/g' > {output.renamed_gff}
        python3 scripts/gff_to_bed.py --gff {output.renamed_gff} --bed {output.renamed_bed}
        bedtools getfasta -fi {input.mag_fa} -bed {output.renamed_bed} -fo {output.renamed_ffn}
        # python3 scripts/gff_to_bed.py --gff {output.renamed_gff} --bed {output.renamed_bed} --rename
        # bedtools getfasta -fi {input.mag_fa} -bed {output.renamed_bed} -fo {output.renamed_ffn} -nameOnly
        """

# rule motupan:
#     input:
#     output:
# checkm_out="motupan_test/generic_completeness.tsv"
# cat ../06_MAG_binning/evaluate_bins/checkm_drep_summary.txt | tr "," "\t" | sed -e 's/.fa//' | sed -e 's/genome/Bin Id/' | sed -e 's/comp/Comp/' | sed -e 's/cont/Cont/' > ${checkm_out}
# input_faa_dir="../database/MAGs_database_Orthofinder/g__Lactobacillus/"
# out_dir="motupan_test/g__Lactobacillus"
# mkdir -p ${out_dir}
# input_orthofile = "/scratch/aprasad/211018_Medgenome_india_samples/database/MAGs_database_Orthofinder/g__Lactobacillus/OrthoFinder/Results_g__Lactobacillus/Orthogroups/Orthogroups.GeneCount.tsv"
# output_cogfile = "/scratch/aprasad/211018_Medgenome_india_samples/15_FurtherProcessing/motupan_test/g__Lactobacillus/motupan_cog.tsv" 
# OG_file="motupan_test/g__Lactobacillus/motupan_cog.tsv" # using scripts/prepare_motupan_OG_file.py
# mOTUpan.py --gene_clusters_file ${OG_file} --boots 100 -o ${out_dir}/mOTUpan.tsv --checkm ${checkm_out} | tee ${out_dir}/mOTUpan.log

# rule summarise_motupan