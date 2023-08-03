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