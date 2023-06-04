#!/usr/bin/env python

"""
name: mag_db_instrain
description: Takes binning results from metabat2 and summarises it then runs checkm, gtdbtk, drep on mags and creates a filtered mag database for instrain (non-redundant) and redundant for mapping to
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
rules:
    - checkm_evaluation
scripts:
    - *.py
targets:
    - ...
"""

rule make_mag_rep_database:
    input:
        collect_mags_marker = "results/09_MAGs_collection/All_mags_sub/MAGs/collect_mags.done",
        mag_metadata_summary = "results/09_MAGs_collection/All_mags_sub/MAGs/mag_metadata_summary.tsv"
    output:
        mag_rep_database = "results/09_MAGs_collection/All_mags_sub/MAGs_rep_db/mag_rep_database.fa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    threads: 4
    log: "results/09_MAGs_collection/All_mags_sub/MAGs_rep_db/mag_rep_database.log"
    benchmark: "results/09_MAGs_collection/All_mags_sub/MAGs_rep_db/mag_rep_database.benchmark"
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python scripts/make_mag_rep_database.py \
                --collect_mags_marker {input.collect_mags_marker} \
                --mag_metadata_summary {input.mag_metadata_summary} \
                --mag_rep_database {output.mag_rep_database}
        """

rule bwa_index_rep_db:
    input:

