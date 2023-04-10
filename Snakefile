#!/usr/bin/env python

"""
name: honeybee-MAGs-pipeline
description: Snakefile to be used to launch the pipeline and run qc, assembly, binning and some downstream steps
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
dependencies:
    - raw data
    - config files
"""

import os
import glob
import itertools

configfile: "config/config.yaml"
# read config information into local variables to improve readability of rules
if config["LocalBackup"]:
    localrules: backup

SAMPLES = config["SAMPLES_KE"] + config["SAMPLES_INDIA"] + config["SAMPLES_MY"]

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

# def get_wildcard_info(path):
#     split_path = [x for x in path.split("/") if "_" in x and x not in ["NGS_data", "general_data"]]
#     run_name = split_path[0]
#     file_name = split_path[-1].strip(".fastq.gz")
#     if "_L" in file_name:
#         lane = [x for x in file_name.split("_") if x.startswith("L") and len(x) == 2][0]
#         read = [x for x in file_name.split("_") if x.startswith("R") and len(x) == 2][0]
#         sample = file_name.split("_L")[0]
#     else:
#         lane = ""
#         read = [x for x in file_name.split("_") if x.startswith("R") and len(x) == 2][0]
#         sample = file_name.split("_R")[0]
#     to_return = {"run": run_name, "read": read, "sample": sample, "lane": lane}    
#     return(to_return)

def get_input_file(sample, read, lane, run):
    """
    runs are idetified by their date (first 8 digits) - if this is different edit function accordingly
    """
    split_path = [x for x in path.split("/") if "_" in x and x not in ["NGS_data", "general_data"]]
    run_name = split_path[0]
    file_name = split_path[-1].strip(".fastq.gz")
    if "_L" in file_name:
        lane = [x for x in file_name.split("_") if x.startswith("L") and len(x) == 2][0]
        read = [x for x in file_name.split("_") if x.startswith("R") and len(x) == 2][0]
        sample = file_name.split("_L")[0]
    else:
        lane = ""
        read = [x for x in file_name.split("_") if x.startswith("R") and len(x) == 2][0]
        sample = file_name.split("_R")[0]
    to_return = {"run": run_name, "read": read, "sample": sample, "lane": lane}    
    return()



onstart:
    shell("python3 scripts/make_reads_list_file.py")

raw_paths_dict = yaml.safe_load(open("config/raw_file_paths.yaml", "r"))

targets:
    html_qc_raw = expand("results/00_rawreads/fastqc/{sample}_{lane}_{read}_{run}_fastqc.html", sample = SAMPLES),
    html_qc_trim = expand("results/00_trimmedreads/fastqc/{sample}_{lane}_{read}_{run}_fastqc.html", sample = SAMPLES),
    trimmed_reads = expand("results/00_trimmedreads/{sample}_{lane}_{read}_{run}.fastq.gz", ...),
    concat_reads = expand("results/01_trimmedconcatreads/{sample}_{read}.fastq.gz", ...)

include: "rules/trim-qc.smk"
include: "rules/motus-profiling.smk"
include: "rules/assemble-qc.smk"

# to do next, 
# binning
# drep and gtdb ..
# rename MAGs (with completeness and genus in the name?) and checkpoint table
# annotation
# phylogenies
# core coverage
# sdp validation
# popcogent