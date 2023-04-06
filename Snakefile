#!/usr/bin/env python

"""
name: honeybee-MAGs-pipeline
description: Snakefile to be used to launch the pipeline and run qc, assembly, binning and some downstream steps
author: Aiswarya Prasad (aiswarya.prasad@unil.ch)
dependencies:
    - trim-qc
"""

import os
import glob
import itertools

configfile: "config/config.yaml"
# read config information into local variables to improve readability of rules
if config["LocalBackup"]:
    localrules: backup

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

def get_all_read_raw_files():
    """
    returns all read files R1, R2 and all lanes
    unmerged for running through fastqc raw
    figure out how to do this with wildcards or
    how to name the fastqc output accordingly.
    Just merge it maybe...
    """
    RawdataPath = config["RawdataPaths"]
    file_paths = []
    for path in RawdataPath:
        file_paths = list(itertools.chain(file_paths, [file for file in glob.glob(os.path.join(path, "*.fastq.gz"))]))
    return(file_paths)


def get_read_raw_files(sample, read):
    """
    raw data comes from multiple lanes and from across multiple runs
    so after fastqc is run on raw reads seperately, the files of the
    same sample from across runs and lanes will be concatenated and
    deposited in the 00_RawData directory. The path to the files will
    be for now, that on the scratch directory (copied in) but next week
    it will be be changed to the path as per their location in the 
    NGS_data directory in nas dcsr
    returns the list of files for each sample across runs and lanes
    sample is name of the sample and read should be either R1 or R2
    """
    RawdataPath = config["RawdataPaths"]
    read_pattern = "_" + read
    if read != "R1" or read != "R2":
        return(None)
    for path in RawdataPath:
        file_paths = [file for file in glob.glob(os.path.join(path, "*.fastq.gz")) if sample in file and read_pattern in file]
    return(file_paths)

targets:
    html_qc_raw=expand("results/fastqc/raw/{sample}_{read}_fastqc.html"),
    html_qc_trim="fastqc/trim/{sample}_{read}_trim_fastqc.html",

include: "rules/time-qc.smk"
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