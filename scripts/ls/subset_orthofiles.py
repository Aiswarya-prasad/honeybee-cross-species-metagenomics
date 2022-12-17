#!/usr/bin/env python3

import os
import argparse

def parse_prokka_MAG_header(header):
    """
    read header with prokka prefix and contig number
    and get MAG name eg. MAG_C1.1_12 from gnl|Prokka|MAG_C1.1_12_1
    """
    "_".join(header.split("|")[-1].split("_")[:-1])


parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--input', metavar="input", required=True, help="Database of concatenated mags", action="store")
requiredNamed.add_argument('--ref_info',metavar="ref_info",required=True, help="Info about MAG status with MAG name in first column and 1 0r 0 denating if it is ref id in the third (tsv format)", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
requiredNamed.add_argument('--output', metavar='output', required=True, help="Space separated list of orthofiles to be filtered", action="store")
requiredNamed.add_argument('--group', metavar='group', required=True, help="Genus group for which orthogroups were inferred and for which to consider rep genomes", action="store")
args = parser.parse_args()

input_orthofile = args.input
ref_info = args.ref_info
output_orthofile = args.output
group_provided = args.group

rep_genomes = set()
with open(ref_info, "r") as ref_info_fh:
    header = ref_info_fh.readline()
    for line in ref_info_fh:
        line = line.strip()
        genome_id = line.split("\t")[0]
        cluster = line.split("\t")[1]
        rep_genome_status = int(line.split("\t")[2])
        group = line.split("\t")[3]
        if rep_genome_status == 1:
            rep_genomes.add(genome_id)

with open(input_orthofile, "r") as fh_orthofile:
    print(f"Writing to {output_file}")
    with open(output_orthofile, "w") as fh_orthofile_out:
        for line in fh_orthofile:
            line = line.strip()
            split_line = line.split()
            OG_id = split_line.pop(0)
            out_string = OG_id
            for gene in split_line:
                genome = gene.split("_")[0]
                if genome in rep_genomes:
                    out_string = out_string + " " + gene
                else:
                    pass
            out_string = out_string + "\n"
            fh_orthofile_out.write(out_string)
