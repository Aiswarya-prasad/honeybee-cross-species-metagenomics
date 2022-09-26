#!/usr/bin/env python3

import os

# Read the vcf-file, and generate a temporary vcf-file ([SDP]"_temp.vcf"), subset for lines corresponding to the SDP indicated by the first input file
vcf_all = os.path.join(os.getcwd(), snakemake.input["vcf"])
vcf_list = snakemake.output["vcfs"]

metafile = os.path.join(os.getcwd(), snakemake.input["metafile"])

g_sdp_dict = dict()
with open(metafile, "r") as fh:
    for line in fh:
        line = line.strip()
        genome = line.split()[0]
        sdp = line.split()[2]
        g_sdp_dict[sdp] = genome

#Read the vcf-file, and create a subset vcf-file, containing data for the SDP of interest

for vcf in vcf_list:
    print("working on "+vcf)
    sdp = vcf.split("/")[1].split("_all_samples.vcf")[0]
    print("SDP: "+sdp)
    outfile = os.path.join(os.getcwd(), vcf)
    print("writing to "+outfile+"...")
    with open(vcf_all, "r") as vcf_fh:
        print("reading from "+vcf_all+"...")
        with open(outfile, "w+") as out_fh:
            for line in vcf_fh:
                if "#" in line:
                    out_fh.write(line)
                else:
                    genome = line.split("\t")[0]
                    if genome == g_sdp_dict[sdp]:
                        out_fh.write(line)
