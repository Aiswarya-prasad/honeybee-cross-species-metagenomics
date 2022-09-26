#!/usr/bin/env python3
import os

# The script will print the total number of candidate SNVs found to be polymorphic, and the number of positions that were skipped due to missing data.

input_vcfs = snakemake.input["vcf"] + snakemake.input["vcf_filt"]
output = os.path.join(os.getcwd(), snakemake.output["summary"])

for vcf in input_vcfs:
    vcf = os.path.join(os.getcwd(), vcf)
    with open(output, "w+") as out_fh:
        with open(vcf, "r") as vcf_fh:
            num_polymorphic_sites = 0
            num_skipped = 0
            for line in vcf_fh:
                if "##" in line:
                    continue
                if line.split("\t")[0] == "#CHROM":
                    samples = line.split("\t")[9:]
                    out_fh.write("SNP_info")
                    for sample in samples:
                        out_fh.write(f"\t{sample}")

                    out_fh.write("\n")
                else:
                    freqs = line.split("\t")[9:]
