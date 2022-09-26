#!/usr/bin/env python3

import os
import sys

# This script is used as part of the pipeline for generating vcf-files for each SDP, subset for samples with at least 20x terminus coverage. The script works as follows:
#
# Step 1: Read the file containing the terminus coverage for all samples ([SDP]"_corecov_coord.txt", and get the names of samples meeting the threshold
#
# Step 2: Read the vcf-file, and generate a temporary vcf-file ([SDP]"_temp.vcf"), subset for lines corresponding to the SDP indicated by the first input file
#
# Step 3: Print the names of the filtered samples to stdout
#
# Thus, the script will output a temporary vcf-file, and provide the input for vcflib.

input_vcf = sys.argv[1]

coord_file = os.path.join(os.getcwd(), sys.argv[2])

high_cov_samples = list()
with open(coord_file, "r") as coord_fh:
    for line in coord_fh:
        if line.startswith("SDP"):
            continue
        else:
            sample = line.split()[1]
            ter_cov = float(line.split()[2])
            if ter_cov >= 20:
                high_cov_samples.append(sample)

# formatted for vcflib's vcfkeepsamples
# for bcftools view,
# comma-seperated list with no spaces
high_cov_samples_string = ""
for sample_raw in high_cov_samples:
    if "_mapped" in sample_raw:
        sample = sample_raw.split("_mapped")[0]
    else:
        sample = sample_raw
    high_cov_samples_string = high_cov_samples_string + " " + sample
print(high_cov_samples_string)
