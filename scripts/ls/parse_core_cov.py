#!/usr/bin/env python3

import os

# parse it to look like
# columns: sdp1; sdp2 ... and so on
# sample in rows

# make 2 output files one sample vs sdp for cov other for PTR

# for orthofile in os.listdir(snakemake.input["orthodir"]):

# true if file is there and in tact else return false
def file_checked(file):
    try:
        fh = open(file, "r")
        fh.close()
        return(True)
    except:
        print("Cant find/open this file: ", file)
        print("Not including in table")
        return(False)

cov_out = os.path.join(os.getcwd(), snakemake.output["cov_csv"])
ptr_out = os.path.join(os.getcwd(), snakemake.output["ptr_csv"])

input_files = [os.path.join(os.getcwd(), x) for x in \
               snakemake.input["txt"]]
# make a list of all the sdps and samples in the analysis

sdps = list()
samples = list()

for file in input_files:
    if file_checked(file):
        with open(file, "r") as fh_in:
            for line in fh_in:
                splits = line.split()
                if splits[0] == "SDP":
                    pass
                else:
                    if splits[0] not in sdps:
                        sdps.append(splits[0])
                    if splits[1] not in samples:
                        samples.append(splits[1])
    else:
        pass

# get ptrs and cov of each of the sdps for all the samples

sample_cov = {x: [] for x in samples}
sample_ptr = {x: [] for x in samples}

for sample in samples:
    sample_cov[sample] = ["NA"] * len(sdps)
    sample_ptr[sample] = ["NA"] * len(sdps)

for file in input_files:
    if file_checked(file):
        with open(file, "r") as fh_in:
            for line in fh_in:
                splits = line.split()
                if splits[0] == "SDP":
                    pass
                else:
                    sdp = splits[0]
                    sample_name = splits[1]
                    cov = splits[2]
                    ptr = splits[3]
                    index = sdps.index(sdp)
                    sample_cov[sample_name][index] = cov
                    sample_ptr[sample_name][index] = ptr
# write to files
with open(cov_out, "w") as fh_cov:
    string = "SDP"
    for sample in samples:
        string = string + "," + sample
    fh_cov.write(f"{string}\n")
    for sdp in sdps:
        string = sdp
        index = sdps.index(sdp)
        for sample in samples:
            string = string + "," + sample_cov[sample][index]
        fh_cov.write(f"{string}\n")

with open(ptr_out, "w") as fh_ptr:
    string = "SDP"
    for sample in samples:
        string = string + "," + sample
    fh_ptr.write(f"{string}\n")
    for sdp in sdps:
        string = sdp
        index = sdps.index(sdp)
        for sample in samples:
            string = string + "," + sample_ptr[sample][index]
        fh_ptr.write(f"{string}\n"))
