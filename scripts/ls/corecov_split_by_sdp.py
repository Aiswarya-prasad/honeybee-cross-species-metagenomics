#!usr/bin/env python3

import os

in_coord_list = snakemake.input["coords"]
in_cov_list = snakemake.input["corecovs"]

for coord in in_coord_list:
    sdps = []
    with open(os.path.join(os.getcwd(), coord), "r") as coord_fh:
        phylotype = coord.split("/")[-1].split("_")[0]
        for line in coord_fh:
            if line.startswith("SDP"):
                continue
            else:
                if line.split("\t")[0] in sdps:
                    pass
                else:
                    sdps.append(line.split("\t")[0])
    for sdp in sdps:
        with open(os.path.join(os.getcwd(), coord), "r") as coord_fh:
            outfile = coord.split(phylotype)[0]+sdp+"_split"+coord.split(phylotype)[1]
            with open(outfile, "w+") as out_fh:
                for line in coord_fh:
                    if line.startswith("SDP"):
                        out_fh.write(line)
                    else:
                        if sdp in line:
                            out_fh.write(line)
                        else:
                            pass

for cov in in_cov_list:
    sdps = []
    with open(os.path.join(os.getcwd(), cov), "r") as cov_fh:
        phylotype = cov.split("/")[-1].split("_")[0]
        for line in cov_fh:
            if line.startswith("SDP"):
                continue
            else:
                if line.split("\t")[0] in sdps:
                    pass
                else:
                    sdps.append(line.split("\t")[0])
    for sdp in sdps:
        with open(os.path.join(os.getcwd(), cov), "r") as cov_fh:
            outfile = cov.split(phylotype)[0]+sdp+"_split"+cov.split(phylotype)[1]
            with open(outfile, "w+") as out_fh:
                for line in cov_fh:
                    if line.startswith("SDP"):
                        out_fh.write(line)
                    else:
                        if sdp in line:
                            out_fh.write(line)
                        else:
                            pass
