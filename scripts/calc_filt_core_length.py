#!/usr/bin/env python3
import os

# This script will read all core bedfiles  (*core_filt.bed), and create a small
# summary table with the total core gene length. Since genes may be overlapping,
# it is not simply a matter of summing the gene-lengths.
# The overlaps are subtracted with the current script.


#Sub-routine
# For line 1 you store start (let's call it tmp_start) and end (tmp_end). Initialize a value=0  variable.
# For line 2, look at start:
# if it is >=tmp_end , then . Then add the difference between tmp_endand tmp_start  to the value variable.
# if it is <tmp_end, then store end in tmp_end.
# Then repeat the line 2 process for the next lines until file is done.
def calc_length(bedfile):
    length = 0
    with open(bedfile, "r") as bf:
        line = bf.readline()
        p_start = int(line.split()[1])
        p_end = int(line.split()[2])
        for line in bf:
            c_start = int(line.split()[1])
            c_end = int(line.split()[2])
            if (c_start >= p_end):
                if c_start == p_end:
                    length += p_end - p_start
                else:
                    length += p_end - p_start + 1
                p_start = c_start
                p_end = c_end
            else:
                p_end = max(c_end, p_end)
        length += p_end - p_start + 1
        return(length)

bed_files = snakemake.input["core_red_beds"]

genomes_dict = snakemake.params["genomes_sdp_dict"]

outfile = os.path.join(os.getcwd(), snakemake.output["lengths"])

#Get the names of the files that are to be analyzed

with open(outfile, "w+") as out_fh:
    out_fh.write("sdp\tgenome\tcore_length\n")
    for sdp in genomes_dict.keys():
        genome = genomes_dict[sdp]
        file = [os.path.join(os.getcwd(), bed_file) for bed_file in bed_files if (genome in bed_file)][0]
        core_length = calc_length(file)
        out_fh.write(f"{sdp}\t{genome}\t{core_length}\n")
