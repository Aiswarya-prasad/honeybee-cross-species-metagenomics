#!/usr/bin/env python3

import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l',metavar="length_threshold",required=True, help="Contigs shorter than the threshold will be removed", action="store")
parser.add_argument('-c', metavar='coverage_threshold', required=True, help="Contigs with coverage less than the threshold will be removed", action="store")
parser.add_argument('-i', metavar='input_contigs', required=True, help="File contiaining original unfiltered contigs", action="store")
parser.add_argument('-o', metavar='output_contigs', required=True, help="File contiaining filtered contigs (will be created and old files replaced)", action="store")

args = vars(parser.parse_args())

length_threshold = int(args['l'])
coverage_threshold = int(args['c'])
f_in = args['i']
f_out = args['o']

print("Filtering contigs from " + f_in + " to " + f_out)

def accept_contig(header):
    length = header.split("_")[3]
    cov = header.split("_")[5]
    if float(length) > length_threshold and float(cov) > coverage_threshold:
        return(True)
    else:
        return(False)

write_line = False
count = 0
count_filt = 0

with open(f_in, "r") as fh_in:
    with open(f_out, "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                count = count + 1
                if accept_contig(line):
                    write_line = True
                    count_filt = count_filt + 1
                else:
                    write_line = False
            if write_line:
                fh_out.write(line)
            else:
                pass
