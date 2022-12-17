#!/usr/bin/env python3

"""
Only argument is the path to the csv (must be xxxx.csv) file to be converted
Converts it into a tsv file with the same name but different extension

USAGE: python3 path/to/script/csv_to_tsv.py path/to/file.csv

expected output - path/to/file.tsv

!!! rewrite using argparse instead of sys!!!
"""

import sys
import csv

file_path = sys.argv[1]
out_path = file_path.split(".csv")[0] + ".tsv"

with open(out_path, "w+") as out_fh:
    with open(file_path, "r") as fh_in:
        for line_list in csv.reader(fh_in):
            line_list = [item.strip() for item in line_list]
            line_list = [item for item in line_list if item != "\n"]
            line = "\t".join(line_list)
            out_fh.write(f"{line}\n")
