#!/usr/bin/env python3
import sys
import re
import argparse

# from KE's github repo https://github.com/kirsten2/Community_profiling/blob/master/bin/filter_bam.py
#samtools view -h Ig13619_sub.bam | python3 filter_bam.py -e 5 -m 50 | samtools sort - -o test.bam

#Usage: samtools view -h Ig13619_sub.bam | python3 filter_bam.py -e 5 -m 50 |

"""
This is a Python script that filters SAM-formatted sequence alignment data by read match length and maximum edit distance. 
It reads in SAM data from standard input, filters the data based on 
filtering parameters (minimum read match length and maximum edit distance), 
and writes the filtered SAM data to standard output. 
The script uses the argparse module to parse command-line arguments and the re module 
to extract matches from the CIGAR string of the SAM data.
"""

parser = argparse.ArgumentParser()
parser.add_argument('-m',required=False, help="Minimum read match length")
parser.add_argument('-e', required=False, help="Maximum edit-distance")
args = vars(parser.parse_args())

#Parse filtering thresholds (check that at least one is provided)
aln_threshold = "NA"
edit_dist = "NA"
if (args['m']):
    aln_threshold = args['m']
if (args['e']):
    edit_dist = args['e']
if aln_threshold == "NA" and edit_dist == "NA":
    print("No filtering parameters were provided, you must provide at least one:")
    parser.print_help(sys.stderr)
    sys.exit(1)

for line in sys.stdin:
    line = line.strip()
    if line.startswith('@'):
        print(line)
    else:
        split_line = line.split("\t")
        genome_match = split_line[2]
        if genome_match == '*': continue #Unmapped read, skip
        edit_dist_status = 0 #Change status to 1 if user-defined threshold is not met
        aln_length_status = 0 #Change status to 1 if user-defined threshold is not met
        if edit_dist != "NA":
            edit_dist_tab = [i for i in split_line if "NM:i" in i]
            split_tab = edit_dist_tab[0].split(':')
            edit_dist_line = split_tab[-1]
            if int(edit_dist_line) >= int(edit_dist):
                edit_dist_status = 1
        if aln_threshold != "NA":
            cigar = split_line[5]
            aln_length = 0
            matches = re.findall("[0-9]+M", cigar)
            if len(matches) == 0:
                aln_length_status = 1
            else:
                for match in matches:
                    match_length = match[0:-1]
                    aln_length += int(match_length)
            if aln_length <= int(aln_threshold):
                    aln_length_status = 1
        if edit_dist_status == 0 and aln_length_status == 0:
            print(line)
