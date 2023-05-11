#!/usr/bin/env python3
import sys
import os.path
from Bio import SeqIO

"""
summary of what the script does
"""
#Usage: python3 xxx.py ...

#Parse input options
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--xxx', metavar='xxx', required=True, help="xxx", action="store")
requiredNamed.add_argument('--xxx', metavar="xxx", required=True, help="xxx", action="store")
requiredNamed.add_argument('--xxx_list', nargs = "+", metavar="xxx",required=True, help="xxx", action="store")
args = parser.parse_args()

xxx = args.xxx
xxx = args.xxx
xxx = args.xxx