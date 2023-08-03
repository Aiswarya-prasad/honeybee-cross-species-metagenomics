import os
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Accept a file containing sequences IDs to filter from a GFF file.
The input file is assumed to be in fasta format with the headers
being the IDs to filter.
usage: python3 scripts/get_filt_gff.py --ffn_in {input.filt_ffn} --gff_in {input.gff} --gff_out {output.filt_gff}
"""

def get_args():
    parser = argparse.ArgumentParser(description="Filter a GFF file based on a list of sequence IDs")
    parser.add_argument("--ffn_in", help="Input fasta file containing sequence IDs to filter")
    parser.add_argument("--gff_in", help="Input GFF file to filter")
    parser.add_argument("--gff_out", help="Output GFF file")
    return parser.parse_args()
ffn_in = "results/06_metagenomicORFs/D9-1/filt_orfs/D9-1.ffn"
gff_in = "results/06_metagenomicORFs/D9-1/prodigal_out/D9-1.gff"
def get_filt_ids(ffn_in):
    """
    Get IDs to filter from a fasta file using bioseq
    """
    filt_ids = []
    for record in SeqIO.parse(ffn_in, "fasta"):
        filt_ids.append(record.id)
    return filt_ids

def get_filt_gff(gff_in, filt_ids, gff_out):
    """
    Filter a GFF file based on a list of sequence IDs
    """
    with open(gff_in, "r") as f:
        with open(gff_out, "w") as out_fh:
            for line in f:
                if line.startswith("#"):
                    out_fh.write(line)
                else:
                    line_header = line.split("\t")[0] + "_" + \
                            line.split("\t")[8].split("ID=")[1].split(";")[0].split("_")[1]
                    if line_header in filt_ids:
                        out_fh.write(line)

def main():
    args = get_args()
    filt_ids = get_filt_ids(args.ffn_in)
    get_filt_gff(args.gff_in, filt_ids, args.gff_out)

if __name__ == "__main__":
    main()
