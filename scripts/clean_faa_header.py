import os
import sys
import argparse

"""
reads a fasta file and print just the headers into standard output
It removes everything after the first space (including hash characters)
in the header unless the --asis option is used 
The rest of the lines (the sequences) are unmodified
if option -i is not given, it reads from standard input

usage: python3 scripts/get_clean_gene_headers.py -i {input.faa} > {output.unfilt_faa}
"""

def main():
    parser = argparse.ArgumentParser(description='reads a fasta file and print just the headers into standard output')
    parser.add_argument('-i', metavar='fasta_file', 
                        help='fasta file')
    parser.add_argument('--asis', action='store_true', default=False,
                        help='print the header as is, without removing anything')
    args = parser.parse_args()
    fasta_file = None
    if args.i:
        fasta_file = args.i
    asis = False
    if args.asis:
        asis = args.asis

    if fasta_file is not None:
        if not os.path.isfile(fasta_file):
            print("Error: file not found: " + fasta_file)
            sys.exit(1)

        with open(fasta_file) as f:
            for line in f:
                if line.startswith('>'):
                    if asis:
                        print(line.strip())
                    else:
                        print(line.strip().split()[0])
                else:
                    print(line.strip())
    else:
        for line in sys.stdin:
            if line.startswith('>'):
                if asis:
                    print(line.strip())
                else:
                    print(line.strip().split()[0])
            else:
                print(line.strip())
    
if __name__ == "__main__":
    main()