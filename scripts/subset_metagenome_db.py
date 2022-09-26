#!/usr/bin/env python3

# assume that it is run from the Snakefile working directory...
# ------ what does best practice dictate here?

# Not part of Snakefile
# used to make updates when needed
# find a way to incorporate in pipeline only if needed
# or identify when needed and print a message saying, use this!

# Orthofile needs to be made from faa using orthofinder - implement this later
# if needed
# implement this and rename this script as update orthofile or something similar
# it will also need to make
# genome_db meta file - ask KE how? or make manually
# bed_files - if not already made for orthofile that comes with genome_db_210402
# DO everything from chosen "faa?" files


# for now, just subset from files that came with genome_db_210402 set up by KE's
# download.py as I am now only concerned about the reduced database for SNVs

import os
from Bio import SeqIO
# import argparse

def check_file(file):
    try:
        fh = open(file)
    except:
        print(f"Cant find/open this file: {file}")
        print("Exiting script!")
        exit()

def run_orthofinder(faa_dir):
    print("not implemented")

def make_metafile(glist):
    with open(os.path.join(os.getcwd(), snakemake.input["metafile"]), "r") as fh_metafile:
        with open(os.path.join(os.getcwd(), snakemake.output["metafile"]), "w") as fh_outfile:
            for line in fh_metafile:
                split_line = line.split("\t")
                genome_id = split_line[0]
                if genome_id in glist:
                    # did not do line strip
                    # so it includes \n
                    fh_outfile.write(line)
                else:
                    pass

def make_db_red(glist):
    input_fasta = os.path.join(os.getcwd(), snakemake.input["db"])
    output_fasta = os.path.join(os.getcwd(), snakemake.output["db"])
    records = (r for r in SeqIO.parse(input_fasta, "fasta") if r.id in glist)
    number_red = SeqIO.write(records, output_fasta, "fasta")
    print(f"adding {number_red} records to reduced database")

def subset_orthofile(glist):
    for orthofile in os.listdir(snakemake.input["orthodir"]):
        if orthofile.endswith("_single_ortho_filt.txt"):
            with open(os.path.join(os.getcwd(), snakemake.input["orthodir"], orthofile), "r") as fh_orthofile:
                out_dir = os.path.join(os.getcwd(), snakemake.output["orthodir"])
                if not os.path.isdir(out_dir):
                    os.makedirs(out_dir)
                output_file = os.path.join(out_dir, orthofile)
                print(f"Writing to {output_file}")
                with open(output_file, "w") as fh_orthofile_out:
                    for line in fh_orthofile:
                        line = line.strip()
                        split_line = line.split()
                        OG_id = split_line.pop(0)
                        out_string = OG_id
                        for gene in split_line:
                            genome = gene.split("_")[0]
                            if genome in glist:
                                # add it to output
                                out_string = out_string + " " + gene
                            else:
                                pass
                        out_string = out_string + "\n"
                        # add this line to output
                        fh_orthofile_out.write(out_string)
                    # move to next orthofile



list_genomes = list()
with open(os.path.join(os.getcwd(), snakemake.input["metafile"]), "r") as fh_metafile:
    for line in fh_metafile:
        if line.split("\t")[3] == "Ref":
            list_genomes.append(line.split("\t")[0])

print(f"making db_red {snakemake.output["db_red"]}")
make_db_red(list_genomes)

print(f"making metafile {snakemake.output["metafile"]}")
make_metafile(list_genomes)

print(f"making orthodir {snakemake.output["orthodir"]}")
subset_orthofile(list_genomes)
