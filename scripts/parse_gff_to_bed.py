import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
makes the genes file for midas as described here
https://github.com/snayfach/MIDAS/blob/master/docs/build_db.md
only tested on prokka (Prodigal:002006) gff files
"""

parser = argparse.ArgumentParser()
parser.add_argument('--gff', action='store', help='gff file as input')
parser.add_argument('--faa', action='store', help='faa file as input')
parser.add_argument('--outfile', action='store', help='Bed file for database')
args = parser.parse_args()
faa = args.faa
gff = args.gff
outfile = args.outfile

genes_present_in_faa = set()
for seq_record in SeqIO.parse(faa, "fasta"):
    locustag = seq_record.id
    genes_present_in_faa.add(locustag)

with open(gff, "r") as gff_fh:
    with open(outfile, "w") as outfile_fh:
        # outfile_fh.write("\t".join(["gene_id","scaffold_id","start","end","strand","gene_type"])+"\n")
        print("skipped gff entries:\n")
        for line in gff_fh:
            if line.startswith("##"):
                continue
            if line.startswith(">"):
                break
            info_dict = {x.split("=")[0] : x.split("=")[1] for x in line.split("\t")[-1].split(";")}
            if "ID" not in info_dict.keys():
                print(line)
                continue
            gene_id = info_dict["ID"]
            scaffold_id = line.split("\t")[0].split("|")[-1]
            start = line.split("\t")[3]
            end = line.split("\t")[4]
            strand = line.split("\t")[6]
            gene_type = line.split("\t")[2]
            if gene_type != "CDS":
                continue
            if gene_id in genes_present_in_faa:
                outfile_fh.write(f"gnl|Prokka|{scaffold_id}\t{start}\t{end}\t{gene_id}\n")
            else:
                print(f"{gene_id} not present in faa file.")
