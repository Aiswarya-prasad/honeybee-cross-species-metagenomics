import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
for a given gene, find the scaffold it is on 
using the gene_id and a given gff or bed file
and optionally the start and end positions and print them to stdout
and depending on what is asked for, info about the contig

example usage:
python scripts/find_gene_info.py --gene_id "gnl|Prokka|NC_000913.3_1" --type gff --file file.gff
"""

parser = argparse.ArgumentParser()
parser.add_argument('--gene_id', action='store', help='id of gene to find', required=True)
parser.add_argument('--type', action='store', help='type of file to search (gff or bed)')
parser.add_argument('--file', action='store', help='file to search')
parser.add_argument('--get', action='store', help='scaffold (default) or position or whokaryote')
args = parser.parse_args()
gene_id = args.gene_id
if args.type:
    type = args.type
if args.file:
    file = args.file
if args.get:
    get = args.get
else:
    get = "scaffold"

# gene_id = "C1-1_NODE_1_length_971086_cov_44.956229_1"
scaffold = ""
start = ""
end = ""
whokaryote = ""
if "NODE" in gene_id:
    scaffold = "_".join(gene_id.split("_")[:-1])
    sample = gene_id.split("_")[0]
    with open(f"results/06_metagenomicORFs/{sample}/{sample}.ffn", "r") as ffn_fh:
        for line in ffn_fh:
            if line.startswith(">"):
                if "_".join(gene_id.split("_")[1:]) in line:
                    start = line.split("# ")[1]
                    end = line.split("# ")[2]
                    break
            else:
                continue
    if get == "whokaryote":
        with open(f"results/05_assembly/contig_fates/whokaryote/{sample}/whokaryote_predictions_S.tsv", "r") as whokaryote_fh:
            for line in whokaryote_fh:
                line = line.strip()
                if line.split("\t")[0] in scaffold:
                    whokaryote = line.split("\t")[1]
                    break

if get == "scaffold":
    print(f"{scaffold}")
elif get == "position":
    print(f"{start}\t{end}")
elif get == "whokaryote":
    print(f"{whokaryote}")
# print(f"{scaffold}\t{start}\t{end}")