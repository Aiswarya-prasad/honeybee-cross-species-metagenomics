import os
import argparse

"""
convert a given gff file to a bed file
"""

parser = argparse.ArgumentParser(description='convert a given gff file to a bed file if --bed is provided it will save it there if not it will print it to stdout. default gff id is scaffold#_gene# and is used as the name in the bed file.')
parser.add_argument('--gff', action='store', help='gff file as input')
parser.add_argument('--bed', action='store', help='faa file as input')
parser.add_argument('--rename', action='store_true', help='if provided the name column will be renamed as genome_xxxx where genomes is already appended to the default prodigal id. where xxxx is the four digit version of the index of the gene in the order that they are listed in the gff file else they will be named as they are in the gff file ID flag. DO NOT use this on unmodified prodigal output. To rename default prodigal output ID use a simple sed command: sed \'s/ID=/ID=prefix_/g\' prodigal.gff > prodigal_renamed.gff')
args = parser.parse_args()

gff = args.gff
bed = args.bed
rename = args.rename

def four_digit(n):
    n = int(n)
    if n < 10:
        return f"000{n}"
    elif n < 100:
        return f"00{n}"
    elif n < 1000:
        return f"0{n}"
    else:
        if n > 9999:
            print("Gene number seems too high. Check your genome!")
        return f"{n}"

if bed:
    bed_fh = open(bed, "w")
gene_iter = 1
with open(gff, "r") as gff_fh:
    # bed format expects no header line!
    # if bed:
    #     bed_fh.write("\t".join(["chrom","start","end","name"])+"\n")
    # else:
    #     print("\t".join(["chrom","start","end","name"]))
    for line in gff_fh:
        line = line.strip()
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            break
        atributes = {x.split("=")[0] : x.split("=")[1] for x in line.split("\t")[-1].split(";") if x}
        if "ID" not in atributes.keys():
            print(line)
            continue
        if rename:
            genome = "_".join(atributes["ID"].split("_")[:-2])
            gene_id = genome + "_" + four_digit(gene_iter)
            gene_iter += 1
        else:
            gene_id = atributes["ID"]
        scaffold_id = line.split("\t")[0]
        start = line.split("\t")[3]
        end = line.split("\t")[4]
        strand = line.split("\t")[6]
        gene_type = line.split("\t")[2]
        if gene_type != "CDS":
            continue
        if bed:
            bed_fh.write(f"{scaffold_id}\t{start}\t{end}\t{gene_id}\n")
        else:
            print(f"{scaffold_id}\t{start}\t{end}\t{gene_id}")

if bed:
    bed_fh.close()