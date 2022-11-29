#!/usr/bin/env python3
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse

#Usage: python3 extract_orthologs.py firm5_single_ortho.txt

def genome_count_line(seq_ids):
    genomes = dict()
    for seq_id in seq_ids:
        split_id = seq_id.split('_')
        genome_id = "_".join(split_id[:-1])
        genomes[genome_id] = genomes.get(genome_id,0) + 1
    return(genomes)

def get_seq_objects(filename):
    seq_record_dict = dict()
    for seq_record in SeqIO.parse(filename, "fasta"):
        if (seq_record.seq[-1] == '*'):  #Sometimes, stop-codons are indicated with an asterisk for amino acid gene sequences, which in turn generates warnings when using "muscle" for alignments. This block of code remove the trailing asterisk when present.
            seq_string = str(seq_record.seq)
            seq_chomp = seq_string[:-1]
            new_record =  SeqRecord(
                Seq(seq_chomp),
                id=seq_record.id,
                name="",
                description=""
                )
            seq_record_dict[seq_record.id] = new_record
        else:
            seq_record_dict[seq_record.id] = seq_record
    return(seq_record_dict)
# to do update annotation
parser = argparse.ArgumentParser()
parser.add_argument('--orthofile', action='store', help='*_ single_ortho.txt as input')
parser.add_argument('--outdir', action='store', help='name of directory to write outputs to')
mutually_exclusive_group = parser.add_mutually_exclusive_group()
mutually_exclusive_group.add_argument('--faaffndir', action='store', help='name of directory where annotation output directories are')
mutually_exclusive_group.add_argument('--faadir', action='store', help='name of directory where faa files aer')
parser.add_argument('--ffndir', action='store', help='name of directory where ffn files are')
args = parser.parse_args()

ortho_single = args.orthofile
if args.faaffndir:
    faaffndir = args.faaffndir
else:
    faaffndir = None
if args.faadir and args.ffndir:
    faadir = args.faadir
    ffndir = args.ffndir
else:
    faadir = None
    ffndir = None
ortho_seq_dir = args.outdir

#Open the ortholog file
if os.path.exists(ortho_single):
    pass
else:
    print('Input ortho-file not found. Exiting script')
    exit()


#Get all the genome-ids present in the ortholog-file, and all the gene-ids associated with each gene-family
print('reading single_ortho file')
genome_ids = dict()
OG_fams = dict()
with open(ortho_single, "r") as fh_ortho_in:
    for line in fh_ortho_in:
        line = line.strip()
        split_line = line.split()
        OG_id = split_line.pop(0)[:-1]
        OG_fams[OG_id] = split_line
        genome_count = genome_count_line(split_line)
        for genome in genome_count:
            genome_ids[genome] = 1

#Construct dictionaries of gene-sequences

ffn_seq_objects = dict()
faa_seq_objects = dict()
for genome in genome_ids.keys():
    if faaffndir:
        ffn_file = os.path.join(faaffndir, genome, genome + '.ffn')
        faa_file = os.path.join(faaffndir, genome, genome + '.faa')
    if faadir and ffndir:
        ffn_file = os.path.join(ffndir, genome + '.ffn')
        faa_file = os.path.join(faadir, genome + '.faa')
    if os.path.exists(ffn_file) and os.path.exists(faa_file):
        pass
    else:
        print('One or both of these files are missing: ')
        print(ffn_file)
        print(faa_file)
        print('Exiting script')
        exit()
    genome_ffn_seq_objects = get_seq_objects(ffn_file)
    ffn_seq_objects.update(genome_ffn_seq_objects)
    genome_faa_seq_objects = get_seq_objects(faa_file)
    faa_seq_objects.update(genome_faa_seq_objects)

#Create the output directory, and move into it
output_dir = ortho_seq_dir

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

#Print multi-fasta files corresponding to each gene-family in the output dir
for OG in OG_fams.keys():
    OG_seq_ids = OG_fams[OG]
    OG_ffn_seq_obj = [ffn_seq_objects[x] for x in OG_seq_ids]
    OG_faa_seq_obj = [faa_seq_objects[x] for x in OG_seq_ids]
    ffn_outfile = os.path.join(output_dir, OG + '.ffn')
    faa_outfile = os.path.join(output_dir, OG + '.faa')
    SeqIO.write(OG_ffn_seq_obj,ffn_outfile,"fasta")
    SeqIO.write(OG_faa_seq_obj,faa_outfile,"fasta")
