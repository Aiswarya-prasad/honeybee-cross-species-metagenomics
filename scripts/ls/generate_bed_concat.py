#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Usage: python generate_bed_concat.py Locustag_Phylotype_SDP_final.txt

#Function for generate concat-file of assembly contigs. Will also return an ordered list of contigs, and their sizes
def make_concat_fna (fna_file, genome):
    seq_strings = list()
    contig_list = list()
    contig_sizes = dict()
    for seq_record in SeqIO.parse(fna_file, "fasta"):
        contig_seq = str(seq_record.seq)
        seq_strings.append(contig_seq)
        contig_list.append(seq_record.id)
        contig_length = len(seq_record.seq)
        contig_sizes[seq_record.id] = contig_length
    concat_seq = "".join(seq_strings)
    concat_filename = genome + '_concat.fna'
    concat_seq_record = SeqRecord(Seq(concat_seq), id=genome, description="")
    SeqIO.write(concat_seq_record, concat_filename, "fasta")
    return(contig_list, contig_sizes)

#Function for extracting gff-lines corresponding to CDS from gff-file, and sorting them by their associated contig. Will return dictionary of contigs, with a list of gff-lines for each contig
def get_contig_gff(gff_file):
    fh_gff_file = open(gff_file)
    contig_gff_lines = dict()
    for line in fh_gff_file:
        line = line.strip()
        if line.startswith("#"): continue
        split_line = line.split("\t")
        contig_id = split_line[0]
        if contig_id not in contig_gff_lines:
            contig_gff_lines[contig_id] = list()
        if split_line[2] == "CDS":
            contig_gff_lines[contig_id].append(line)
    return(contig_gff_lines)

#Function for re-calculating CDS positions relative to concatenated assembly-file. Will return a list of new positions for each CDS, in bed-format. Will also check that the CDS is contained within the faa-file.
def get_new_cds_pos(contig_gff_lines, contig_sizes, filt_tags):
    pos = 0
    contig_bed_lines = list()
    for contig in contig_gff_lines.keys():
        if (len(contig_gff_lines[contig]) > 0):
            gff_lines = contig_gff_lines[contig]
            for line in gff_lines:
                split_line = line.split("\t")
                desc_tab = split_line[8]
                split_desc = desc_tab.split(';')
                tag_list = [i for i in split_desc if "locus_tag=" in i]
                split_tag = tag_list[0].split('=')
                locustag = split_tag[1]
                if locustag in filt_tags:
                    start_pos = split_line[3]
                    end_pos = split_line[4]
                    new_start_pos = pos + int(start_pos)
                    new_end_pos = pos + int(end_pos)
                    bed_line = [genome, str(new_start_pos), str(new_end_pos), locustag]
                    bed_line_str = "\t".join(bed_line)
                    contig_bed_lines.append(bed_line_str)
        pos += contig_sizes[contig]
    return(contig_bed_lines)

#Check that the requisite directories with genome-files are present in the run dir. Prepare the bed-file dir
cwd = os.getcwd()
faa_file_dir = 'faa_files'
fna_file_dir = 'fna_files'
gff_file_dir = 'gff_files'
bed_file_dir = 'bed_files'
genome_dirs = [faa_file_dir, fna_file_dir, gff_file_dir]
for dir in genome_dirs:
    if os.path.isdir(dir): continue
    else:
        print("The following directory doesnt exist in the run-dir:", dir)
        print("Exiting script")
        exit()
if not (os.path.isdir(bed_file_dir)):
    os.mkdir(bed_file_dir)

#Read the meta-file containing locustag identifiers for all genomes in db
fh_metafile = open(sys.argv[1])
genomes = dict()
for line in fh_metafile:
    line = line.strip()
    split_line = line.split('\t')
    locustag = split_line[0]
    genomes[locustag] = 1
fh_metafile.close()

#Process the genome-files
for genome in genomes.keys():
    print("Processing:", genome)
    faa_file = genome + '.faa'
    fna_file = genome + '.fna'
    gff_file = genome + '.gff'
    #Go to faa-dir, and get the locus-tags for the protein-coding genes
    genome_locustags = dict()
    os.chdir(faa_file_dir)
    for seq_record in SeqIO.parse(faa_file, "fasta"):
        locustag = seq_record.id
        genome_locustags[locustag] = 1
    os.chdir(cwd)
    #Go to fna-dir, and get contig-ids in order, plus their sizes
    os.chdir(fna_file_dir)
    genome_contig_list, genome_contig_sizes = make_concat_fna(fna_file, genome)
    os.chdir(cwd)
    #Go to gff-dir, and get new positions for CDS contained within faa-file (generate bed-lines)
    os.chdir(gff_file_dir)
    genome_gff_lines = get_contig_gff(gff_file)
    genome_bed_lines = get_new_cds_pos(genome_gff_lines, genome_contig_sizes, genome_locustags)
    os.chdir(cwd)
    #Print bed-file to bed-file dir
    os.chdir(bed_file_dir)
    bed_file = genome + '.bed'
    fh_bed_file = open(bed_file, 'w')
    for line in genome_bed_lines:
        fh_bed_file.write(line + "\n")
    fh_bed_file.close()
    os.chdir(cwd)
