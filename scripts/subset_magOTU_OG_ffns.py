#!/usr/bin/env python3
import sys
import os
import argparse
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def parse_prokka_get_MAG_name(header):
    """
    read header with prokka prefix and contig number
    and get MAG name eg. MAG_C1.1_12 from 'gnl|Prokka|MAG_C1.1_12_1'
    it removes the final number after the underscore which goes from 1 to n
    where n is the number of contigs in that MAG numbered by prokka
    """
    if "|" in header:
        parsed_header = "_".join(header.split("|")[-1].split("_")[:-1])
    else:
        parsed_header = header
    return(parsed_header)

def parse_MAG_name_from_gene_header(header):
    """
    read header with gene id in ffn files and get MAG name
    eg. MAG_C1.5_18 from '>MAG_C1.5_18_01186 hypothetical protein'
    it removes the final number after the underscore which goes from 1 to n
    where n is the gene id numbered by prokka
    """
    if ">" in header:
        parsed_header = "_".join(header.split(" ")[0].split(">")[1].split("_")[:-1])
    else:
        parsed_header = "_".join(header.split(" ")[0].split("_")[:-1])
    return(parsed_header)

def get_magOTUs_of_group(group_provided, ref_info):
    magOTU_list = []
    with open(ref_info, "r") as ref_info_fh:
        for line in ref_info_fh:
            line = line.strip()
            if line.startswith("ID"):
                header = line+"\n"
                continue
            genome_id = line.split("\t")[0]
            group_read = line.split("\t")[18]
            magOTU = line.split("\t")[11]
            if group_read == group_provided:
                if not magOTU in magOTU_list:
                    magOTU_list.append(magOTU)
            else:
                continue
    return(magOTU_list)

def get_mags_of_magOTU(magOTU_provided, ref_info):
    mags_list = []
    with open(ref_info, "r") as ref_info_fh:
        for line in ref_info_fh:
            line = line.strip()
            if line.startswith("ID"):
                header = line+"\n"
                continue
            genome_id = line.split("\t")[0]
            magOTU_read = line.split("\t")[11]
            if magOTU_read == magOTU_provided:
                if not genome_id in mags_list:
                    mags_list.append(genome_id)
            else:
                continue
    return(mags_list)

def get_seqs_from_ffn(ffn_file_path):
    seq_record_dict = dict()
    for seq_record in SeqIO.parse(ffn_file_path, "fasta"):
        gene_header = seq_record.id
        genome_id = parse_MAG_name_from_gene_header(gene_header)
        seq_record_dict[genome_id]=seq_record
    return(seq_record_dict)


# def is_valid_magOTU_OG(seq_records):
#     """
#     OG is valid if it is present in at least 1 mag for each magOTU
#     """
#     for mag in get()

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--group', metavar="group", required=True, help="Name of genus or group that the script is being run for (historically, phylotype)", action="store")
requiredNamed.add_argument('--ref_info',metavar="ref_info",required=True, help="Info about MAG status with MAG name in first column and other information accroding to checkpoint output", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
requiredNamed.add_argument('--input_seq_dir', metavar="input_seq_dir", required=True, help="Directory containing genes, aligned and back-translated genes", action="store")
requiredNamed.add_argument('--magOTU_seqs_dir_path', metavar="magOTU_seqs_dir_path", required=True, help="Path to store magOTU-wise seperated sequences (inside which magOTU_seqs_dir(s) exist)", action="store")
requiredNamed.add_argument('--log_path', metavar="log_path", required=True, help="Path to store magOTU-wise seperated sequences (inside which magOTU_seqs_dir(s) exist)", action="store")


args = parser.parse_args()
group = args.group
ref_info = args.ref_info
input_seq_dir = args.input_seq_dir
magOTU_seqs_dir_path = args.magOTU_seqs_dir_path
log_path = args.log_path

log_fh = open(log_path, "w")

magOTU_list = get_magOTUs_of_group(group, ref_info)

magOTU_seqs_dir_dict = {}
for magOTU in magOTU_list:
    magOTU_seqs_dir = os.path.join(magOTU_seqs_dir_path, magOTU)
    magOTU_seqs_dir_dict[magOTU] = magOTU_seqs_dir
    if os.path.exists(magOTU_seqs_dir):
        pass
    else:
        os.makedirs(magOTU_seqs_dir)

ffn_suffix = '*ffn'
ffn_files = glob.glob(os.path.join(input_seq_dir, ffn_suffix))
print(f"For {group}, {len(ffn_files)} OG groups were found")

if len(ffn_files) == 0:
        print(f"No fasta-files with the suffix '*ffn' were found in {input_seq_dir}. Check the directory:\n")
        sys.exit()

count_file_progress = 0
for ffn_file in ffn_files:
    seq_records = get_seqs_from_ffn(ffn_file)
    count_file_progress += 1
    for magOTU in magOTU_list:
        mags = get_mags_of_magOTU(magOTU, ref_info)
        for mag in mags:
            if mag not in seq_records.keys():
                log_fh.write(f"{mag} from {magOTU} not in {os.path.basename(ffn_file)}\n")
        magOTU_seq_records = {key: seq_records[key] for key in mags if key in seq_records.keys()}
        og_ffn = os.path.basename(ffn_file)
        magOTU_ffn_outfile =  os.path.join(magOTU_seqs_dir_dict[magOTU], og_ffn)
        if len(magOTU_seq_records) > 0:
            SeqIO.write(magOTU_seq_records.values(), magOTU_ffn_outfile, "fasta")
        else:
            print(f"skipping {og_ffn} for magOTU {magOTU} as none of its mags have it")
    if (count_file_progress % 100 == 0):
        print(f"Done sub-setting {count_file_progress} ffn-files..")
    log_fh.write(f"{ffn_file} has {len(seq_records)} MAGs or single-copy genes\n")

log_fh.close()
