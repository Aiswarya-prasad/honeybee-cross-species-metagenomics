#!/usr/bin/env python3
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO

#Usage: python3 aln_aa_to_dna.py OG0000319_aln.fasta OG0000319.ffn

def get_seq_dict(seq_obj):
    seq_dict = dict()
    for record in seq_obj:
        seq_id = record.id
        seq = record.seq
        seq_dict[seq_id] = seq
    return(seq_dict)

def aa_to_nuc(aln_seq, nuc_seq):
    nuc_aln = str()
    for aa in aln_seq:
        if (aa == '-'):
            nuc_aln = nuc_aln + '---'
        else:
            nuc_codon = nuc_seq[0:3]
            nuc_aln = nuc_aln + nuc_codon
            nuc_seq = nuc_seq[3:]
    nuc_aln_str = str(nuc_aln)
    return(nuc_aln_str)

#Read the input-files into objects, generate dictionaries (seq_id - seq) for both files
try:
    aln_in = sys.argv[1]
    ffn_in = sys.argv[2]
except:
    print('Provide an input amino acid alignment file in fasta format, and a corresponding fasta-file with nucleotide sequences. Exiting script')
    exit()

aln_aa = AlignIO.read(aln_in, "fasta")
aln_aa_dict = get_seq_dict(aln_aa)
ffn = list(SeqIO.parse(ffn_in, "fasta"))
ffn_dict = get_seq_dict(ffn)

#Translate each aa-seq from alignment file to the corresponding nuc-seq, generate new seq-object, and save in list
aln_nuc_seq_records = list()  
for seq_id in aln_aa_dict.keys():
    aln_nuc = aa_to_nuc(aln_aa_dict[seq_id], ffn_dict[seq_id])
    new_record = SeqRecord(
        Seq(aln_nuc),
        id = seq_id,
        description = ""
    )
    aln_nuc_seq_records.append(new_record)

#Print the translated sequences to file
OG_id = sys.argv[1][0:9]
outfile_name = OG_id + '_aln_nuc.fasta'
SeqIO.write(aln_nuc_seq_records, outfile_name, "fasta")