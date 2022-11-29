#!/usr/bin/env python3
import sys
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import numpy as np

#Usage: python3 trim_aln.py OG0000319_aln_nuc.fasta

def trimmed_col_indices(aln):
    nb_seqs = len(aln)
    nb_col = aln.get_alignment_length()
    i = 0
    gap_char = '-'
    selected_columns = list()
    while (i < nb_col):
        aln_col = aln[:,i]
        nb_gaps = aln_col.count(gap_char)
        gap_fraction = round((nb_gaps/nb_seqs),3)
        if (gap_fraction < 0.5):
            selected_columns.append(i)
        i += 1
    return(selected_columns)

def trim_aln_np(aln, col_ids):
    aln_np = np.array([list(rec) for rec in aln], dtype=str)
    nb_seqs = len(aln)
    aln_np_trim = np.empty((nb_seqs,0),dtype=str)
    for i in col_ids:
        col_seq = aln_np[:,i].tolist()
        aln_np_trim = np.append(aln_np_trim,np.array([col_seq]).transpose(),axis=1)
    return(aln_np_trim)

def aln_np_to_bio(aln_np,seq_ids):
    nb_seqs = len(seq_ids)
    i = 0
    seq_records = list()
    while(i < nb_seqs):
        row_seq = aln_np[i,:].tolist()
        row_seq_str = ''.join(row_seq)
        row_seq_id = seq_ids[i]
        new_seq_record = SeqRecord(Seq(row_seq_str), id = row_seq_id, description="")
        seq_records.append(new_seq_record)
        i += 1
    new_aln = MultipleSeqAlignment(seq_records)
    return(new_aln)

#Open alignment file, read-in as alignment object with BioPython, store seq-ids in list 
try:
    aln_in = sys.argv[1]
except:
    print('Provide an input alignment file in fasta format. Exiting script')
    exit()
aln = AlignIO.read(aln_in, "fasta")
seq_ids = list()
for record in aln:
    seq_ids.append(record.id)
    
#Filter alignment columns, and generate trimmed alignment object
trim_columns = trimmed_col_indices(aln)
trimmed_aln_np = trim_aln_np(aln, trim_columns)
trimmed_aln = aln_np_to_bio(trimmed_aln_np,seq_ids)

#Print trimmed alignment to file
infile_name_split = aln_in.split('_')
outfile = infile_name_split[0] + '_aln_trim.fasta'
AlignIO.write(trimmed_aln,outfile,"fasta")
