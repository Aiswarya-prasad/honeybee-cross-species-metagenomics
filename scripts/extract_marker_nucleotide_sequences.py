from os import makedirs, system, getcwd, path
import numpy as np
import pandas as pd
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


'''
To calculate the substitution rate we use the bac120 marker gene nucleotide sequences
collected and aligned using a codon-aware aligner kcalign (https://github.com/davebx/kc-align)
First get the names of the marker genes that need to be collected and get the name of the gene
in the respective MAG corresponding to the marker gene and then extract the nucleotide sequence
from the MAG fasta file to write into the fasta file corresponding to the marker gene
'''

# Get the names of the marker genes
marker_genes = [path.basename(x).split('.fa')[0] for x in glob('results/09_MAGs_collection/gtdb_output/identify/intermediate_results/single_copy_fasta/bac120/*.fa')]
mag_gene_files = {path.basename(x).split('_protein.fna')[0]: x for x in glob('results/09_MAGs_collection/gtdb_output/identify/intermediate_results/marker_genes/*/*_protein.fna')}

makedirs('results/11_phylogenies/05_MAG_bac120_nucleotide_trees', exist_ok=True)
makedirs('results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences', exist_ok=True)

#  in retrospect, it would have been better to nest the marker loop inside the mag loop
for j, marker in enumerate(marker_genes):
    # print(f'-----working on {marker}-----\n')
    for i, mag in enumerate(mag_gene_files):
        print(f'         Marker: {j}/120 --- MAG {i}/{len(mag_gene_files)}        ', end='\r')
        # Get the name of the gene in the respective MAG corresponding to the marker gene
        # from pfam and tigrfam _tophit.tsv files
        gene_name = ''
        if 'PF' in marker:
            pfam_file = f'results/09_MAGs_collection/gtdb_output/identify/intermediate_results/marker_genes/{mag}/{mag}_pfam_tophit.tsv'
            pfam_df = pd.read_csv(pfam_file, sep='\t', header=0)
            if not pfam_df.empty:
                # split top hits column by , and the first is marker
                pfam_df['marker'] = pfam_df['Top hits (Family id,e-value,bitscore)'].str.split(',').str[0]
                if marker in pfam_df['marker'].values:
                    # get the marker gene name
                    gene_name = pfam_df.loc[pfam_df['marker'] == marker]['Gene Id'].values[0]
                else:
                    continue
            else:
                continue
        if 'TIGR' in marker:
            tigrfam_file = f'results/09_MAGs_collection/gtdb_output/identify/intermediate_results/marker_genes/{mag}/{mag}_tigrfam_tophit.tsv'
            tigrfam_df = pd.read_csv(tigrfam_file, sep='\t', header=0)
            if not tigrfam_df.empty:
                # split top hits column by , and the first is marker
                tigrfam_df['marker'] = tigrfam_df['Top hits (Family id,e-value,bitscore)'].str.split(',').str[0]
                if marker in tigrfam_df['marker'].values:
                    # get the marker gene name
                    gene_name = tigrfam_df.loc[tigrfam_df['marker'] == marker]['Gene Id'].values[0]
                else:
                    continue
            else:
                continue
        # Get the nucleotide sequence from the MAG fasta file
        mag_fasta_file = mag_gene_files[mag]
        mag_fasta = SeqIO.parse(mag_fasta_file, 'fasta')
        for record in mag_fasta:
            if gene_name  == record.id:
                # Write the nucleotide sequence to the marker gene fasta file
                with open(f'results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences/{marker}.fa', 'a') as f:
                    success = f.write(f'>{mag}\n{record.seq}\n')
                break