import os
import glob
import pandas as pd
from Bio import SeqIO

'''
make a list of every scaffold and the bin they belong to for each sample
add unbinned if they are not binned mention the completeness and contamination of the bin
along with gtdb taxonomy and drep cluster if available
'''

outfile = 'results/09_MAGs_collection/all_scaffold_to_bin.tsv'
metadata_df = pd.read_csv('results/09_MAGs_collection/MAGs_metadata_summary.tsv', sep='\t')
metadata_df.columns
mags = glob.glob('results/09_MAGs_collection/MAGs/*.fa')

with open(outfile, 'w') as out_fh:
    out_fh.write('scaffold\tmag\tsize\tcompleteness\tcontamination\tgenus\tspecies\tmagotu\n')

for mag_file in mags:
    print(f'working on {mag_file}')
    with open(mag_file, 'r') as mag_fh:
        with open(outfile, 'a') as out_fh:
            mag = os.path.basename(mag_file).split('.fa')[0]
            size = metadata_df.loc[metadata_df['ID'] == os.path.basename(mag), 'Size'].values[0]
            completeness = metadata_df.loc[metadata_df['ID'] == os.path.basename(mag), 'Completeness'].values[0]
            contamination = metadata_df.loc[metadata_df['ID'] == os.path.basename(mag), 'Contamination'].values[0]
            genus = metadata_df.loc[metadata_df['ID'] == os.path.basename(mag), 'Genus'].values[0]
            species = metadata_df.loc[metadata_df['ID'] == os.path.basename(mag), 'Species'].values[0]
            magotu = metadata_df.loc[metadata_df['ID'] == os.path.basename(mag), 'magOTU'].values[0]
            mag_records = SeqIO.parse(mag_file, 'fasta')
            for record in mag_records:
                out_fh.write(f'{record.id}\t{mag}\t{size}\t{completeness}\t{contamination}\t{genus}\t{species}\t{magotu}\n')