import os
import sys
import argparse
import pandas as pd
import numpy as np

"""
python3 make_phylo_metadata.py \
            --mag_metadata {input.mag_metadata} \
            --Isolate_metadata {input.Isolate_metadata} \
            --phylo_metadata {output.phylo_metadata}
"""

parser = argparse.ArgumentParser(description='make a summary table of the MAGs and isolates to use to make phylogenetic trees')
parser.add_argument('--mag_metadata', help='mag metadata file')
parser.add_argument('--Isolate_metadata', help='isolate metadata file')
parser.add_argument('--phylo_metadata', help='output file name')

args = parser.parse_args()

mag_metadata = args.mag_metadata
# Isolate_metadata = "config/Outgroup_isolate_genomes.tsv"
Isolate_metadata = args.Isolate_metadata
phylo_metadata = args.phylo_metadata

# mag_metadata = "results/09_MAGs_collection/All_mags_sub/All_mags_sub_metadata_summary.tsv"
# Isolate_metadata = "config/Outgroup_isolate_genomes.tsv"

# read in the mag metadata file
mag_metadata = pd.read_csv(mag_metadata, sep='\t', header=0)
# only include genera wirth >2 mags of high or medium quality
mag_metadata = mag_metadata.groupby('Genus').filter(lambda x: len(x) > 2)
mag_metadata = mag_metadata[(mag_metadata['Quality'] == 'high') | (mag_metadata['Quality'] == 'medium')]
sample_host_dict = {'M' : 'Apis mellifera',
                    'C' : 'Apis cerana',
                    'D' : 'Apis dorsata',
                    'F' : 'Apis florea',
                    'A' : 'Apis andreniformis'
                    }
mag_metadata['Host'] = [sample_host_dict[x] for x in mag_metadata['ID'].str[0]]
mag_metadata['Phylogeny'] = mag_metadata['Genus']

# read isolate genome metadata file
Isolate_metadata_df = pd.read_csv(Isolate_metadata, sep='\t', header=0)

Isolate_metadata_df['Phylogeny'] = Isolate_metadata_df['Group']
# if strain name is not in species name append it with a space
Isolate_metadata_df['Species'] = [species if strain in species else species + ' ' + strain for strain, species in zip(Isolate_metadata_df['Strain'], Isolate_metadata_df['Species'])]

# make output dataframe

output_df = pd.DataFrame(columns=['ID', 'Species', 'magOTU', 'Type', 'Completeness', 'Host', 'Phylogeny_group'])
output_df['ID'] = mag_metadata['ID'].to_list() + Isolate_metadata_df['Name'].to_list()
output_df['Species'] = mag_metadata['Species'].to_list() + Isolate_metadata_df['Species'].to_list()
output_df['magOTU'] = mag_metadata['magOTU'].to_list() + ['NA'] * len(Isolate_metadata_df)
output_df['Type'] = ['MAG' for n in mag_metadata['magOTU']] + ['Isolate' for n in Isolate_metadata_df['Strain']]
output_df['Completeness'] = mag_metadata['Completeness'].to_list() + Isolate_metadata_df['Completeness'].to_list()
output_df['Host'] = mag_metadata['Host'].to_list() + Isolate_metadata_df['Host'].to_list()
output_df['Phylogeny_group'] = mag_metadata['Phylogeny'].to_list() + Isolate_metadata_df['Phylogeny'].to_list()


# write to output file
output_df.to_csv(phylo_metadata, sep='\t', header=True, index=False)

with open("results/11_phylogenies/assembly_summary.txt", 'r') as f:
    for line in f:
        if line.startswith('#') and 'ftp_path' in line:
            ind_link = line.strip().split('\t').index('ftp_path')
        else:
            line_split = line.strip().split('\t')

