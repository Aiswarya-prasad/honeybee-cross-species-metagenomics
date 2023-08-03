import os
import sys
import re
import glob
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objs as go
import plotly.express as px
import plotly.io as pio
import plotly.figure_factory as ff
import plotly.offline as offline
from plotly.subplots import make_subplots
from plotly.offline import plot, iplot
from plotly.graph_objs import *
from plotly import tools

"""
This script takes the output of whokaryote and kaiju and summarises the contig fates
in terms of a final table (for each sample) where each row is a contig and each column
summarises the following:
    - Contig ID (with sample )
    - Contig ID in sample (contig_id)
    - Contig length (contig_length)

Usage:
python scripts/summarize_contig_fates.py \
            --scaffolds {input.scaffolds} \
            --whokaryote_out {input.whokaryote_out} \
            --kaiju_out {input.kaiju_names} \
            --kaiju_names {input.kaiju_names} \
            --kaiju_names_full {input.kaiju_names_full} \
            --bins_directory {input.bins_directory} \
            --gtdb_bac {input.gtdb_bac} \
            --gtdb_ar {input.gtdb_ar} \
            --sample {wildcards.sample} \
            --outfile {output.contig_fates}
"""

def parse_args():
    parser = argparse.ArgumentParser(description='Summarise contig fates')
    parser.add_argument('--scaffolds', help='Scaffolds file containing all scaffolds for each sample')
    parser.add_argument('--whokaryote_out', help='whokaryote output file')
    parser.add_argument('--kaiju_out', help='kaiju output file')
    parser.add_argument('--kaiju_names', help='kaiju output with taxonomic names')
    parser.add_argument('--kaiju_names_full', help='kaiju output with full taxonomic names')
    parser.add_argument('--bins_directory', help='Bins directory containing all bins from the sample')
    parser.add_argument('--gtdb_bac', help='GTDB taxonomy for bacteria for MAGs')
    parser.add_argument('--gtdb_ar', help='GTDB taxonomy for archaea for MAGs')
    parser.add_argument('--sample', help='GTDB taxonomy for archaea for MAGs')
    parser.add_argument('--outfile', help='Output file')
    return parser.parse_args()

args = parse_args()

scaffolds = args.scaffolds
whokaryote_out = args.whokaryote_out
kaiju_out = args.kaiju_out
kaiju_names = args.kaiju_names
kaiju_names_full = args.kaiju_names_full
bins_directory = args.bins_directory
gtdb_bac = args.gtdb_bac
gtdb_ar = args.gtdb_ar
sample = args.sample
outfile = args.outfile

# sample = "A1-1"
# scaffolds = f"/scratch/aprasad/20230313_apis_species_comparison/results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta"
# whokaryote_out = f"/scratch/aprasad/20230313_apis_species_comparison/results/05_assembly/contig_fates/whokaryote/{sample}/whokaryote_predictions_S.tsv"
# kaiju_out = f"/scratch/aprasad/20230313_apis_species_comparison/results/05_assembly/contig_fates/kaiju/nr/{sample}.kaiju"
# kaiju_names = f"/scratch/aprasad/20230313_apis_species_comparison/results/05_assembly/contig_fates/kaiju/nr/{sample}_names.txt"
# kaiju_names_full = f"/scratch/aprasad/20230313_apis_species_comparison/results/05_assembly/contig_fates/kaiju/nr/{sample}_fullnames.txt"
# bins_directory = f"/scratch/aprasad/20230313_apis_species_comparison/results/07_MAG_binng_QC/02_bins/{sample}/"
# gtdb_bac = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/gtdb_output/classify/All_mags_sub.bac120.summary.tsv"
# gtdb_ar = "/scratch/aprasad/20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/gtdb_output/classify/All_mags_sub.ar53.summary.tsv"
# outfile = f"/scratch/aprasad/20230313_apis_species_comparison/results/05_assembly/contig_fates/{sample}/contig_fates_{sample}.tsv"

contig_ids = []
contig_ids_dict = {}
contig_length_dict = {}
contig_bins_dict = {}
contig_whokaryote_dict = {}
contig_kaiju_dict = {}
contig_kaiju_full_dict = {}
mag_gtdb_dict = {}
contig_gtdb_dict = {}

with open(scaffolds, 'r') as scaffolds_fh:
    for line in scaffolds_fh:
        if line.startswith('>'):
            line = line.strip()
            contig_id = line.split(' ')[0].replace('>', '')
            node_id = contig_id.split('NODE_')[1].split('_')[0]
            contig_ids.append(contig_id)
            contig_ids_dict[node_id] = contig_id
            contig_length_dict[node_id] = contig_id.split('_')[3]

with open(whokaryote_out, 'r') as whokaryote_out_fh:
    header = whokaryote_out_fh.readline()
    for line in whokaryote_out_fh:
        line = line.strip().split('\t')
        contig_id = line[0]
        node_id = contig_id.split('NODE_')[1].split('_')[0]
        contig_whokaryote = line[1]
        if contig_whokaryote == 'Unclassified':
            contig_whokaryote = 'whokaryote_unclassified'
        contig_whokaryote_dict[node_id] = contig_whokaryote

with open(kaiju_names, 'r') as kaiju_names_fh:
    for line in kaiju_names_fh:
        line = line.strip().split('\t')
        status = line[0]
        contig_id = line[1]
        node_id = contig_id.split('NODE_')[1].split('_')[0]
        if status == 'C':
            kaiju_tax = line[7]
            contig_kaiju_dict[node_id] = kaiju_tax
        else:
            contig_kaiju_dict[node_id] = 'Unassigned'

with open(kaiju_names_full, 'r') as kaiju_names_full_fh:
    for line in kaiju_names_full_fh:
        line = line.strip().split('\t')
        status = line[0]
        contig_id = line[1]
        node_id = contig_id.split('NODE_')[1].split('_')[0]
        if status == 'C':
            try:
                kaiju_tax_full = line[7]
            except IndexError:
                if line[2] == '1':
                    kaiju_tax_full = 'root'
            contig_kaiju_full_dict[node_id] = kaiju_tax_full
        else:
            contig_kaiju_full_dict[node_id] = 'Unassigned'

with open(gtdb_bac, 'r') as gtdb_bac_fh:
    header = gtdb_bac_fh.readline().strip().split('\t')
    for line in gtdb_bac_fh:
        line = line.strip().split('\t')
        mag = line[0]
        gtdb_tax = line[1]
        if gtdb_tax != 'N/A':
            if mag in mag_gtdb_dict.keys():
                mag_gtdb_dict[mag] = gtdb_tax
            else:
                mag_gtdb_dict[mag] = gtdb_tax

with open(gtdb_ar, 'r') as gtdb_ar_fh:
    header = gtdb_ar_fh.readline().strip().split('\t')
    for line in gtdb_ar_fh:
        line = line.strip().split('\t')
        mag = line[0]
        gtdb_tax = line[1]
        if gtdb_tax != 'N/A':
            if mag in mag_gtdb_dict.keys():
                mag_gtdb_dict[mag] = gtdb_tax
            else:
                mag_gtdb_dict[mag] = gtdb_tax

# traverse os.walk(bins_directory) and find the mag in which contig_id is present
for root, dirs, files in os.walk(bins_directory):
    for file in files:
        if file.endswith('.fa'):
            mag = file.replace('.fa', '')
            mag = mag.replace('.', '_')
            with open(os.path.join(root, file), 'r') as mag_fh:
                for line in mag_fh:
                    if line.startswith('>'):
                        contig_id = line.strip().replace('>', '')
                        node_id = contig_id.split('NODE_')[1].split('_')[0]
                        contig_bins_dict[node_id] = mag

for node in contig_bins_dict.keys():
    if contig_bins_dict[node] in mag_gtdb_dict.keys():
        contig_gtdb_dict[node] = mag_gtdb_dict[contig_bins_dict[node]]
    else:
        contig_gtdb_dict[node] = 'gtdb_unclassified'

length_all_contigs = 0
length_binned_contigs = 0
length_unbinned_contigs = 0
num_binned_contigs = 0
num_unbinned_contigs = 0

for contig_id in contig_ids:
    node_id = contig_id.split('NODE_')[1].split('_')[0]
    length = contig_length_dict[node_id]
    length_all_contigs += int(length)
    if 'unbinned' in contig_bins_dict[node_id]:
        length_unbinned_contigs += int(length)
        num_unbinned_contigs += 1
    else:
        length_binned_contigs += int(length)
        num_binned_contigs += 1

for contig_id in contig_ids:
    node_id = contig_id.split('NODE_')[1].split('_')[0]
    contig_ids_dict[node_id] = contig_id

#  write summary file in the same path as outfile with different suffix
summary_file = outfile.replace('.tsv', '_summary.tsv')
with open(summary_file, 'w') as summary_file_fh:
    summary_file_fh.write('sample\tnum_contigs\tnum_whokaryote\tnum_kaiju\tnum_binned\n')
    summary_file_fh.write(f'{sample}\t{len(contig_ids)}\t{len(contig_whokaryote_dict)}\t{len(contig_kaiju_dict)}\t{num_binned_contigs}\n')
    # print('sample\tnum_contigs\tnum_whokaryote\tnum_kaiju\tnum_binned\n')
    # print(f'{sample}\t{len(contig_ids)}\t{len(contig_whokaryote_dict)}\t{len(contig_kaiju_dict)}\t{num_binned_contigs}\n')

with open(outfile, "w+") as outfile_fh:
    outfile_fh.write('contig_id\tcontig_length\tcontig_whokaryote\tcontig_kaiju\tmag_gtdb\tcontig_bins\n')
    for contig_id in contig_ids:
        node_id = contig_id.split('NODE_')[1].split('_')[0]
        contig_length = contig_length_dict[node_id]
        if node_id in contig_whokaryote_dict.keys():
            contig_whokaryote = contig_whokaryote_dict[node_id]
        else:
            contig_whokaryote = 'whokaryote_Unclassified'
        # if node_id in contig_kaiju_dict.keys():
        #     contig_kaiju = contig_kaiju_dict[node_id]
        # else:
        #     contig_kaiju = 'Unclassified'
        if node_id in contig_kaiju_full_dict.keys():
            contig_kaiju = contig_kaiju_full_dict[node_id]
        else:
            contig_kaiju = 'kaiju_Unclassified'
        if node_id in contig_gtdb_dict.keys():
            contig_gtdb = contig_gtdb_dict[node_id]
        else:
            contig_gtdb = 'gtdb_Unclassified'
        if node_id in contig_bins_dict.keys():
            contig_bins = contig_bins_dict[node_id]
        else:
            contig_bins = 'metabat_Unclassified'
        outfile_fh.write(f'{contig_id}\t{contig_length}\t{contig_whokaryote}\t{contig_kaiju}\t{contig_gtdb}\t{contig_bins}\n')

df = pd.read_csv(outfile, sep='\t', index_col=False)

filtered_df = df
filtered_df = filtered_df[filtered_df['contig_length'].astype(int) > 10000]
# filtered_df = filtered_df[filtered_df['contig_whokaryote'] != 'eukaryote']
filtered_df['mag_gtdb'] =filtered_df['mag_gtdb'].str.extract(r'^(?:.*;)?(?:s__|g__|f__)?([^;]+)$')
filtered_df['contig_kaiju'] = filtered_df['contig_kaiju'].str.rsplit(';', n=2, expand=True)[1].str.strip()


links = []
nodes = []
node_ids = {}
node_x = []
# node_y = []

for _, row in filtered_df.iterrows():
    source_node = row['contig_whokaryote']
    target_node = row['contig_bins']
    kaiju_node = row['contig_kaiju']

    if source_node not in nodes:
        nodes.append(source_node)
        node_ids[source_node] = len(nodes) - 1
        node_x.append(0)
        # node_y.append(len(nodes) - 1)

    if target_node not in nodes:
        nodes.append(target_node)
        node_ids[target_node] = len(nodes) - 1
        node_x.append(1)
        # node_y.append(len(nodes) - 1)

    if kaiju_node not in nodes:
        nodes.append(kaiju_node)
        node_ids[kaiju_node] = len(nodes) - 1
        node_x.append(2)
        # node_y.append(len(nodes) - 1)

    links.append(
        dict(
            source=node_ids[source_node],
            target=node_ids[target_node],
            value=row['contig_length'],
            hovertemplate='Source: %{source.label}<br>' +
                          'Target: %{target.label}<br>' +
                          'Value: %{value}<br>' +
                          'Contig Length: %{value}',
        )
    )
    links.append(
        dict(
            source=node_ids[target_node],
            target=node_ids[kaiju_node],
            value=row['contig_length'],
            hovertemplate='Source: %{source.label}<br>' +
                          'Target: %{target.label}<br>' +
                          'Value: %{value}<br>' +
                          'Contig Length: %{value}',
        )
    )

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=10,
        line=dict(color='black', width=0.5),
        label=nodes,
        x=node_x
        # y=node_y
    ),
    link=dict(
        source=[link['source'] for link in links],
        target=[link['target'] for link in links],
        value=[link['value'] for link in links],
        hovertemplate='Source: %{source.label}<br>' +
                      'Target: %{target.label}<br>' +
                      'Value: %{value}<br>' +
                      'Contig Length: %{value}',
    ),
)])

fig.update_layout(
    title='Contig Whokaryote, Contig Bins, and Contig Kaiju Sankey Diagram (Contig Length > 5000)',
    height=2000,
    font=dict(size=8)
)

fig.show()
fig.write_html(outfile.replace('.tsv', '_contig_classification_sankey.html'))

filtered_df = df
# filtered_df = df[df['contig_length'].astype(int) > 5000]  # Filter contigs longer than 10,000
filtered_df = filtered_df[filtered_df['contig_whokaryote'] == 'prokaryote']  # Filter contigs longer than 10,000
filtered_df = filtered_df[~filtered_df['contig_bins'].str.contains('unbinned')]
filtered_df['mag_gtdb'] =filtered_df['mag_gtdb'].str.extract(r'^(?:.*;)?(?:s__|g__|f__)?([^;]+)$')
filtered_df['contig_kaiju'] = filtered_df['contig_kaiju'].str.rsplit(';', n=2, expand=True)[1].str.strip()


links = []
nodes = []
node_ids = {}
node_x = []
# node_y = []

for _, row in filtered_df.iterrows():
    source_node = row['contig_whokaryote']
    target_node = row['contig_bins']
    kaiju_node = row['contig_kaiju']

    if source_node not in nodes:
        nodes.append(source_node)
        node_ids[source_node] = len(nodes) - 1
        node_x.append(0)
        # node_y.append(len(nodes) - 1)

    if target_node not in nodes:
        nodes.append(target_node)
        node_ids[target_node] = len(nodes) - 1
        node_x.append(1)
        # node_y.append(len(nodes) - 1)

    if kaiju_node not in nodes:
        nodes.append(kaiju_node)
        node_ids[kaiju_node] = len(nodes) - 1
        node_x.append(2)
        # node_y.append(len(nodes) - 1)

    links.append(
        dict(
            source=node_ids[source_node],
            target=node_ids[target_node],
            value=row['contig_length'],
            hovertemplate='Source: %{source.label}<br>' +
                          'Target: %{target.label}<br>' +
                          'Value: %{value}<br>' +
                          'Contig Length: %{value}',
        )
    )
    links.append(
        dict(
            source=node_ids[target_node],
            target=node_ids[kaiju_node],
            value=row['contig_length'],
            hovertemplate='Source: %{source.label}<br>' +
                          'Target: %{target.label}<br>' +
                          'Value: %{value}<br>' +
                          'Contig Length: %{value}',
        )
    )

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=10,
        line=dict(color='black', width=0.5),
        label=nodes,
        x=node_x
        # y=node_y
    ),
    link=dict(
        source=[link['source'] for link in links],
        target=[link['target'] for link in links],
        value=[link['value'] for link in links],
        hovertemplate='Source: %{source.label}<br>' +
                      'Target: %{target.label}<br>' +
                      'Value: %{value}<br>' +
                      'Contig Length: %{value}',
    ),
)])

fig.update_layout(
    title='Contig Whokaryote, Contig Bins, and Contig Kaiju Sankey Diagram (Contig Length > 5000)',
    height=800,
    font=dict(size=8)
)

fig.show()
fig.write_html(outfile.replace('.tsv', '_prokaryote_contig_classification_sankey.html'))

# Create the sunburst plot
df = pd.read_csv(outfile, sep='\t', index_col=False)
filtered_df = df
filtered_df = df[df['contig_length'].astype(int) > 5000]  # Filter contigs longer than 10,000
filtered_df['mag_gtdb'] =filtered_df['mag_gtdb'].str.extract(r'^(?:.*;)?(?:s__|g__|f__)?([^;]+)$')
filtered_df['contig_kaiju'] = filtered_df['contig_kaiju'].str.rsplit(';', n=2, expand=True)[1].str.strip()
filtered_df['mag_gtdb'].fillna('gtdb_Unclassified', inplace=True)
filtered_df['contig_kaiju'].fillna('kaiju_Unclassified', inplace=True)

fig_sun = px.sunburst(
    filtered_df,
    path=['contig_whokaryote', 'contig_bins', 'mag_gtdb', 'contig_kaiju'],
    values='contig_length',
    color='contig_bins',
    hover_data={
        'mag_gtdb': True,
        'contig_kaiju': True,
        'contig_whokaryote': True,
        'contig_id': True,
        'contig_length': True
    }
)

# Customize the appearance of the sunburst plot
fig_sun.update_layout(
    title='Contig Classification Sunburst Plot (Contig Length > 5000)',
    height=800  # Adjust the height as needed
)

# Display the sunburst plot
fig_sun.show()
fig_sun.write_html(outfile.replace('.tsv', '_contig_classification_sunburst.html'))