import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from Bio import SeqIO
from pprint import pprint

'''

confirm that the max covered host genome corresponds to the host species that the sample was collected from

coverage out put is at:
./results/03_host_mapping/{sample}_coverage.tsv

looks like this:

#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
GL575021.1	1	4900	4	41	0.836735	0.0167347	35.5	27.5
GL575022.1	1	2968	15	189	6.36792	0.204515	36.2	12.8
GL575023.1	1	3590	15	154	4.28969	0.257382	34.6	2.6
GL575024.1	1	3301	7	228	6.907	0.3196	35.6	0

chromosome names are like this:
>AUPE... dorsata (contig)
>GL57... florea
>JANR... cerana
>KI2... dorsata (scaffold)
>NC... mellifera
>NW... mellifera
>NC_001566.1 Apis mellifera ligustica mitochondrion, complete genome

'''

samples = []

for sample in [os.path.basename(x).split('_cov')[0] for x in glob.glob('results/03_host_mapping/*_coverage.tsv')]:
    samples.append(sample)

tag_host_dict = {'AU': 'Apis dorsata',
                 'GL': 'Apis florea',
                 'JA': 'Apis cerana',
                 'KI': 'Apis dorsata',
                 'NC': 'Apis mellifera',
                 'NW': 'Apis mellifera_mito'}

host_tag_dict = {'Apis dorsata': 'AU',
                 'Apis florea': 'GL',
                 'Apis cerana': 'JA',
                 'Apis dorsata': 'KI',
                 'Apis mellifera': 'NC',
                 'Apis mellifera_mito': 'NW'}

lengths = {}
for value in tag_host_dict.values():
    lengths[value] = np.array([])

# summarise number of and length of contigs

for record in SeqIO.parse('data/host_database/apis_bees_db.fasta', 'fasta'):
    tag = record.id.split()[0]
    host = tag_host_dict[tag[:2]]
    lengths[host] = np.append(lengths[host], len(record.seq))

# number of contigs for each host species

for key, value in lengths.items():
    print(f'{key},\n number: {len(value)},\n median length: {np.median(value)}, \n total length: {np.sum(value):,}')

depths = {}
covs = {}

for sample in samples:
    covs[sample] = {}
    depths[sample] = {}

for sample in samples:
    df = pd.read_csv(f'results/03_host_mapping/{sample}_coverage.tsv', sep='\t')
    # get cov and depth for each host species
    # mutate to add a host column based on the chromosome name
    df['host'] = df['#rname'].str[:2].apply(lambda x: tag_host_dict[x])
    # get depth for each host species
    for host in tag_host_dict.values():
        depths[sample][host] = df[df['host'] == host]['meandepth'].median()
        covs[sample][host] = df[df['host'] == host]['coverage'].median()
    
os.makedirs('results/figures/host_validation')

for sample in samples:
    fig, ax = plt.subplots()
    ax.bar([x-0.2 for x in range(len(depths[sample]))], depths[sample].values(), color='#e41a1c', width=0.45)
    ax.set_ylabel('depth', color='#e41a1c')
    ax.set_xticks(range(len(depths[sample])))
    ax.set_xticklabels(depths[sample].keys(), rotation=45, ha='right')
    ax.tick_params(axis='y')
    # plot coverage
    ax2 = ax.twinx()
    ax2.bar([x+0.2 for x in range(len(covs[sample]))], covs[sample].values(), color='#377eb8', width=0.45)
    # angle x tick labels
    ax2.set_ylabel('coverage', color='#377eb8')
    ax2.tick_params(axis='y')
    ax2.set_ylim(0, 100)
    ax.set_title(f'{sample}')
    plt.tight_layout()
    plt.savefig(f'results/figures/host_validation/{sample}_depth_cov.png')
    plt.close()

# merge depths and covs to one df with one row for each sample
depths_df = pd.DataFrame(depths)
covs_df = pd.DataFrame(covs)
depths_df = depths_df.transpose()
covs_df = covs_df.transpose()
# add a columns to each called max_host one of the species exclude Apis mellifera_mito
depths_df['max_host'] = depths_df.drop('Apis mellifera_mito', axis=1).idxmax(axis=1)
covs_df['max_host'] = covs_df.drop('Apis mellifera_mito', axis=1).idxmax(axis=1)
# append _cov to all the column names
depths_df.columns = [x + '_depth' for x in depths_df.columns]
covs_df.columns = [x + '_cov' for x in covs_df.columns]
df_to_write_to_csv = pd.concat([depths_df, covs_df], axis=1)
df_to_write_to_csv.to_csv('results/figures/host_validation.csv')



