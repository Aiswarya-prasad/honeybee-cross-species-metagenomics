import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm
import glob
import pickle
from collections import Counter
from itertools import combinations
from PyComplexHeatmap import *

figpath = 'results/figures/07-Orthogroup_heatmaps'

og_genes_dict = {}
'''
og_genes_dict is a dictionary with the following structure:
{og_id: {genome_id: [gene_id, gene_id, ...], genome_id: [gene_id, gene_id, ...], ...}, ...}
'''

# Load orthogroup data
mags_info = pd.read_csv('results/09_MAGs_collection/MAGs_metadata_summary.tsv', sep='\t')
mags_info['Strain'] = mags_info['magOTU']
mags_info['Type'] = 'MAG'
ref_info = pd.read_csv('config/Outgroup_isolate_genomes.tsv', sep='\t')
ref_info['Genus'] = ref_info['Group']
ref_info['Type'] = 'reference'
ref_info['ID'] = ref_info['Strain']
combined_info = pd.concat([mags_info, ref_info], axis=0)

for genus in [os.path.basename('/'.join(x.split('/')[:-1])) for x in glob.glob('results/11_phylogenies/02_orthofinder_results/*/')]:
    if not os.path.exists(f'results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Orthogroups/Orthogroups.txt'):
        continue
    og_counts = pd.read_csv(f'results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Orthogroups/Orthogroups.GeneCount.tsv', sep='\t')
    og_counts = og_counts[og_counts['Total'] >= 10]
    # at least 5 non-zero values
    og_counts = og_counts[(og_counts != 0).sum(axis=1) >= 25]
    # remove column for Total
    og_counts.drop(columns=['Total'], inplace=True)
    # add some info about mags after melting
    og_counts = og_counts.melt(id_vars=['Orthogroup'], var_name='MAG', value_name='GeneCount')
    # add some info about reference genomes
    og_counts = og_counts.merge(combined_info[['Genus', 'ID', 'magOTU', 'Completeness', 'Type']], how='left', left_on='MAG', right_on='ID')
    # make heatmap of mag vs OG
    og_counts = og_counts.pivot(index='Orthogroup', columns='MAG', values='GeneCount')
    # remove rows with all NaNs
    og_counts = og_counts.dropna(how='all')
    plt.figure(figsize=(10, 10))
    # add annotation to heatmap completeness and magotu
    sns.clustermap(og_counts, cmap='viridis', col_cluster=False, row_cluster=True, yticklabels=True, xticklabels=True, figsize=(10, 10))
    
    # shrink axis text
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.title(f'Orthogroup gene content for {genus}')
    plt.tight_layout()
    plt.savefig(f'{figpath}/{genus}_orthogroup_heatmap.png', dpi=300)



## within each genus  
# how many orthogroups in mags of each host species?
# Do they share the same orthogroups regardless?
# are there host-species specific orthogroups?
# are some of these orthogroups annotates in the cazyme database?

# anno with ghostkoala and then check more

