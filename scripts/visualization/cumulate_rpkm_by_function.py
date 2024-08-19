import os
import shutil
import sys
import io
import json
import gzip
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pickle
from collections import Counter
from pprint import pprint
from Bio import SeqIO
import random
from itertools import product, combinations, chain
from scipy import stats
import numpy as np
import random as r

samples_to_exclude = ["F4-5", "F5-1", "M6-2", "D9-5", "F7-5"]

samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]
samples = [x for x in samples if x not in samples_to_exclude]

'''
This script writes matrices than can be used for PCoA of KOs/GHs etc of a certain family of functions (eg. carbon metabolism)
per mag or other unit when needed. It writes it into a file that can be used for PCoA in R located in the directory
results/figures/08-summarize_functions/rpkm_matrix_by_function/
'''



'''
make a matrix to be used for PCoA of KOs of a certain family of functions (eg. carbon metabolism)
with bacterial species (and host species) information
each column is a KO, and each row is a MAG
values within are rpkm
'''

# consider KOs included in "carbon metabolism" pathways in the DRAM distill output

# make a matrix with MAGs as rows and KOs as columns with rpkm within
# for each sample, get the rpkm of the KOs in the KOs_of_interest list
# and make a row with the MAG name as the index
# and the rpkm as the values
# then concatenate all of these into a matrix
# then add the species and genus information as columns
# then add the host species information as a column
# then save this matrix

host_species = ['Apis mellifera', 'Apis cerana',
                'Apis dorsata', 'Apis florea', 'Apis andreniformis']

def host_of_sample(sample_name):
    if sample_name:
        if sample_name.startswith('M') or sample_name.startswith('Dr') or sample_name.startswith('Gr') or sample_name.startswith('Am'):
            return 'Apis mellifera'
        if sample_name.startswith('C') or sample_name.startswith('Ac'):
            return 'Apis cerana'
        if sample_name.startswith('D'):
            return 'Apis dorsata'
        if sample_name.startswith('F'):
            return 'Apis florea'
        if sample_name.startswith('A'):
            return 'Apis andreniformis'

os.makedirs('results/figures/08-summarize_functions/rpkm_matrix_by_function/', exist_ok=True)

# make the matrix for all KOs
df_collected_rpkms = pd.DataFrame()
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus', 'cazyme']]
    df_joined = df_info.merge(df_covs, on='gene')
    df_joined = df_joined[df_joined['ko'].notnull()]
    # get summed count per ko
    df_joined = df_joined.groupby(['ko', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined = df_joined.pivot(index='mag', columns='ko', values='rpkm').reset_index()
    df_collected_rpkms = pd.concat([df_collected_rpkms, df_joined])
df_collected_rpkms = df_collected_rpkms.fillna(0)
# remove columns with all 0s
df_collected_rpkms = df_collected_rpkms.loc[:, (df_collected_rpkms != 0).any(axis=0)]
df_collected_rpkms.to_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_all_kos.csv', index=False)

# normalize rpkm for absolute abundance by multiplying by the norm factor
qpcr_df = pd.read_csv('results/figures/qpcr_plot_df.csv')
min_copies = qpcr_df['copy_num'].min()
# # if there is no value for copy number for a sample, use the median for its host
# [x/y for x,y in zip(qpcr_df['host_median'], qpcr_df['norm_factor'])]
qpcr_df['norm_factor'] = qpcr_df['copy_num'] / min_copies
qpcr_df['host_median'] = qpcr_df.groupby('Host')['norm_factor'].transform('median')
qpcr_df['norm_factor'] = qpcr_df['norm_factor'].fillna(qpcr_df['host_median'])
df_collected_rpkms = pd.DataFrame()
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus', 'cazyme']]
    df_joined = df_info.merge(df_covs, on='gene')
    df_joined = df_joined[df_joined['ko'].notnull()]
    df_joined['sample'] = sample
    if sample in qpcr_df['Sample.Name'].values:
        norm_factor = qpcr_df[qpcr_df['Sample.Name'] == sample]['norm_factor'].values[0]
    else:
        norm_factor = qpcr_df[qpcr_df['Host'] == host_of_sample(sample)]['host_median'].values[0]
    df_joined['norm_factor'] = norm_factor
    df_joined['rpkm'] = df_joined['rpkm'] * df_joined['norm_factor']
    # get summed count per ko
    df_joined = df_joined.groupby(['ko', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined = df_joined.pivot(index='mag', columns='ko', values='rpkm').reset_index()
    df_collected_rpkms = pd.concat([df_collected_rpkms, df_joined])
df_collected_rpkms = df_collected_rpkms.fillna(0)
# remove columns with all 0s
df_collected_rpkms = df_collected_rpkms.loc[:, (df_collected_rpkms != 0).any(axis=0)]
df_collected_rpkms.to_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_abs_all_kos.csv', index=False)



# make a matrix of only the KOs of interest (carbon metabolism - as defined in DRAM)
DRAM_modules = pd.read_csv('scripts/visualization/DRAM_KOs_list_from_excel.txt', sep = '\t', header = None) # copied and paseted from DRAM metabolism output
DRAM_modules.columns = ['ko', 'enzyme', 'module', 'module_type', 'category']
# only keep those with module type in the list of interesting ones
DRAM_modules = DRAM_modules[DRAM_modules['module_type'].isin(['central carbon'])]
KOs_of_interest = [x for x in DRAM_modules['ko'].unique() if x.startswith('K') and len(x) == 6]

df_collected_rpkms = pd.DataFrame()
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus', 'cazyme']]
    df_joined = df_info.merge(df_covs, on='gene')
    df_joined = df_joined[df_joined['ko'].notnull()]
    df_joined = df_joined[df_joined['ko'].isin(KOs_of_interest)]
    # get summed count per ko
    df_joined = df_joined.groupby(['ko', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined = df_joined.pivot(index='mag', columns='ko', values='rpkm').reset_index()
    # df_joined['sample'] = sample
    df_collected_rpkms = pd.concat([df_collected_rpkms, df_joined])
df_collected_rpkms = df_collected_rpkms.fillna(0)
# remove columns with all 0s
df_collected_rpkms = df_collected_rpkms.loc[:, (df_collected_rpkms != 0).any(axis=0)]
df_collected_rpkms.to_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_central_carbon_metabolism.csv', index=False)




df_collected_rpkms = pd.DataFrame()
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus', 'cazyme']]
    df_joined = df_info.merge(df_covs, on='gene')
    df_joined = df_joined[df_joined['cazyme'].notnull()]
    # get summed count per cazyme
    df_joined = df_joined.groupby(['cazyme', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined = df_joined.pivot(index='mag', columns='cazyme', values='rpkm').reset_index()
    df_collected_rpkms = pd.concat([df_collected_rpkms, df_joined])
df_collected_rpkms = df_collected_rpkms.fillna(0)
# remove columns with all 0s
df_collected_rpkms = df_collected_rpkms.loc[:, (df_collected_rpkms != 0).any(axis=0)]
df_collected_rpkms.to_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_all_cazy.csv', index=False)



DRAM_modules = pd.read_csv('scripts/visualization/DRAM_KOs_list_from_excel.txt', sep = '\t', header = None) # copied and paseted from DRAM metabolism output
DRAM_modules.columns = ['ko', 'enzyme', 'module', 'module_type', 'category']
# only keep those with module type in the list of interesting ones
DRAM_modules = DRAM_modules[DRAM_modules['category'].isin(['Transporters'])]
KOs_of_interest = [x for x in DRAM_modules['ko'].unique() if x.startswith('K') and len(x) == 6]

df_collected_rpkms = pd.DataFrame()
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus', 'cazyme']]
    df_joined = df_info.merge(df_covs, on='gene')
    df_joined = df_joined[df_joined['ko'].notnull()]
    df_joined = df_joined[df_joined['ko'].isin(KOs_of_interest)]
    # get summed count per ko
    df_joined = df_joined.groupby(['ko', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined = df_joined.pivot(index='mag', columns='ko', values='rpkm').reset_index()
    # df_joined['sample'] = sample
    df_collected_rpkms = pd.concat([df_collected_rpkms, df_joined])
df_collected_rpkms = df_collected_rpkms.fillna(0)
# remove columns with all 0s
df_collected_rpkms = df_collected_rpkms.loc[:, (df_collected_rpkms != 0).any(axis=0)]
df_collected_rpkms.to_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_central_transporters.csv', index=False)



# get all the KOs associated with secretion systems
# 2 possible ways, one is "grepping" for secretion system or the like in the DRAM annotation output
# the other is getting it from minpath data.. try both and check the differeces

collect_sec_kos = set()
for sample in samples:
    annatations_file = f'results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv'
    with open(annatations_file, 'r') as f:
        for line in f:
            if 'secretion system' in line:
                ko = line.split('\t')[3]
                if ko.startswith('K'):
                    collect_sec_kos.add(ko)

# RESUME HERE
pathway_info = pd.read_csv('scripts/visualization/kegg_pathway_categories.tsv', sep='\t')
# read and make into a dictionary
pathway_numbers_dict = pd.read_csv('scripts/MinPath/data/KEGG-pathway.txt', sep = '\t', header = None).set_index(0).to_dict()[1]
# for each pathway num make a list of KOs in it
pathway_ko_map_df = pd.read_csv('scripts/MinPath/data/KEGG-mapping.txt', sep = '\t', header = None)
pathway_ko_map = {}
for i, row in pathway_ko_map_df.iterrows():
    name_pathway = pathway_numbers_dict[row[0]]
    if name_pathway not in pathway_ko_map.keys():
        pathway_ko_map[name_pathway] = []
    pathway_ko_map[name_pathway].append(row[1])
# All "carbon metabolism" genes
# get the KOs for carbon metabolism from minpath
functions_of_interest = ['1.1 Carbohydrate metabolism']
# functions_of_interest = ['1.2 Energy metabolism']
pathways_of_interest = pathway_info[pathway_info['subsection'].isin(functions_of_interest)]['path_name'].unique()
KOs_of_interest = []
for path in pathways_of_interest:
    if path not in pathway_ko_map.keys():
        print(f'{path} not in pathway_ko_map')
        continue
    KOs_of_interest = KOs_of_interest + pathway_ko_map[path]


# KOs_of_interest = 

df_collected_rpkms = pd.DataFrame()
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus', 'cazyme']]
    df_joined = df_info.merge(df_covs, on='gene')
    df_joined = df_joined[df_joined['ko'].notnull()]
    df_joined = df_joined[df_joined['ko'].isin(KOs_of_interest)]
    # get summed count per ko
    df_joined = df_joined.groupby(['ko', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined = df_joined.pivot(index='mag', columns='ko', values='rpkm').reset_index()
    # df_joined['sample'] = sample
    df_collected_rpkms = pd.concat([df_collected_rpkms, df_joined])
df_collected_rpkms = df_collected_rpkms.fillna(0)
# remove columns with all 0s
df_collected_rpkms = df_collected_rpkms.loc[:, (df_collected_rpkms != 0).any(axis=0)]
df_collected_rpkms.to_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_central_sec_systems.csv', index=False)

# set(DRAM_modules['module_type'])
# DRAM_modules = DRAM_modules[DRAM_modules['module_type'] == 'central carbon']

# pathway_info = pd.read_csv('scripts/visualization/kegg_pathway_categories.tsv', sep='\t')
# # read and make into a dictionary
# pathway_numbers_dict = pd.read_csv('scripts/MinPath/data/KEGG-pathway.txt', sep = '\t', header = None).set_index(0).to_dict()[1]
# # for each pathway num make a list of KOs in it
# pathway_ko_map_df = pd.read_csv('scripts/MinPath/data/KEGG-mapping.txt', sep = '\t', header = None)
# pathway_ko_map = {}
# for i, row in pathway_ko_map_df.iterrows():
#     name_pathway = pathway_numbers_dict[row[0]]
#     if name_pathway not in pathway_ko_map.keys():
#         pathway_ko_map[name_pathway] = []
#     pathway_ko_map[name_pathway].append(row[1])
# # All "carbon metabolism" genes
# # get the KOs for carbon metabolism from minpath
# functions_of_interest = ['1.1 Carbohydrate metabolism']
# # functions_of_interest = ['1.2 Energy metabolism']
# pathways_of_interest = pathway_info[pathway_info['subsection'].isin(functions_of_interest)]['path_name'].unique()
# KOs_of_interest = []
# for path in pathways_of_interest:
#     if path not in pathway_ko_map.keys():
#         print(f'{path} not in pathway_ko_map')
#         continue
#     KOs_of_interest = KOs_of_interest + pathway_ko_map[path]

# Organic Nitrogen
# Transporters
# carbon utilization
# Energy
# carbon utilization (Woodcroft)
# Organic Nitrogen
