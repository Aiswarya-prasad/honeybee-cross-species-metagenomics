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
from itertools import combinations_with_replacement
from PyComplexHeatmap import *


'''
read the abundance matrices per magOTU and perform various analyses
results/figures/magOTU_abundance_matrix_coverage.csv
results/figures/magOTU_abundance_matrix_presence_absence.pdf

samples_info

get host from the first letters of the samples
Gr, Dr, Am, M are mellifera
Ac, C are cerana
D are dorsata
F are florea
A are andreniformis

for location
Dr, Gr Switzerland
Am, Ac from Japan

M[1] from India
M[2-7] from Malaysia

C[1-3] from India
C[4-9] from Malaysia

D[1-3] from India
D[4-9] from Malaysia

F[1-3] from India
F[4-9] from Malaysia

A[1-6] fom Malaysia
'''

def host_spec(sample_name):
    if sample_name[0:2] in ['Gr', 'Dr', 'Am', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7']:
        return 'mellifera'
    elif sample_name[0:2] in ['Ac', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']:
        return 'cerana'
    elif sample_name[0:2] in ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9']:
        return 'dorsata'
    elif sample_name[0:2] in ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9']:
        return 'florea'
    elif sample_name[0:2] in ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']:
        return 'andreniformis'
    else:
        return 'unknown'

def location(sample_name):
    if sample_name[0:2] in ['Gr', 'Dr']:
        return 'Switzerland'
    elif sample_name[0:2] in ['Am', 'Ac']:
        return 'Japan'
    elif sample_name[0:2] in ['M1']:
        return 'India'
    elif sample_name[0:2] in ['M2', 'M3', 'M4', 'M5', 'M6', 'M7']:
        return 'Malaysia'
    elif sample_name[0:2] in ['C1', 'C2', 'C3']:
        return 'India'
    elif sample_name[0:2] in ['C4', 'C5', 'C6', 'C7', 'C8', 'C9']:
        return 'Malaysia'
    elif sample_name[0:2] in ['D1', 'D2', 'D3']:
        return 'India'
    elif sample_name[0:2] in ['D4', 'D5', 'D6', 'D7', 'D8', 'D9']:
        return 'Malaysia'
    elif sample_name[0:2] in ['F1', 'F2', 'F3']:
        return 'India'
    elif sample_name[0:2] in ['F4', 'F5', 'F6', 'F7', 'F8', 'F9']:
        return 'Malaysia'
    elif sample_name[0:2] in ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']:
        return 'Malaysia'
    else:
        return 'unknown'

ab_mat = pd.read_csv('results/figures/magOTU_abundance_matrix_coverage.csv', index_col=0)
pa_mat = pd.read_csv('results/figures/magOTU_abundance_matrix_presence_absence.csv', index_col=0)

pa_mat['g__Bifidobacterium_6_11-F7-3_10']['A1-2']

samples_measured = pa_mat.index.tolist()
print(f'a total of {samples} are included in the analysis')
print(f'from mellifera we have {len([sample for sample in samples_measured if host_spec(sample) == "mellifera"])} samples')
print(f'from cerana we have {len([sample for sample in samples_measured if host_spec(sample) == "cerana"])} samples')
print(f'from dorsata we have {len([sample for sample in samples_measured if host_spec(sample) == "dorsata"])} samples')
print(f'from florea we have {len([sample for sample in samples_measured if host_spec(sample) == "florea"])} samples')
print(f'from andreniformis we have {len([sample for sample in samples_measured if host_spec(sample) == "andreniformis"])} samples')

# only keep samples from Malaysia for now
samples_measured = [sample for sample in samples_measured if location(sample) == 'Malaysia']
print(f'after filtering for Malaysia, we have {len(samples_measured)} samples')
print(f'from mellifera we have {len([sample for sample in samples_measured if host_spec(sample) == "mellifera"])} samples')
print(f'from cerana we have {len([sample for sample in samples_measured if host_spec(sample) == "cerana"])} samples')
print(f'from dorsata we have {len([sample for sample in samples_measured if host_spec(sample) == "dorsata"])} samples')
print(f'from florea we have {len([sample for sample in samples_measured if host_spec(sample) == "florea"])} samples')
print(f'from andreniformis we have {len([sample for sample in samples_measured if host_spec(sample) == "andreniformis"])} samples')
# subset indices in samples_measured
pa_mat_malaysia = pa_mat.loc[samples_measured, ]

ab_mat_malaysia = ab_mat.loc[samples_measured, ]

'''
for each magotu (column) we want to make a plot of the proportion of
sample pairs of the same host species that share the magotu on the y axis
and the proportion of sample pairs of different host species that share the magotu on the x axis
the host species of each sample (row) is determined by the first two letters of the sample name
'''

# make a dictionary of host species for each sample
host_dict = {}
for sample in samples_measured:
    host_dict[sample] = host_spec(sample)


host_combos = []
for host1, host2 in combinations(set(host_dict.values()), 2):
    host_combos.append('_'.join(sorted([host1, host2])))
for host in set(host_dict.values()):
    host_combos.append('_'.join([host, host]))

final_dict = {}
for magotu in pa_mat_malaysia.columns:
    shared_dict = {}
    totals_dict = {}
    for comb in host_combos:
        shared_dict[comb] = 0
        totals_dict[comb] = 0
    
    combinations_seen = set()
    # populate it with sample sample combinations for all samples in the samples_measured
    for sample in samples_measured:
        combinations_seen.add(frozenset([sample, sample]))
    
    final_dict[magotu] = {}

    for sample1, sample2 in combinations(samples_measured, 2):
        if frozenset([sample1, sample2]) in combinations_seen:
            continue
        combinations_seen.add(frozenset([sample1, sample2]))
        comb_i = '_'.join(sorted([host_dict[sample1], host_dict[sample2]]))
        totals_dict[comb_i] += 1
        if pa_mat_malaysia.loc[sample1, magotu] == pa_mat_malaysia.loc[sample2, magotu] == 1:
            shared_dict[comb_i] += 1
    
    print(f'{magotu}')
    for comb in host_combos:
        print(f'{comb}: {shared_dict[comb]} / {totals_dict[comb]}')
        final_dict[magotu][comb] = shared_dict[comb] / totals_dict[comb]

shared_df_plot = pd.DataFrame(final_dict).T
shared_df_plot.to_csv('results/figures/magotu_shared_same_diff.csv')

final_dict_ab = {}
for magotu in pa_mat_malaysia.columns:
    shared_dict = {}
    totals_dict = {}
    for comb in host_combos:
        shared_dict[comb] = 0
        totals_dict[comb] = 0
    
    combinations_seen = set()
    # populate it with sample sample combinations for all samples in the samples_measured
    for sample in samples_measured:
        combinations_seen.add(frozenset([sample, sample]))
    
    final_dict_ab[magotu] = {}

    for sample1, sample2 in combinations(samples_measured, 2):
        if frozenset([sample1, sample2]) in combinations_seen:
            continue
        combinations_seen.add(frozenset([sample1, sample2]))
        comb_i = '_'.join(sorted([host_dict[sample1], host_dict[sample2]]))
        totals_dict[comb_i] += 1
        # if the log of the ratio of higher one and lower one is
        # smaller than 1, we count them as shared
        if np.log2(ab_mat_malaysia.loc[sample1, magotu] / ab_mat_malaysia.loc[sample2, magotu]) < 1:
            shared_dict[comb_i] += 1
        elif np.log2(ab_mat_malaysia.loc[sample2, magotu] / ab_mat_malaysia.loc[sample1, magotu]) < 1:
            shared_dict[comb_i] += 1
    
    print(f'{magotu}')
    for comb in host_combos:
        print(f'{comb}: {shared_dict[comb]} / {totals_dict[comb]}')
        final_dict_ab[magotu][comb] = shared_dict[comb] / totals_dict[comb]

shared_df_plot_ab = pd.DataFrame(final_dict_ab).T
shared_df_plot_ab.to_csv('results/figures/magotu_shared_same_diff_ab.csv')
