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
from scipy.stats import rankdata


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

ab_mat = pd.read_csv('results/figures/magOTU_abundance_matrix_copies.csv', index_col=0)
rel_mat = pd.read_csv('results/figures/magOTU_abundance_matrix_rel.csv', index_col=0)
pa_mat = pd.read_csv('results/figures/magOTU_abundance_matrix_presence_absence.csv', index_col=0)

samples = pa_mat.index.tolist()
print(f'a total of {len(samples)} are included in the analysis')
print(f'from mellifera we have {len([sample for sample in samples if host_spec(sample) == "mellifera"])} samples')
print(f'from cerana we have {len([sample for sample in samples if host_spec(sample) == "cerana"])} samples')
print(f'from dorsata we have {len([sample for sample in samples if host_spec(sample) == "dorsata"])} samples')
print(f'from florea we have {len([sample for sample in samples if host_spec(sample) == "florea"])} samples')
print(f'from andreniformis we have {len([sample for sample in samples if host_spec(sample) == "andreniformis"])} samples')

# using data from India and Malaysia only

# # only keep samples from Malaysia for now
# samples = [sample for sample in samples if location(sample) == 'Malaysia']
# print(f'after filtering for Malaysia, we have {len(samples)} samples')
# print(f'from mellifera we have {len([sample for sample in samples if host_spec(sample) == "mellifera"])} samples')
# print(f'from cerana we have {len([sample for sample in samples if host_spec(sample) == "cerana"])} samples')
# print(f'from dorsata we have {len([sample for sample in samples if host_spec(sample) == "dorsata"])} samples')
# print(f'from florea we have {len([sample for sample in samples if host_spec(sample) == "florea"])} samples')
# print(f'from andreniformis we have {len([sample for sample in samples if host_spec(sample) == "andreniformis"])} samples')

'''
for each magotu (column) we want to make a plot of the proportion of
sample pairs of the same host species that share the magotu on the y axis
and the proportion of sample pairs of different host species that share the magotu on the x axis
the host species of each sample (row) is determined by the first two letters of the sample name
'''

# make a dictionary of host species for each sample
host_dict = {}
for sample in samples:
    host_dict[sample] = host_spec(sample)


host_combos = []
for host1, host2 in combinations(set(host_dict.values()), 2):
    host_combos.append('_'.join(sorted([host1, host2])))
for host in set(host_dict.values()):
    host_combos.append('_'.join([host, host]))

final_dict = {}
for magotu in pa_mat.columns:
    shared_dict = {}
    totals_dict = {}
    for comb in host_combos:
        shared_dict[comb] = 0
        totals_dict[comb] = 0
    
    combinations_seen = set()
    # populate it with sample sample combinations for all samples in the samples
    for sample in samples:
        combinations_seen.add(frozenset([sample, sample]))
    
    final_dict[magotu] = {}

    for sample1, sample2 in combinations(samples, 2):
        if frozenset([sample1, sample2]) in combinations_seen:
            continue
        combinations_seen.add(frozenset([sample1, sample2]))
        comb_i = '_'.join(sorted([host_dict[sample1], host_dict[sample2]]))
        totals_dict[comb_i] += 1
        if pa_mat.loc[sample1, magotu] == pa_mat.loc[sample2, magotu] == 1:
            shared_dict[comb_i] += 1
    
    print(f'{magotu}')
    for comb in host_combos:
        print(f'{comb}: {shared_dict[comb]} / {totals_dict[comb]}')
        final_dict[magotu][comb] = shared_dict[comb] / totals_dict[comb]

shared_df_plot = pd.DataFrame(final_dict).T
shared_df_plot.to_csv('results/figures/magotu_shared_same_diff.csv')

final_dict_ab = {}
for magotu in pa_mat.columns:
    shared_dict = {}
    totals_dict = {}
    for comb in host_combos:
        shared_dict[comb] = 0
        totals_dict[comb] = 0
    
    combinations_seen = set()
    # populate it with sample sample combinations for all samples in the samples
    for sample in samples:
        combinations_seen.add(frozenset([sample, sample]))
    
    final_dict_ab[magotu] = {}

    for sample1, sample2 in combinations(samples, 2):
        if frozenset([sample1, sample2]) in combinations_seen:
            continue
        combinations_seen.add(frozenset([sample1, sample2]))
        comb_i = '_'.join(sorted([host_dict[sample1], host_dict[sample2]]))
        totals_dict[comb_i] += 1
        # if the log of the ratio of higher one and lower one is
        # smaller than 1, we count them as shared
        if np.log2(ab_mat.loc[sample1, magotu] / ab_mat.loc[sample2, magotu]) < 1:
            shared_dict[comb_i] += 1
        elif np.log2(ab_mat.loc[sample2, magotu] / ab_mat.loc[sample1, magotu]) < 1:
            shared_dict[comb_i] += 1
    
    print(f'{magotu}')
    for comb in host_combos:
        print(f'{comb}: {shared_dict[comb]} / {totals_dict[comb]}')
        final_dict_ab[magotu][comb] = shared_dict[comb] / totals_dict[comb]

shared_df_plot_ab = pd.DataFrame(final_dict_ab).T
shared_df_plot_ab.to_csv('results/figures/magotu_shared_same_diff_ab.csv')


# estimating rhodes std. index using prevalence
'''
Based on https://hal.sorbonne-universite.fr/hal-03246067/document
rohde's standard index of host specificity
Specificity of species i (there are j = 5 host species)
min_Si = 1/5*(1/1 + 1/2 + 1/3 + 1/4 + 1/5)
S_i_ratio = sum_j(x/nh) / sum_j(x/n)
# if p is prevalence,
S_i_ratio = sum_j(p/h) / sum_j(p)
S_i = (S_i_ratio - min_Si) / (1 - min_Si)
x is the number of individuals of j that have i, n is the number of individuals of j measured
and h is the rank in terms of prevalence of i in j of host j
'''

#  make a dictonary of magotu and rank for each host species, if two hosts have the same rank,

rank_dict = {}
prevalence_in_host = {}
for magotu in pa_mat.columns:
    rank_dict[magotu] = {}
    prevalence_in_host[magotu] = {}
    for host in set(host_dict.values()):
        # get the prevalence for each host using host_spec function
        prevalences = pa_mat.loc[[sample for sample in samples if host_spec(sample) == host], magotu]
        # get the rank of the prevalence
        prevalence_in_host[magotu][host] = sum(prevalences) / len(prevalences)
    # ties are broken by assigning the same rank to each of the tied values, with the next value(s) receiving the immediately following rank
    # this is to ensure that they contribute equally to the index
    rank_dict[magotu] = dict(zip(prevalence_in_host[magotu].keys(), rankdata([-i for i in prevalence_in_host[magotu].values()], method='min')))

min_Si = 1/5*(1/1 + 1/2 + 1/3 + 1/4 + 1/5)
rohdes_dict = {}
for magotu in pa_mat.columns:
    rohdes_dict[magotu] = None
    S_i_ratio_numerator = 0
    S_i_ratio_denominator = 0
    for host in set(host_dict.values()):
        p = prevalence_in_host[magotu][host]
        h = rank_dict[magotu][host]
        S_i_ratio_numerator += p/h
        S_i_ratio_denominator += p
    S_i_ratio = S_i_ratio_numerator / S_i_ratio_denominator
    S_i = (S_i_ratio - min_Si) / (1 - min_Si)
    rohdes_dict[magotu] = S_i

ranks_df = pd.DataFrame.from_dict(rank_dict, orient='index')
ranks_df.to_csv('results/figures/magotu_prev_ranks.csv')

# estimating rhodes std. index using intensity (abundance)
'''
Specificity of species i (there are j = 5 host species)
min_Si = 1/5*(1/1 + 1/2 + 1/3 + 1/4 + 1/5)
S_i_ratio = sum_j(x/nh) / sum_j(x/n)
S_i = (S_i_ratio - min_Si) / (1 - min_Si)
x is the number of individuals i in host j (take the mean?) that have i,
n is the number of individuals of host j measured (x/n is the <mean>intensity...? of infection)
and h is the rank in terms of prevalence of i in j of host j
## use mean not median - median makes 0s in many cases where only a few individuals are infected and makes those look like a value of 2 even though they are individuals of different host species##
'''

rank_dict = {}
intensity_in_host = {}
for magotu in ab_mat.columns:
    rank_dict[magotu] = {}
    intensity_in_host[magotu] = {}
    for host in set(host_dict.values()):
        # get the intensity for each host using host_spec function
        intensitys = ab_mat.loc[[sample for sample in samples if host_spec(sample) == host], magotu]
        # get the rank of the intensity
        intensity_in_host[magotu][host] = np.mean(intensitys)
    # ties are broken by assigning the same rank to each of the tied values, with the next value(s) receiving the immediately following rank
    # this is to ensure that they contribute equally to the index
    rank_dict[magotu] = dict(zip(intensity_in_host[magotu].keys(), rankdata([-i for i in intensity_in_host[magotu].values()], method='min')))

min_Si = 1/5*(1/1 + 1/2 + 1/3 + 1/4 + 1/5)
rohdes_dict_ab = {}
for magotu in ab_mat.columns:
    rohdes_dict_ab[magotu] = None
    S_i_ratio_numerator = 0
    S_i_ratio_denominator = 0
    for host in set(host_dict.values()):
        x = intensity_in_host[magotu][host]
        n = len(ab_mat.loc[[sample for sample in samples if host_spec(sample) == host], magotu])
        h = rank_dict[magotu][host]
        S_i_ratio_numerator += x/n/h
        S_i_ratio_denominator += x/n
    S_i_ratio = S_i_ratio_numerator / S_i_ratio_denominator
    S_i = (S_i_ratio - min_Si) / (1 - min_Si)
    rohdes_dict_ab[magotu] = S_i

ranks_df = pd.DataFrame.from_dict(rank_dict, orient='index')
ranks_df.to_csv('results/figures/magotu_ab_ranks.csv')

# estimating rhodes std. index using relative intensity (abundance)
'''
Specificity of species i (there are j = 5 host species)
min_Si = 1/5*(1/1 + 1/2 + 1/3 + 1/4 + 1/5)
S_i_ratio = sum_j(x/nh) / sum_j(x/n)
S_i = (S_i_ratio - min_Si) / (1 - min_Si)
x is the number of individuals i in host j (take the mean?) that have i,
n is the number of individuals of host j measured (x/n is the <mean>intensity...? of infection)
and h is the rank in terms of prevalence of i in j of host j
'''

rank_dict = {}
rel_intensity_in_host = {}
for magotu in rel_mat.columns:
    rank_dict[magotu] = {}
    rel_intensity_in_host[magotu] = {}
    for host in set(host_dict.values()):
        # get the intensity for each host using host_spec function
        intensitys = rel_mat.loc[[sample for sample in samples if host_spec(sample) == host], magotu]
        # get the rank of the intensity
        rel_intensity_in_host[magotu][host] = np.mean(intensitys)
    # ties are broken by assigning the same rank to each of the tied values, with the next value(s) receiving the immediately following rank
    # this is to ensure that they contribute equally to the index
    rank_dict[magotu] = dict(zip(rel_intensity_in_host[magotu].keys(), rankdata([-i for i in rel_intensity_in_host[magotu].values()], method='min')))

min_Si = 1/5*(1/1 + 1/2 + 1/3 + 1/4 + 1/5)
rohdes_dict_rel = {}
for magotu in rel_mat.columns:
    rohdes_dict_rel[magotu] = None
    S_i_ratio_numerator = 0
    S_i_ratio_denominator = 0
    for host in set(host_dict.values()):
        x = rel_intensity_in_host[magotu][host]
        n = len(rel_mat.loc[[sample for sample in samples if host_spec(sample) == host], magotu])
        h = rank_dict[magotu][host]
        S_i_ratio_numerator += x/n/h
        S_i_ratio_denominator += x/n
    S_i_ratio = S_i_ratio_numerator / S_i_ratio_denominator
    S_i = (S_i_ratio - min_Si) / (1 - min_Si)
    rohdes_dict_rel[magotu] = S_i

ranks_df = pd.DataFrame.from_dict(rank_dict, orient='index')
ranks_df.to_csv('results/figures/magotu_rel_ranks.csv')

# plot the rhodes std. index for each type with each magotu
# make a df with rel, ab and prev indices as type
rohdes_df = pd.DataFrame.from_dict(rohdes_dict, orient='index', columns=['prev'])
rohdes_df['ab'] = pd.Series(rohdes_dict_ab)
rohdes_df['rel'] = pd.Series(rohdes_dict_rel)
rohdes_df.to_csv('results/figures/magotu_rohdes_index.csv')