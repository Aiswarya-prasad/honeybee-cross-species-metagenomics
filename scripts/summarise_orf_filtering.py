import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

'''
This script summarises the results of the ORF filtering step.
It creates a summary table and plots the results into
results/figures/summary_ORF_filtering.png
It does not need any input files and can be run as follows:
python3 scripts/summarise_orf_filtering.py
as the result paths are hard coded in here
At the moment it was run by hand, but might be integrated into the pipeline later
'''

df_info = pd.DataFrame()
for file in glob.glob('results/06_metagenomicORFs/*/*_prodigal_filt_orfs.log.cluster.o'):
    with open(file, 'r') as fh:
        for line in fh:
            if line.startswith('sample'):
                line = fh.readline()
                line = fh.readline()
                info = line.strip().split()
                sample = info[0]
                total_records = int(info[1])
                records_filt_length = int(info[2])
                records_filt_eukaryote = int(info[3])
                records_filt = int(info[4])
                fraction = float(info[5])
                df_info = df_info._append({'sample': sample, 'total_records': total_records, 'records_filt_length': records_filt_length, 'records_filt_eukaryote': records_filt_eukaryote, 'records_filt': records_filt, 'fraction': fraction}, ignore_index=True)
# df_info.to_csv('results/06_metagenomicORFs/summary_ORF_filtering.csv', index=False)
# plot summary
df_info_plt = pd.DataFrame()
df_info_plt = df_info.melt(id_vars=['sample'], value_vars=['total_records', 'records_filt_length', 'records_filt_eukaryote', 'records_filt'], var_name='type', value_name='count')

plt.figure(figsize=(10, 15))
sns.barplot(y='sample', x='count', hue='type', data=df_info_plt[df_info_plt['sample'].str.contains('M')])
plt.xticks(rotation=90)
plt.savefig('results/figures/summary_ORF_filtering-M.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(10, 15))
sns.barplot(y='sample', x='count', hue='type', data=df_info_plt[df_info_plt['sample'].str.contains('C')])
plt.xticks(rotation=90)
plt.savefig('results/figures/summary_ORF_filtering-C.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(10, 15))
sns.barplot(y='sample', x='count', hue='type', data=df_info_plt[df_info_plt['sample'].str.contains('D')])
plt.xticks(rotation=90)
plt.savefig('results/figures/summary_ORF_filtering-D.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(10, 15))
sns.barplot(y='sample', x='count', hue='type', data=df_info_plt[df_info_plt['sample'].str.contains('F')])
plt.xticks(rotation=90)
plt.savefig('results/figures/summary_ORF_filtering-F.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(10, 15))
sns.barplot(y='sample', x='count', hue='type', data=df_info_plt[df_info_plt['sample'].str.contains('A')])
plt.xticks(rotation=90)
plt.savefig('results/figures/summary_ORF_filtering-A.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(10, 25))
sns.barplot(y='sample', x='count', hue='type', data=df_info_plt[df_info_plt['sample'].str.contains('1-')], order = df_info_plt[df_info_plt['sample'].str.contains('1-')].sort_values('sample')['sample'])
plt.xticks(rotation=90)
plt.savefig('results/figures/summary_ORF_filtering-1.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(10, 25))
sns.barplot(y='sample', x='count', hue='type', data=df_info_plt[df_info_plt['sample'].str.contains('6-')], order = df_info_plt[df_info_plt['sample'].str.contains('6-')].sort_values('sample')['sample'])
plt.xticks(rotation=90)
plt.savefig('results/figures/summary_ORF_filtering-6.png', dpi=300, bbox_inches='tight')