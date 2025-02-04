import os
import shutil
import sys
import io
import json
import gzip
import math
import json
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
# import statsmodels.sandbox.stats.multicomp as sm
# from enrichm.draw_plots import Plot
# from enrichm.databases import Databases
# from enrichm.module_description_parser import ModuleDescription
# from enrichm.parser import Parser, ParseAnnotate
# from enrichm.writer import Writer
# from enrichm.synteny_searcher import SyntenySearcher

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
        
samples_to_exclude = ["F4-5", "F5-1", "M6-2", "D9-5", "F7-5"]

samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]
samples = [x for x in samples if x not in samples_to_exclude]

# minpath tells you which pathways are found (the families is not a measure of completeness! That is inferred from microbeannotator..)
'''
The goal of using minpath here is to infer a table with samples as rows and pathways as columns with counts of each pathway
inferred as the minimum number of hits of the KOs in the pathway for the sample
First create a directory for each sample and write the KOs to a file
also include a file for each genus and each species only for the samples that have that genus or species
the minpath input needs to be of the form
dir structure:
    results/figures/summarize/minpath_inputs
        by_sample
            sample1_kos.csv
            sample2_kos.csv
            ...
        by_genus
            sample1_genus1_kos.csv
            sample1_genus2_kos.csv
            ...
            sample2_genus1_kos.csv
            ...
        by_species
            sample1_species1_kos.csv
            sample1_species2_kos.csv
            ...
            sample2_species1_kos.csv
            ...
within each file there are 2 columns (no header) with the first column being the gene id 
and the second column being the KO id (this can be done with EC number for metacyc if needed?)
'''
os.makedirs('results/figures/08-summarize_functions/minpath_KO_collections', exist_ok=True)s

# get the gene counts from gff quantifier output
for i, sample in enumerate(samples):
    print(f'working on {sample}')
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    # also make an accompanying file with the gene id and the read coverage of that gene that can later be
    # used to infer the coverage of the pathway by using minimum coverage of the KOs in the pathway
    # and KO coverage is the summed coverage of the genes that are annotated with that KO
    print(f'collecting genes in {sample}')
    genes_df = pd.read_csv(f'results/08_gene_content/01_profiling_bowtie2/{sample}_gene_coverage.txt', sep = '\t', header = None)
    genes_df.columns = ['scaffold', 'start', 'end', 'id', 'num_reads', 'num_bases', 'length', 'cov_fraction']
    genes_df['gene'] = [f'{sample}_']*len(genes_df['scaffold']) + genes_df['scaffold'].apply(lambda x: str(x)) + ['_']*len(genes_df['id']) + genes_df['id'].apply(lambda x: x.split('_')[1])
    # only consider a gene detected if there are >5 reads (1000/150) and at least 90% of the bases in the gene are covered
    genes_df = genes_df.loc[(genes_df['num_reads'] > 5) & (genes_df['cov_fraction'] > 0.90)]
    # now genes df should match what df_detected_genes_info has
    # next to it, write {sample}_df_detected_genes_cov.csv file
    genes_df = genes_df[['gene', 'start', 'end', 'num_reads', 'cov_fraction']]
    # bedtools coverage might have double counted some reads becuase they overlap with multiple genes - we tolerate this
    # rpkm = reads / ((gene_length/1000) * total_reads * 1000000)
    # cpm = reads / total_reads * 1000000
    total_reads = genes_df['num_reads'].sum()
    genes_df['rpkm'] = genes_df['num_reads'] / ( ((genes_df['end'] - genes_df['start'])/1000) * (total_reads/1000000) )
    genes_df['cpm'] = genes_df['num_reads'] / total_reads * 1000000
    genes_df.to_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv', index=False)

os.makedirs('results/figures/08-summarize_functions/ko_counts', exist_ok=True)
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, cazy and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus']]
    df_joined = df_info.merge(df_covs, on='gene')
    # get summed count per cazy
    df_joined = df_joined.groupby(['ko', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined.to_csv(f'results/figures/08-summarize_functions/ko_counts/{sample}_ko_counts.csv', index=False)

# write tables to be read for ko cumulative curve
genus_names = []
kos_in_sample = {}
for sample in samples:
    df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    kos_in_sample[sample] = set(df['ko'].unique())
    genus_names = genus_names + list(df['genus'].unique())
genus_names = list(set(genus_names))
pickle.dump(kos_in_sample, open('results/figures/08-summarize_functions/kos_in_sample.p', 'wb'))
pickle.dump(genus_names, open('results/figures/08-summarize_functions/genus_names.p', 'wb'))

os.makedirs('results/figures/08-summarize_functions/cumulative_curves', exist_ok=True)
iterations = 20
with open('results/figures/08-summarize_functions/cumulative_curves/ko_curves.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tkos\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(iterations):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/11 and sample size {j}/{len(samples_host)}', end = '\r')
                cum_kos = 0
                samples_iter = r.sample(samples_host, sample_size)
                kos_in_sample_iter = set()
                for sample in samples_iter:
                    kos_in_sample_iter.update(kos_in_sample[sample])
                cum_kos = len(kos_in_sample_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_kos}\n')

kos_in_sample_genus = {}
for genus in genus_names:
    kos_in_sample_genus[genus] = {}
    for sample in samples:
        df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        df = df[df['genus'] == genus]
        kos_in_sample_genus[genus][sample] = set(df['ko'].unique())
pickle.dump(kos_in_sample_genus, open('results/figures/08-summarize_functions/kos_in_sample_genus.p', 'wb'))

for genus in genus_names:
    print(f'working on {genus}')
    with open(f'results/figures/08-summarize_functions/cumulative_curves/ko_curves_{genus}.tsv', 'w+') as out_fh:
        out_fh.write(f'host\titeration\tsize\tkos\n')
        for i, host in enumerate(host_species):
            print(f'working on {host} {i}/5 ')
            samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
            for iteration in range(iterations):
                for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                    print(f'working om iteration {iteration}/11 and sample size {j}/{len(samples_host)}', end = '\r')
                    cum_kos = 0
                    samples_iter = r.sample(samples_host, sample_size)
                    kos_in_sample_iter = set()
                    for sample in samples_iter:
                        kos_in_sample_iter.update(kos_in_sample_genus[genus][sample])
                    cum_kos = len(kos_in_sample_iter)
                    n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_kos}\n')

# cpm and rpkm normalize for total reads but the pm part makes it per million reads and in essense results in 
# a relative value which will be higher for samples with fewer genes detected in total
os.makedirs('results/figures/08-summarize_functions/cazyme_counts', exist_ok=True)
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, cazy and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'cazyme', 'mag', 'species', 'genus']]
    df_joined = df_info.merge(df_covs, on='gene')
    # get summed count per cazy
    df_joined = df_joined.groupby(['cazyme', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined.to_csv(f'results/figures/08-summarize_functions/cazyme_counts/{sample}_cazyme_counts.csv', index=False)

os.makedirs('results/figures/08-summarize_functions/minpath_KO_collections/by_sample/', exist_ok=True)
os.makedirs('results/figures/08-summarize_functions/minpath_KO_collections/by_genus/', exist_ok=True)
os.makedirs('results/figures/08-summarize_functions/minpath_KO_collections/by_species/', exist_ok=True)

for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    # write the list of kos withouth the header with first column as gene id and second column as ko as tsv
    df_sample[['gene', 'ko']].to_csv(f'results/figures/08-summarize_functions/minpath_KO_collections/by_sample/{sample}_kos.tsv', sep='\t', index=False, header=False)
    for genus in df_sample['genus'].unique():
        os.makedirs(f'results/figures/08-summarize_functions/minpath_KO_collections/by_genus/{genus}', exist_ok=True)
        df_genus = df_sample[df_sample['genus'] == genus]
        # only write the list of kos withouth the header
        df_genus[['gene', 'ko']].to_csv(f'results/figures/08-summarize_functions/minpath_KO_collections/by_genus/{genus}/{sample}_kos.tsv', sep='\t', index=False, header=False)
        print(f'printed {genus} kos', end='\r')
    for species in df_sample['species'].unique():
        os.makedirs(f'results/figures/08-summarize_functions/minpath_KO_collections/by_species/{species}', exist_ok=True)
        df_species = df_sample[df_sample['species'] == species]
        # only write the list of kos withouth the header
        df_species[['gene', 'ko']].to_csv(f'results/figures/08-summarize_functions/minpath_KO_collections/by_species/{species}/{sample}_kos.tsv', sep='\t', index=False, header=False)
        print(f'printed {species} kos', end='\r')

# run mipath for the by_sample files
# python ../MinPath.py -ko ../../../results/figures/visualize_temp/KOs.csv -report KOs_test_gil.ko.minpath -details KOs_test_gil.ko.minpath.details
# Sinteractive -m50G
# conda activate 20230313_scripts_env
# python3
os.makedirs('results/figures/08-summarize_functions/minpath_outputs', exist_ok=True)
os.makedirs('results/figures/08-summarize_functions/minpath_outputs/by_sample', exist_ok=True)
for sample in samples:
    print(f'running minpath for {sample}')
    if os.path.exists(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.minpath.details'):
        print(f'already done with {sample}')
        continue
    command = f'python scripts/MinPath/MinPath.py -ko \
        results/figures/08-summarize_functions/minpath_KO_collections/by_sample/{sample}_kos.tsv \
        -report results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.minpath \
        -details results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.minpath.details'
    os.system(command)
    print(f'done with {sample}')

# get list of kos per mag and write to file per mag
df_mag_kos = pd.read_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_abs_all_kos.csv')
df_mag_kos = df_mag_kos.set_index('mag')
for mag in df_mag_kos.index:
    i = 0
    kos = df_mag_kos.loc[mag].dropna().index
    # only keep kos for whom rpkm > 0
    kos = [x for x in kos if df_mag_kos.loc[mag][x] > 0]
    with open(f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/{mag}_kos.tsv', 'w+') as out_fh:
        for ko in kos:
            out_fh.write(f'{mag}_{i}\t{ko}\n')
            i += 1
os.makedirs('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag', exist_ok=True)


all_mags = [os.path.basename(x).split('_kos.tsv')[0] for x in glob.glob('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/*')]

for n, mag in enumerate(all_mags):
    print(f'working on {mag} {n+1}/{len(all_mags)}           ', end = '\r')
    if os.path.exists(f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath.details'):
        continue
    command = f'python3 scripts/MinPath/MinPath.py -ko \
        additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/{mag}_kos.tsv \
        -report additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath \
        -details additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath.details'
    os.system(command)
    # print(f'done with {mag}')

os.makedirs(f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag_renamed/', exist_ok=True)
df_mag_details = pd.read_csv('results/figures/magOTU_handmade_names.csv')
# rename the MAGs with species names
for mag in all_mags:
    species_name = df_mag_details[df_mag_details['ID'] == mag]['MAG_species_name_final'].values[0]+f'--{mag}'
    # replace spaces with underscore
    species_name = species_name.replace(' ', '_')
    os.system(f'cp additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath.details additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag_renamed/{species_name}.minpath.details')
    os.system(f'cp additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag_renamed/{species_name}.minpath')

for mag in all_mags:
    species_name = df_mag_details[df_mag_details['ID'] == mag]['MAG_species_name_final'].values[0]+f'--{mag}'
    # replace spaces with underscore
    species_name = species_name.replace(' ', '_')
    minpath_file = f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag_renamed/{species_name}.minpath'
    minpath_parsed = f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag_renamed/{species_name}.minpath.parsed'
    completed_mags = set([x.split('by_mag/')[1].split('_kos')[0] for x in pd.read_csv('additional_analyses/results/functional_comparison/08-summarize_functions/completed_mags', sep = " ", header = None)[5].values])
    if not os.path.exists(minpath_parsed) or mag not in completed_mags:
        print(f'{species_name} {mag} not completed')
        os.remove(minpath_parsed)
        continue
    else:
        with open(minpath_file, 'r') as in_fh:
            with open(minpath_parsed, 'w+') as out_fh:
                for line in in_fh:
                    line = line.replace('  ', '\t')
                    line = line.replace('naive ', 'naive\t')
                    line = line.replace('minpath ', 'minpath\t')
                    out_fh.write(line)

with open('additional_analyses/results/functional_comparison/08-summarize_functions/all_mag_names.csv', 'w+') as out_fh:
    for mag in all_mags:
        species_name = df_mag_details[df_mag_details['ID'] == mag]['MAG_species_name_final'].values[0]
        species_name_no_space = species_name.replace(' ', '_')
        file_name = species_name_no_space+f'--{mag}'
        genus = df_mag_details[df_mag_details['ID'] == mag]['MAG_species_name_final'].values[0].split(' ')[0]
        # print(f'{mag},{species_name},{file_name}\n')
        out_fh.write(f'{mag},{species_name},{file_name},{genus}\n')

# df_pathway_presence = pd.DataFrame()
# df_pathway_fam_found = pd.DataFrame()
# # summarize into pathway completeness or presence absence at first
# for n, mag in enumerate(all_mags):
#     mag = 'M7-5_7'
#     minpath_out = f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath'
#     with open(minpath_out, 'r') as in_fh:
#         for line in in_fh:
#             pathway_id = line.split('  ')[0].split(' ')[1]
#             pathway_name = line.split('#')[1].strip()

# compare input ad output of "C4-5_2", "D5-1_10", "D3-1_33"

pd.read_csv('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/C4-5_2_kos.tsv', sep='\t', header=None)
pd.read_csv('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/D5-1_10_kos.tsv', sep='\t', header=None)
pd.read_csv('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/D3-1_33_kos.tsv', sep='\t', header=None)

col_1_df = pd.read_csv('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/C4-5_2_kos.minpath', sep='  ', header=None)
col_2_df = pd.read_csv('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/D5-1_10_kos.minpath', sep='  ', header=None)
col_3_df = pd.read_csv('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/D3-1_33_kos.minpath', sep='  ', header=None)



# mag = 'C4-5_2'
mag = 'D5-1_10'

command = f'python3 scripts/MinPath/MinPath.py -ko \
        additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/{mag}_kos.tsv \
        -report additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath \
        -details additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath.details'
os.system(command)

# parse the minpath main output and details file
# the main output file has the pathway id and the pathway name
# the details file has the pathway id and the KO id and the number of hits and only includes pathways minpath identified as present
# details file looks like:
'''
path 00130 fam0 59 fam-found 16 # Ubiquinone and other terpenoid-quinone biosynthesis
   K00355 hits 5 # NQO1; NAD(P)H dehydrogenase (quinone) [EC:1.6.5.2]
   K00568 hits 5 # ubiG; 2-polyprenyl-6-hydroxyphenyl methylase / 3-demethylubiquinone-9 3-methyltransferase [EC:2.1.1.222 2.1.1.64]
   K01661 hits 2 # menB; naphthoate synthase [EC:4.1.3.36]
   K01911 hits 1 # menE; o-succinylbenzoate---CoA ligase [EC:6.2.1.26]
'''

for i, sample in enumerate(samples):
    # if os.path.exists(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.path.rpkm'):
    #     print(sample)
    #     continue
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko']]
    df_joined = df_info.merge(df_covs, on='gene')
    # get summed count per ko
    df_joined = df_joined.groupby(['ko']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    with open(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.minpath.details', 'r') as in_fh:
        path_ko_dict = {}
        for line in in_fh:
            if line.startswith('path'):
                # print(f'Working on: {line}')
                path_id = line.split()[1]
                path_name = line.split('#')[1].strip()
                path_ko_dict[f'path {path_id}'] = {'name': path_name, 'kos': []}
            if line.strip().startswith('K'):
                ko_id = line.split()[0]
                path_ko_dict[f'path {path_id}']['kos'].append(ko_id)
        pathway_counts = {x: {'name': path_ko_dict[x]['name'], 'rpkm': 0, 'cpm': 0} for x in path_ko_dict.keys()}
        for path_iter in path_ko_dict.keys():
            counts = []
            rpkms = []
            for ko_id in path_ko_dict[path_iter]['kos']:                
                # get number as cpm or rpkm
                count = df_joined[df_joined['ko'] == ko_id]['cpm'].values[0]
                rpkm = df_joined[df_joined['ko'] == ko_id]['rpkm'].values[0]
                counts.append(count)
                rpkms.append(rpkm)
            pathway_counts[path_iter]['cpm'] = min(counts)
            pathway_counts[path_iter]['rpkm'] = min(rpkms)
    with open(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.path.rpkm', 'w+') as out_fh:
        success = out_fh.write('pathway_id\tpathway_name\trpkm\n')
        for path in pathway_counts.keys():
            success = out_fh.write(f'{path}\t{pathway_counts[path]["name"]}\t{pathway_counts[path]["rpkm"]}\n')
    with open(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.path.cpm', 'w+') as out_fh:
        success = out_fh.write('pathway_id\tpathway_name\tcpm\n')
        for path in pathway_counts.keys():
            success = out_fh.write(f'{path}\t{pathway_counts[path]["name"]}\t{pathway_counts[path]["cpm"]}\n')

for i, sample in enumerate(samples):
    if os.path.exists(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.path.rpkm'):
        print(sample)
        continue
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko']]
    df_joined = df_info.merge(df_covs, on='gene')
    # get summed count per ko
    df_joined = df_joined.groupby(['ko']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    with open(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.minpath.details', 'r') as in_fh:
        path_ko_dict = {}
        for line in in_fh:
            if line.startswith('path'):
                # print(f'Working on: {line}')
                path_id = line.split()[1]
                path_name = line.split('#')[1].strip()
                path_ko_dict[f'path {path_id}'] = {'name': path_name, 'kos': []}
            if line.strip().startswith('K'):
                ko_id = line.split()[0]
                path_ko_dict[f'path {path_id}']['kos'].append(ko_id)
        pathway_counts = {x: {'name': path_ko_dict[x]['name'], 'rpkm': 0, 'cpm': 0} for x in path_ko_dict.keys()}
        for path_iter in path_ko_dict.keys():
            counts = []
            rpkms = []
            for ko_id in path_ko_dict[path_iter]['kos']:                
                # get number as cpm or rpkm
                count = df_joined[df_joined['ko'] == ko_id]['cpm'].values[0]
                rpkm = df_joined[df_joined['ko'] == ko_id]['rpkm'].values[0]
                counts.append(count)
                rpkms.append(rpkm)
            pathway_counts[path_iter]['cpm'] = min(counts)
            pathway_counts[path_iter]['rpkm'] = min(rpkms)
    with open(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.path.rpkm', 'w+') as out_fh:
        success = out_fh.write('pathway_id\tpathway_name\trpkm\n')
        for path in pathway_counts.keys():
            success = out_fh.write(f'{path}\t{pathway_counts[path]["name"]}\t{pathway_counts[path]["rpkm"]}\n')
    with open(f'results/figures/08-summarize_functions/minpath_outputs/by_sample/{sample}_kos.path.cpm', 'w+') as out_fh:
        success = out_fh.write('pathway_id\tpathway_name\tcpm\n')
        for path in pathway_counts.keys():
            success = out_fh.write(f'{path}\t{pathway_counts[path]["name"]}\t{pathway_counts[path]["cpm"]}\n')

os.makedirs('results/figures/08-summarize_functions/minpath_outputs/by_genus', exist_ok=True)
# go through each genus and do the same if there are >5 samples with that genus
# get genus names from the by_genus dir
genus_names = [x.split('/')[-1] for x in glob.glob('results/figures/08-summarize_functions/minpath_KO_collections/by_genus/*')]
for genus in genus_names:
    # skip if genus is nan
    if 'nan' in genus:
        continue
    if len(glob.glob(f'results/figures/08-summarize_functions/minpath_KO_collections/by_genus/{genus}/*')) < 5:
        continue
    print(f'working on {genus}')
    samples = [x.split('/')[-1].split('_kos.tsv')[0] for x in glob.glob(f'results/figures/08-summarize_functions/minpath_KO_collections/by_genus/{genus}/*')]
    samples = [x for x in samples if x not in samples_to_exclude]
    os.makedirs(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}', exist_ok=True)
    for sample in samples:
        print(f'working on {sample}')
        if os.path.exists(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/{sample}_kos.minpath.details'):
            # print(f'already done with {sample}')
            # continue
            pass
        else:
            command = f'python scripts/MinPath/MinPath.py -ko \
                results/figures/08-summarize_functions/minpath_KO_collections/by_genus/{genus}/{sample}_kos.tsv \
                -report results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/{sample}_kos.minpath \
                -details results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/{sample}_kos.minpath.details'
            os.system(command)
            print(f'done with {sample}')

for genus in genus_names:
    if len(glob.glob(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/*')) < 5:
        continue
    samples = [x.split('/')[-1].split('_kos.minpath.details')[0] for x in glob.glob(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/*')]
    samples = [x for x in samples if x not in samples_to_exclude]
    for sample in samples:
        print(f'working on {sample}')
        df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
        df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        # join them by gene and get the columns gene, ko and count
        df_covs = df_covs[['gene', 'rpkm', 'cpm']]
        df_info = df_info[['gene', 'ko']]
        df_joined = df_info.merge(df_covs, on='gene')
        # get summed count per ko
        df_joined = df_joined.groupby(['ko']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
        if not os.path.exists(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/{sample}_kos.minpath.details'):
            print(f'no minpath details file for {sample}')
            continue
        with open(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/{sample}_kos.minpath.details', 'r') as in_fh:
            path_ko_dict = {}
            for line in in_fh:
                if line.startswith('path'):
                    # print(f'Working on: {line}')
                    path_id = line.split()[1]
                    path_name = line.split('#')[1].strip()
                    path_ko_dict[f'path {path_id}'] = {'name': path_name, 'kos': []}
                if line.strip().startswith('K'):
                    ko_id = line.split()[0]
                    path_ko_dict[f'path {path_id}']['kos'].append(ko_id)
            pathway_counts = {x: {'name': path_ko_dict[x]['name'], 'rpkm': 0, 'cpm': 0} for x in path_ko_dict.keys()}
            for path_iter in path_ko_dict.keys():
                counts = []
                rpkms = []
                for ko_id in path_ko_dict[path_iter]['kos']:
                    # get number as cpm or rpkm
                    count = df_joined[df_joined['ko'] == ko_id]['cpm'].values[0]
                    rpkm = df_joined[df_joined['ko'] == ko_id]['rpkm'].values[0]
                    counts.append(count)
                    rpkms.append(rpkm)
                pathway_counts[path_iter]['cpm'] = min(counts)
                pathway_counts[path_iter]['rpkm'] = min(rpkms)
        with open(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/{sample}_kos.path.rpkm', 'w+') as out_fh:
            success = out_fh.write('pathway_id\tpathway_name\trpkm\n')
            for path in pathway_counts.keys():
                success = out_fh.write(f'{path}\t{pathway_counts[path]["name"]}\t{pathway_counts[path]["rpkm"]}\n')
        with open(f'results/figures/08-summarize_functions/minpath_outputs/by_genus/{genus}/{sample}_kos.path.cpm', 'w+') as out_fh:
            success = out_fh.write('pathway_id\tpathway_name\tcpm\n')
            for path in pathway_counts.keys():
                success = out_fh.write(f'{path}\t{pathway_counts[path]["name"]}\t{pathway_counts[path]["cpm"]}\n')

'''
make a matrix to be used for PCoA of KOs of a certain family of functions (eg. carbon metabolism)
with bacterial species (and host species) information
each column is a KO, and each row is a MAG
values within are rpkm
'''
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
# try KOs_of_interest from DRAM
DRAM_modules = pd.read_csv('scripts/visualization/DRAM_KOs_list_from_excel.txt', sep = '\t', header = None)
DRAM_modules.columns = ['ko', 'enzyme', 'module', 'module_type']
# set(DRAM_modules['module_type'])
# DRAM_modules = DRAM_modules[DRAM_modules['module_type'] == 'central carbon']
DRAM_modules_cazy = DRAM_modules[DRAM_modules['module_type'] == 'CAZY']
cazymes_of_interest = DRAM_modules_cazy['ko'].unique()
KOs_of_interest = [x for x in DRAM_modules['ko'].unique() if x.startswith('K') and len(x) == 6]
# make a matrix with MAGs as rows and KOs as columns with rpkm within
# for each sample, get the rpkm of the KOs in the KOs_of_interest list
# and make a row with the MAG name as the index
# and the rpkm as the values
# then concatenate all of these into a matrix
# then add the species and genus information as columns
# then add the host species information as a column
# then save this matrix
os.makedirs('results/figures/08-summarize_functions/ko_matrix_by_function/', exist_ok=True)
df_collected_rpkms = pd.DataFrame()
for i, sample in enumerate(samples):
    print(f'working on {i}/{len(samples)}', end = '\r')
    df_covs = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv')
    df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # join them by gene and get the columns gene, ko and count
    df_covs = df_covs[['gene', 'rpkm', 'cpm']]
    df_info = df_info[['gene', 'ko', 'mag', 'species', 'genus', 'cazyme']]
    df_joined = df_info.merge(df_covs, on='gene')
    # select kos interest and cazymes of interest together
    df_joined = df_joined[df_joined['ko'].notnull()]
    # df_joined = df_joined[df_joined['ko'].isin(KOs_of_interest)]
    # df_joined = df_joined[df_joined['ko'].isin(KOs_of_interest) | df_joined['cazyme'].isin(cazymes_of_interest)]
    # in rows where cazyme is not null replace ko value by cazyme name
    # df_joined.loc[df_joined['cazyme'].notnull(), 'ko'] = df_joined.loc[df_joined['cazyme'].notnull(), 'cazyme']
    # df_joined = df_joined[df_joined['cazyme'].notnull()]
    # # rename cazyme to ko
    # df_joined['ko'] = df_joined['cazyme']
    # get summed count per ko
    df_joined = df_joined.groupby(['ko', 'mag', 'species', 'genus']).agg({'rpkm': 'sum', 'cpm': 'sum'}).reset_index()
    df_joined = df_joined.pivot(index='mag', columns='ko', values='rpkm').reset_index()
    # df_joined['sample'] = sample
    df_collected_rpkms = pd.concat([df_collected_rpkms, df_joined])
df_collected_rpkms = df_collected_rpkms.fillna(0)
# remove columns with all 0s
df_collected_rpkms = df_collected_rpkms.loc[:, (df_collected_rpkms != 0).any(axis=0)]
df_collected_rpkms.to_csv('results/figures/08-summarize_functions/ko_matrix_by_function/ko_matrix_all_kos.csv', index=False)

# parse information about BGCs from antismash output in results/08_gene_content/09_antismash/{sample} and update gene_info_tables by writing a 
# new file appeneded in name with antismash listing all the features as individual rows with scaffold start and end information as this would
# not correspond to a given id, coverage would have to be extracted from bam files if needed!


all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')
handmade_spec_names = pd.read_csv('results/figures/handmade_species_names.csv')
handmade_spec_names['magotu'] = handmade_spec_names['cluster'].astype(str)
all_scaffolds_info = all_scaffolds_info.merge(handmade_spec_names, on='magotu', how='left')

mag_genus_dict = all_scaffolds_info.set_index('mag').to_dict()['genus']
mag_species_dict = all_scaffolds_info.set_index('mag').to_dict()['MAG_species_name_final']

binned_scaffolds = set(all_scaffolds_info['scaffold'])
scaffold_mag_dict = all_scaffolds_info.set_index('scaffold').to_dict()['mag']

def mag_from_scaffold(scaffold):
    scaffold = scaffold.strip('.')
    if scaffold not in scaffold_mag_dict.keys():
        return 'not_binned'
    return scaffold_mag_dict[scaffold]

def species_of_mag(mag):
    if mag == 'not_binned':
        return 'not_binned'
    return mag_species_dict[mag]

scaffold = 'D1-1_NODE_37_length_157980_cov_41.581086'
species_of_mag(mag_from_scaffold(scaffold))

for sample in samples:
    if not os.path.exists(f'results/08_gene_content/09_antismash/{sample}/{sample}.gbk'):
        print(f'{sample} not found')
        continue
    with open(f'results/figures/08-summarize_functions/09-antismash_parsed/{sample}_regions_summary.txt', 'w+') as out_fh:
        out_fh.write('scaffold\tstart\tend\tmag\tspecies\tproduct\tedge\n')
        with open (f'results/08_gene_content/09_antismash/{sample}/{sample}.gbk', 'r') as in_fh:
            success = gbk_records = SeqIO.parse(in_fh, 'genbank')
            for i, record in enumerate(gbk_records):
                # if i == 4:
                #     break
                scaffold = record.annotations['structured_comment']['antiSMASH-Data']['Original ID']
                # get mag of scaffold and species of mag
                mag = mag_from_scaffold(scaffold)
                species = species_of_mag(mag)
                for feature in record.features:
                    if feature.type == 'region':
                        start = feature.location.start.position
                        end = feature.location.end.position
                        print(f'working on {scaffold} in {mag} {species} at {start} {end}')
                        products = feature.qualifiers['product']
                        if len(products) > 1:
                            for i, product in enumerate(products):
                                product = feature.qualifiers['product'][i]
                                edge = feature.qualifiers['contig_edge'][0]
                                string = f'{scaffold}\t{start}\t{end}\t{mag}\t{species}\t{product}\t{edge}'
                                out_fh.write(string + '\n')
                                # feature.qualifiers['rules'][i]
                        else:
                            product = feature.qualifiers['product'][0]
                            edge = feature.qualifiers['contig_edge'][0]
                            string = f'{scaffold}\t{start}\t{end}\t{mag}\t{species}\t{product}\t{edge}'
                            out_fh.write(string + '\n')
                            # feature.qualifiers['rules'][0]
                        print('----------------------------')




# features do not correspond exactly to genes in the gene_info_tables...
        # for feature in record.features:
        #     start = feature.location.start.position
        #     end = feature.location.end.position
        #     if (scaffold, start, end) in df_genes_dict.keys():
        #         gene = df_genes_dict[(scaffold, start, end)]['gene']
        #         print(f'gene on {scaffold} found in df_genes_dict')
        #     else:
        #         # print(f'gene on {scaffold} not found in df_genes_dict')
        #         continue

# os.makedirs('results/figures/08-summarize_functions/minpath_outputs/by_species', exist_ok=True)
# # go through each species and do the same if there are >5 samples with that species
# # get species names from the by_species dir
# species_names = [x.split('/')[-1] for x in glob.glob('results/figures/08-summarize_functions/minpath_KO_collections/by_species/*')]
# for species in species_names:


# make a list of KOs with the title as the "genome name"
# run with cluster both
# this genome name can be any kind of grouping that we want to compare
# for example, we can use the host species as the genome name and
# compare the same microbial species and / genus across different hosts


# prepare the dictionary to get kegg ortholog information
current_a = ''
current_b = ''
current_c = ''
current_d = ''
kos_discovered = {}
with open('results/figures/08-summarize_functions/kegg_info_dict.txt', 'w+') as f_out:
    header = 'ko\tko_uniq\tA\tB\tC\tD\n'
    success = f_out.write(header)
    with open('data/KEGG_info/ko00001.keg', 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('+') or line.startswith('!'):
                continue
            first_B = True
            if line.startswith('A'):
                current_a = ' '.join(line.strip().split(' ')[1:])
                print(current_a)
                first_B = True
            # B for some reason has a B followed by an empty line and then the next line is a B with the info
            if line.startswith('B'):
                if line == 'B\n':
                    continue
                line = line.strip()
                current_b = ' '.join(line.split('B  ')[1].split(' ')[1:])
            if line.startswith('C'):
                line = line.strip()
                current_c = ' '.join(line.split('C    ')[1].split(' ')[1:])
                # current_c_num = current_c.split(' [PATH:')[0]
            if line.startswith('D'):
                line = line.strip()
                ko = line.split('D      ')[1].split(' ')[0]
                current_d = ' '.join(line.split('D      ')[1].split(' ')[1:])
                if ko not in kos_discovered:
                    kos_discovered[ko] = 0
                kos_discovered[ko] += 1
                success = f_out.write(f'{ko}\t{ko}_{kos_discovered[ko]}\t{current_a}\t{current_b}\t{current_c}\t{current_d}\n')

# with open('results/figures/08-summarize_functions/kegg_info_dict.txt', 'w+') as f:
#     header = 'ko\tA\tB\tC\n'
#     success = f.write(header)
#     for ko in kegg_info_dict.keys():
#         line = f'{ko}\t{kegg_info_dict[ko]["A"]}\t{kegg_info_dict[ko]["B"]}\t{kegg_info_dict[ko]["C"]}\n'
#         success = f.write(line)

# the above was how this was originally done and several heatmaps and other visualizations were made
# but the data was not used in the final manuscript
# now to address one of the reviewers comments we will make a heatmap and I noticed that there is
# an issue with the logic above so I have redone it here (above).




# form cumulate_.... script
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

# rpkm_matrix_abs_all_kos has duplicates of kos, names as K00001_1, K00001_2 etc.
# make a new rpkm matrix with the duplicates included and the values copied over
# avoid <stdin>:4: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`

df = pd.read_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_abs_all_kos.csv')
df_duped = pd.DataFrame()
kos_uniq_list = [x for x in pd.read_csv('results/figures/08-summarize_functions/kegg_info_dict.txt', sep = '\t')['ko_uniq'] if x.split('_')[0] in df.columns]
for ko_uniq in kos_uniq_list:
    print(f'working on {ko_uniq}    ', end = '\r')
    ko_non_uniq = ko_uniq.split('_')[0]
    rpkm_values = df[ko_non_uniq].values
    df_duped = pd.concat([df_duped, pd.DataFrame({ko_uniq: rpkm_values})], axis=1)
df_duped['mag'] = df['mag']
df_duped.to_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_abs_all_kos_duped.csv', index=False)