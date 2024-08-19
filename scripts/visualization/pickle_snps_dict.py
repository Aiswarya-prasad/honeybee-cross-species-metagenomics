import inStrain
import inStrain.SNVprofile
import pandas as pd
import numpy as np
from collections import defaultdict 

import glob
import os
import pickle
import gzip
import sys

all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')
handmade_spec_names = pd.read_csv('results/figures/handmade_species_names.csv')
handmade_spec_names['magotu'] = handmade_spec_names['cluster'].astype(str)
all_scaffolds_info = all_scaffolds_info.merge(handmade_spec_names, on='magotu', how='left')
mag_genus_dict = all_scaffolds_info.set_index('mag').to_dict()['genus']
mag_species_dict = all_scaffolds_info.set_index('mag').to_dict()['MAG_species_name_final']
scaffold_mag_dict = all_scaffolds_info.set_index('scaffold').to_dict()['mag']

# only includes detected species entries - 0.5 breadth cutoff
df_species_ab = pd.read_csv('results/figures/species_abundance_table.csv')

def species_in_sample(species_a, sample_a):
    '''
    check is the species is detected in that sample
    return true or false make sure df_species_ab is read
    '''
    species_detected = set(df_species_ab[df_species_ab['sample'] == sample_a]['MAG_species_name_final'])
    if species_a in species_detected:
        return True
    else:
        return False

samples = [x.split('/')[-2] for x in glob.glob('results/10_instrain/02_instrain_profile/*/')]
samples_to_exclude = ['F4-5', 'F5-1', 'M6-2', 'D9-5', 'F7-5']
samples_in_my = [x for x in samples if 'Dr' not in x and 'Gr' not in x and 'Am' not in x and 'Ac' not in x and x not in samples_to_exclude]

samples = [x.split('/')[-2] for x in glob.glob('results/10_instrain/02_instrain_profile/*/')]
# os.makedirs('results/figures/10-instrain_SNP_summaries', exist_ok=True)

species_list = set([x for x in mag_species_dict.values() if not pd.isna(x)])

for species in species_list:
    os.makedirs('results/figures/10-instrain_SNP_summaries/' + species, exist_ok=True)
    
df_out = pd.DataFrame(columns=['sample', 'Positions', 'SNV', 'SNS', 'con_SNV', 'pop_SNV'])
dfs_out = {x: df_out.copy() for x in species_list}
for i, sample in enumerate(samples_in_my):
    print(f'Processing sample {i+1}/{len(samples_in_my)}', end='\r')
    IS = inStrain.SNVprofile.SNVprofile('results/10_instrain/02_instrain_profile/' + sample)
    snv_table = IS.get_nonredundant_snv_table()
    snv_table['mag'] = snv_table['scaffold'].map(scaffold_mag_dict)
    snv_table['species'] = snv_table['mag'].map(mag_species_dict)
    for spec in species_list:
        if species_in_sample(spec, sample):
            total_positions = 0
            total_SNV = 0
            total_SNS = 0
            total_con_SNV = 0
            total_pop_SNV = 0
            snv_table_spec = snv_table[snv_table['species'] == spec]
            scaffolds = set(snv_table_spec['scaffold'])
            lengths = [int(x.split('_length_')[1].split('_cov')[0]) for x in scaffolds]
            total_positions = sum(lengths)
            total_SNV = snv_table_spec[snv_table_spec['class'] == 'SNV'].shape[0]
            total_SNS = snv_table_spec[snv_table_spec['class'] == 'SNS'].shape[0]
            total_con_SNV = snv_table_spec[snv_table_spec['class'] == 'con_SNV'].shape[0]
            total_pop_SNV = snv_table_spec[snv_table_spec['class'] == 'pop_SNV'].shape[0]
            dfs_out[spec] = dfs_out[spec]._append({'sample': sample, 'Positions': total_positions, 'SNV': total_SNV, 'SNS': total_SNS, 'con_SNV': total_con_SNV, 'pop_SNV': total_pop_SNV}, ignore_index=True)

# snv_table['class'].value_counts()

for spec in species_list:
    dfs_out[spec].to_csv('results/figures/10-instrain_SNP_summaries/' + spec + '/SNP_summary.tsv', sep='\t', index=False)

pickle.dump(dfs_out, open('results/figures/10-instrain_SNP_summaries/dfs_out.pkl', 'wb'))


# samples_read = []
# # # big dict of everything
# # if os.path.exists(f'results/figures/10-instrain_SNP_summaries/mega_snps_dict.pkl'):
# #     mega_snps_dict = pickle.load(open('results/figures/10-instrain_SNP_summaries/mega_snps_dict.pkl', 'rb'))
# # else:
# #     mega_snps_dict = {sample: {} for sample in samples_in_my}
# mega_snps_dict = {sample: {} for sample in samples_in_my}
# for i, sample in enumerate(samples_in_my):
#     if sample in samples_read:
#         continue
#     print(f'reading {i}/{len(samples_in_my)}')
#     zipped = False
#     snv_info_file = f'results/10_instrain/02_instrain_profile/{sample}/output/{sample}_SNVs.tsv'
#     if not os.path.exists(snv_info_file):
#         if os.path.exists(f'{snv_info_file}.gz'):
#             zipped = True
#     if zipped:
#         f = gzip.open(f'{snv_info_file}.gz', 'rt')
#     else:
#         f = open(snv_info_file, 'r')
#     for n, line in enumerate(f):
#         print(f'{sample} line #{n}', end = '\r')
#         if line.startswith(f'scaffold\t'):
#             continue
#         scaffold = line.split('\t')[0]
#         position = line.split('\t')[1]
#         allele_count = int(line.split('\t')[3])
#         if allele_count > 1:
#             spec_of_scaffold = mag_species_dict[scaffold_mag_dict[scaffold]]
#             if species_in_sample(spec_of_scaffold, sample):
#                 if spec_of_scaffold in mega_snps_dict[sample]:
#                     mega_snps_dict[sample][spec_of_scaffold].add(f'{scaffold}@{position}')
#                 else:
#                     mega_snps_dict[sample][spec_of_scaffold] = set([f'{scaffold}@{position}'])
#     print(f'size of mega dict is {sys.getsizeof(mega_snps_dict)/1000000}MB')
#     samples_read.append(sample)
#     pickle.dump(mega_snps_dict, open('results/figures/10-instrain_SNP_summaries/mega_snps_dict.pkl', 'wb'))


samples_read = []
# # big dict of everything
# if os.path.exists(f'results/figures/10-instrain_SNP_summaries/mega_snps_dict.pkl'):
#     mega_snps_dict = pickle.load(open('results/figures/10-instrain_SNP_summaries/mega_snps_dict.pkl', 'rb'))
# else:
#     mega_snps_dict = {sample: {} for sample in samples_in_my}
mega_snps_dict = {sample: {} for sample in samples}
for i, sample in enumerate(samples):
    if sample in samples_read:
        continue
    print(f'reading {i}/{len(samples)}')
    zipped = False
    snv_info_file = f'results/10_instrain/02_instrain_profile/{sample}/output/{sample}_SNVs.tsv'
    if not os.path.exists(snv_info_file):
        if os.path.exists(f'{snv_info_file}.gz'):
            zipped = True
    if zipped:
        f = gzip.open(f'{snv_info_file}.gz', 'rt')
    else:
        f = open(snv_info_file, 'r')
    for n, line in enumerate(f):
        print(f'{sample} line #{n}', end = '\r')
        if line.startswith(f'scaffold\t'):
            continue
        scaffold = line.split('\t')[0]
        position = line.split('\t')[1]
        allele_count = int(line.split('\t')[3])
        if allele_count > 1:
            spec_of_scaffold = mag_species_dict[scaffold_mag_dict[scaffold]]
            if species_in_sample(spec_of_scaffold, sample):
                if spec_of_scaffold in mega_snps_dict[sample]:
                    mega_snps_dict[sample][spec_of_scaffold].add(f'{scaffold}@{position}')
                else:
                    mega_snps_dict[sample][spec_of_scaffold] = set([f'{scaffold}@{position}'])
    print(f'size of mega dict is {sys.getsizeof(mega_snps_dict)/1000000}MB')
    samples_read.append(sample)
    pickle.dump(mega_snps_dict, open('results/figures/10-instrain_SNP_summaries/mega_snps_dict_all.pkl', 'wb'))