'''
pip3 installed instrain from github (31/10/23) into 20230313_scripts_env
###! Successfully installed asteval-0.9.31 biopython-1.74 future-0.18.3 h5py-3.10.0 inStrain-1.8.0 iniconfig-2.0.0 lmfit-1.2.2 pluggy-1.3.0 pysam-0.22.0 pytest-7.4.3 uncertainties-3.1.7
# redid later and had to edit conda by hand to fix issue with "alphabet" in biopython
'''

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

# # get variable sites for A2-2
# # IS = inStrain.SNVprofile.SNVprofile('results/10_instrain/02_instrain_profile/A2-2')
# snv_table = IS.get_nonredundant_snv_table()
# snv_table = snv_table[snv_table['class'] != 'SNS']
# snv_table['class'].value_counts()
# snv_table['position'].value_counts()
# # identify positions uniquely as a combination of scaffold and position - consider the format scaffold@position
# # call this position_id
# snv_table['position_id'] = snv_table['scaffold'] + '@' + snv_table['position'].astype(str)

all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')
handmade_spec_names = pd.read_csv('results/figures/handmade_species_names.csv')
handmade_spec_names['magotu'] = handmade_spec_names['cluster'].astype(str)
all_scaffolds_info = all_scaffolds_info.merge(handmade_spec_names, on='magotu', how='left')
mag_genus_dict = all_scaffolds_info.set_index('mag').to_dict()['genus']
mag_species_dict = all_scaffolds_info.set_index('mag').to_dict()['MAG_species_name_final']
# binned_scaffolds = set(all_scaffolds_info['scaffold'])
scaffold_mag_dict = all_scaffolds_info.set_index('scaffold').to_dict()['mag']
# phylo_metadata = pd.read_csv('results/11_phylogenies/phylo_genomes_metadata.tsv', sep='\t')


def host_of_sample(sample_name):
    if sample_name[0:2] in ['Gr', 'Dr', 'Am', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7']:
        return 'Apis mellifera'
    elif sample_name[0:2] in ['Ac', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']:
        return 'Apis cerana'
    elif sample_name[0:2] in ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9']:
        return 'Apis dorsata'
    elif sample_name[0:2] in ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9']:
        return 'Apis florea'
    elif sample_name[0:2] in ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']:
        return 'Apis andreniformis'
    else:
        return 'unknown'

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
    

'''
for each scaffold the recovered snv_table info about the variant sites
we can the count the number of variant sites per scaffold and the length
of the scaffold is in its name
>>> snv_table.columns
Index(['scaffold', 'position', 'ref_base', 'A', 'C', 'T', 'G', 'con_base',
       'var_base', 'allele_count', 'class', 'cryptic', 'position_coverage',
       'var_freq', 'con_freq', 'ref_freq'],
      dtype='object')
we make a file for each species containing for each sample one row an the
number of SNPs in that sample for that species as the sum across all the scaffolds
and the summed length of the scaffolds
also median and sd of the position_coverage values
'''

# run by the batch script scripts/visualization/pickle_snps_dict.sh
# code to make pickle it is in scripts/visualization/pickle_snps_dict.py
dfs_out = pickle.load(open('results/figures/10-instrain_SNP_summaries/dfs_out.pkl', 'rb'))

# prepare to make cumulative curves per species for SNPs
# we need to make a table with the following columns:
# iteration, number of samples, number of SNPs
# we need to first get the positions per scaffold that have SNPs in a given sample 
# store this in a set as scaffold_position and for subsequent samples add to it the 
# positions not already present in the set
# each time before adding the next sample put in the df the number of SNPs in the set

# run by the batch script scripts/visualization/pickle_snps_dict.sh
# code to make pickle it is in scripts/visualization/pickle_snps_dict.py
mega_snps_dict = pickle.load(open('results/figures/10-instrain_SNP_summaries/mega_snps_dict.pkl', 'rb'))

num_iters = 20
samples = [x.split('/')[-2] for x in glob.glob('results/10_instrain/02_instrain_profile/*/')]
# os.makedirs('results/figures/10-instrain_SNP_summaries', exist_ok=True)
samples_to_exclude = ['F4-5', 'F5-1', 'M6-2', 'D9-5', 'F7-5']
samples_in_my = [x for x in samples if 'Dr' not in x and 'Gr' not in x and 'Am' not in x and 'Ac' not in x and x not in samples_to_exclude]
species_list = set([x for x in mag_species_dict.values() if not pd.isna(x)])
for spec in species_list:
    with open(f'results/figures/10-instrain_SNP_summaries/{spec}/cumulative_SNPs.tsv', 'w') as f_out:
        success = f_out.write('iteration\thost\tnumber_of_samples\tnumber_of_SNPs\n')
        # f_out.write('iteration\thost\tnumber_of_samples\tnumber_of_SNPs\tmedian_position_coverage\tsamples_selected\n')
for iteration in range(num_iters):
    for s, spec in enumerate(species_list):
        samples_of_spec = [x for x in mega_snps_dict.keys() if spec in mega_snps_dict[x].keys()]
        print(f'iteration {iteration+1}/{num_iters} species {spec}')
        for host in ['Apis mellifera', 'Apis cerana', 'Apis dorsata', 'Apis florea', 'Apis andreniformis']:
            samples_of_host = [x for x in samples_of_spec if host_of_sample(x) == host]
            snps_seen = set()
            for sample_num in range(len(samples_of_host)):
                samples_selected = np.random.choice(samples_of_host, sample_num+1, replace=False)
                for sample in samples_selected:
                    snps_seen = snps_seen.union(mega_snps_dict[sample][spec])
                with open(f'results/figures/10-instrain_SNP_summaries/{spec}/cumulative_SNPs.tsv', 'a') as f_out:
                    success = f_out.write(f'{iteration+1}\t{host}\t{len(samples_selected)}\t{len(snps_seen)}\n')
                    print(f'spec: {s}/{len(species_list)}, iteration {iteration+1}/{num_iters}, host {host}, samples {len(samples_selected)}, number of SNPs: {len(snps_seen)}')

# also do this for all species and then by location instead of host
mega_snps_dict_all = pickle.load(open('results/figures/10-instrain_SNP_summaries/mega_snps_dict_all.pkl', 'rb'))
num_iters = 20
samples = [x.split('/')[-2] for x in glob.glob('results/10_instrain/02_instrain_profile/*/')]
# os.makedirs('results/figures/10-instrain_SNP_summaries', exist_ok=True)
samples_to_exclude = ['F4-5', 'F5-1', 'M6-2', 'D9-5', 'F7-5']
samples_in_my = [x for x in samples if 'Dr' not in x and 'Gr' not in x and 'Am' not in x and 'Ac' not in x and x not in samples_to_exclude]
species_list = set([x for x in mag_species_dict.values() if not pd.isna(x)])
for spec in species_list:
    with open(f'results/figures/10-instrain_SNP_summaries/{spec}/cumulative_SNPs_all.tsv', 'w') as f_out:
        success = f_out.write('iteration\thost\tnumber_of_samples\tnumber_of_SNPs\n')
        # f_out.write('iteration\thost\tnumber_of_samples\tnumber_of_SNPs\tmedian_position_coverage\tsamples_selected\n')
for iteration in range(num_iters):
    for s, spec in enumerate(species_list):
        samples_of_spec = [x for x in mega_snps_dict_all.keys() if spec in mega_snps_dict_all[x].keys()]
        print(f'iteration {iteration+1}/{num_iters} species {spec}')
        for host in ['Apis mellifera', 'Apis cerana', 'Apis dorsata', 'Apis florea', 'Apis andreniformis']:
            samples_of_host = [x for x in samples_of_spec if host_of_sample(x) == host]
            snps_seen = set()
            for sample_num in range(len(samples_of_host)):
                samples_selected = np.random.choice(samples_of_host, sample_num+1, replace=False)
                for sample in samples_selected:
                    snps_seen = snps_seen.union(mega_snps_dict_all[sample][spec])
                with open(f'results/figures/10-instrain_SNP_summaries/{spec}/cumulative_SNPs_all.tsv', 'a') as f_out:
                    success = f_out.write(f'{iteration+1}\t{host}\t{len(samples_selected)}\t{len(snps_seen)}\n')
                    print(f'spec: {s}/{len(species_list)}, iteration {iteration+1}/{num_iters}, host {host}, samples {len(samples_selected)}, number of SNPs: {len(snps_seen)}')
        



'''
IS = inStrain.SNVprofile.SNVprofile('results/10_instrain/02_instrain_profile/C2-2')
genes_SNP_count_df = IS.get('genes_SNP_count')
genes_SNP_count_df.columns

genes_SNP_count_df = genes_SNP_count_df[genes_SNP_count_df['mm'] <= 5]
genes_SNP_count_df = genes_SNP_count_df.sort_values('mm')\
            .drop_duplicates(subset=['gene'], keep='last')\
            .sort_index().drop(columns=['mm'])


# aim to get a cumulative curve of SNPs for each species across samples
x-axis: number of samples
y-axis: number of SNPs
For this we need a table with the following columns:
iteration, number of samples, number of SNPs
we need to first get the positions that have SNPs in a given sample and count them
then add to it the number of positions not already counted in previous samples
having snps
'''