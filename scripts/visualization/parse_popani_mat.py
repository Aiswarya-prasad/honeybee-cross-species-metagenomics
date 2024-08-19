import argparse
import os
import pandas as pd
import numpy as np
from pprint import pprint

'''
called from R as:
system(paste0("source ~/.bashrc && conda activate 20230313_scripts_env && python3 scripts/visualization/parse_popani_mat.py -i results/figures/instrain_comp_mat/", spec_formatted, "_instrain_comp_info.csv"))

it reads in the instrain_comp_info.csv file which contains the popANI and conANI values
for each pair of samples. It then parses this into a square matrix with the popANI values
and writes it to a new csv file.

popANI is a measure of how similar two genomes are. It is calculated by aligning the two
genomes and counting the number of nucleotides that are identical. The input is already
expected to be filtered to only include pairs of samples that have enough overlap to
calculate popANI. The output is a square matrix with dissimilarity values (1 - popANI)
for each pair of samples. If there is no popANI value for a pair of samples, the value
is set to None.

'''

def parse_args():
    parser = argparse.ArgumentParser(description='Parse into a filled square matrix')
    parser.add_argument('-i', help='Scaffolds file containing all scaffolds for each sample')
    return parser.parse_args()

args = parse_args()

# input_file = 'results/figures/instrain_comp_mat/Bombilactobacillus mellis_instrain_comp_info.csv'
input_file = args.i
output_file = input_file.split('.csv')[0] + '_parsed.csv'

df_info = pd.read_csv(input_file)
names = list(set(df_info['name1']).union(set(df_info['name1'])))
popani_dict = df_info.groupby(['name1', 'name2']).agg(popANI = ('popANI', 'min')).to_dict()['popANI']
conani_dict = df_info.groupby(['name1', 'name2']).agg(conANI = ('conANI', 'min')).to_dict()['conANI']

mat = np.zeros((len(names), len(names)))
# make a square matrix using these names
for i, name_i in enumerate(names):
    for j, name_j in enumerate(names):
        if i == j:
            mat[i,j] = 0
        else:
            popani = None
            try:
                popani = popani_dict[(name_i, name_j)]
            except KeyError:
                try:
                    popani = popani_dict[(name_j, name_i)]
                except KeyError:
                    popani = None
            if popani is None:
                mat[i,j] = mat[j,i] = None
            else:
                mat[i,j] = mat[j,i] = 1 - float(popani)


pd.DataFrame(mat, index=names, columns=names).to_csv(output_file, index=True, header=True, sep=',')
print('done')