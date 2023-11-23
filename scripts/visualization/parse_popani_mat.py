import argparse
import os
import pandas as pd
import numpy as np
from pprint import pprint

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