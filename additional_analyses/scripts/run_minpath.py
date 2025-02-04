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


# # get list of kos per mag and write to file per mag
# df_mag_kos = pd.read_csv('results/figures/08-summarize_functions/rpkm_matrix_by_function/rpkm_matrix_abs_all_kos.csv')
# df_mag_kos = df_mag_kos.set_index('mag')
# for mag in df_mag_kos.index:
#     i = 0
#     kos = df_mag_kos.loc[mag].dropna().index
#     # only keep kos for whom rpkm > 0
#     kos = [x for x in kos if df_mag_kos.loc[mag][x] > 0]
#     with open(f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/{mag}_kos.tsv', 'w+') as out_fh:
#         for ko in kos:
#             out_fh.write(f'{mag}_{i}\t{ko}\n')
#             i += 1

os.makedirs('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag', exist_ok=True)

all_mags = [os.path.basename(x).split('_kos.tsv')[0] for x in glob.glob('additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/*')]

for n, mag in enumerate(all_mags):
    # print(f'working on {mag} {n+1}/{len(all_mags)}           ', end = '\r')
    if os.path.exists(f'additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath.details'):
        continue
    command = f'python3 scripts/MinPath/MinPath.py -ko \
        additional_analyses/results/functional_comparison/08-summarize_functions/minpath_KO_collections/by_mag/{mag}_kos.tsv \
        -report additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath \
        -details additional_analyses/results/functional_comparison/08-summarize_functions/minpath_outputs/by_mag/{mag}_kos.minpath.details'
    os.system(command)


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
