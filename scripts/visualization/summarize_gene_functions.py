import os
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
import statsmodels.sandbox.stats.multicomp as sm
from enrichm.draw_plots import Plot
from enrichm.databases import Databases
from enrichm.module_description_parser import ModuleDescription
from enrichm.parser import Parser, ParseAnnotate
from enrichm.writer import Writer
from enrichm.synteny_searcher import SyntenySearcher

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
        
samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]

df = pd.read_csv(f'results/figures/all_sample_genes_df.csv')

# make a list of KOs with the title as the "genome name"
# run with cluster both
# this genome name can be any kind of grouping that we want to compare
# for example, we can use the host species as the genome name and
# compare the same microbial species and / genus across different hosts

# make a list of KOs with the title as the "genome name"
# run with cluster both
# this genome name can be any kind of grouping that we want to compare
os.makedirs('results/figures/visualize_temp/KO_collections', exist_ok=True)


# minpath tells you which pathways are found (the families is not a measure of completeness! That is inferred from microbeannotator)
# 
"python ../MinPath.py -ko ../../../results/figures/visualize_temp/KOs.csv -report KOs_test_gil.ko.minpath -details KOs_test_gil.ko.minpath.details"

'''
community-level: across samples, not subset by microbial species
'''
for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    os.makedirs(f'results/figures/visualize_temp/KO_collections/by_sample/', exist_ok=True)
    # only write the list of kos withouth the header
    df_sample['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_sample/{sample}_kos.csv', index=False, header=False)
    # run annotator
    print(f'printed {sample} kos')

os.makedirs(f'results/figures/visualize_temp/microbeannotator_out/', exist_ok=True)
input_list = ' '.join([f'results/figures/visualize_temp/KO_collections/by_sample/{sample}_kos.csv' for sample in samples])
print(f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_sample')

for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    for genus in df_sample['genus'].unique():
        df_genus = df_sample[df_sample['genus'] == genus]
        os.makedirs(f'results/figures/visualize_temp/KO_collections/by_genus/{genus}', exist_ok=True)
        # only write the list of kos withouth the header
        df_genus['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_genus/{genus}/{sample}_kos.csv', index=False, header=False)
        # run annotator
        print(f'printed {sample} kos')
    print(f'printed {sample} kos')

# this tells you the module completeness of each pathway
"python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i results/figures/visualize_temp/KOs.csv -p test_C3-2_snod"


'''
summarize functions of MAGs from samples
'''

# first get the list of KOs and name the file to have the handmade species name and MAG name (do not add any file extensions)