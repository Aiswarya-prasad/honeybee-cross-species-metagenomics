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
import scipy.stats as stats
import statsmodels.api as sm
import glob
import pickle
from collections import Counter
from itertools import combinations
from pprint import pprint
from Bio import SeqIO
import random as r
from itertools import permutations, product, chain
from scipy.spatial.distance import pdist, squareform
# from sklearn.decomposition import PCA
# from Bio.KEGG import REST
# from Bio.KEGG.KGML import KGML_parser
# os.chdir('...<project_dir_path>...')

'''
The goal of this script is to summarize the gene content of the samples
starting from collecting detected genes and their annotations and then
making a matrix of samples as rows and genes as columns with the values
and other data about the genes that can then be used by the script 11_gene_content.R
and others for further analysis and visualization
There is also a section that adds extra annotation to the certain output files coming
from the scripts 11_gene_content.R (masslin2 and the kruskal test randomforest analyses)
'''
# 1. Collect the detected genes and their annotations

# 1a. Collect the detected genes from the bowtie2 output and create a pickle file containing the detected genes

samples = [x.split('/')[-1].split('_gene_coverage')[0] for x in glob.glob('results/08_gene_content/01_profiling_bowtie2/*_gene_coverage.txt')]
'''
os.makedirs('results/figures/08-summarize_functions/genes_detected', exist_ok=True)
genes_detected = {sample: set() for sample in samples}
for i, sample in enumerate(samples):
        # get a list of genes that are detected in each sample
        # this comes from either a countin tool
        # or selecting the genes that have a good enough
        # mapping score, coverage and breadth from the read against whole
        # assembly read mapping. Then see how much of the "detected" genes
        # are annotated by DRAM with a KEGG ID and how many of those are
        # in the gene catalog / cd-hit clustering output anything that is
        # not in the gene catalog is not considered because it must have
        # been in a contig that was filtered out (by whokaryote or kaiju)
        # so the list of detected genes should be restricted to the gene
        # catalog explore the genes that were detected but not in the gene
        # catalog later ..
        # to do this we can count just the genes specified in the output of
        # prodigal_filt_orfs which is located in:
        # results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.gff
        # The output in results/08_gene_content/01_profiling_bowtie2/{sample}_gene_coverage.txt'
        # contains the columns scaffold, start, end, <NODE-ID>_<gene_id>, number of alignments in that region,
        # number of non-zero covered bases, length of region and fraction covered
    # write a table with the number of genes per sample
    print(f'collecting genes in {sample}')
    genes_df = pd.read_csv(f'results/08_gene_content/01_profiling_bowtie2/{sample}_gene_coverage.txt', sep = '\t', header = None)
    genes_df.columns = ['scaffold', 'start', 'end', 'id', 'num_reads', 'num_bases', 'length', 'cov_fraction']
    genes_df['gene'] = [f'{sample}_']*len(genes_df['scaffold']) + genes_df['scaffold'].apply(lambda x: str(x)) + ['_']*len(genes_df['id']) + genes_df['id'].apply(lambda x: x.split('_')[1])
    genes_df.to_csv(f'results/figures/08-summarize_functions/genes_detected/{sample}_df_detected_genes_info.csv')
    # only consider a gene detected if there are >5 reads (1000/150) and at least 90% of the bases in the gene are covered
    genes_detected[sample] = set(genes_df.loc[(genes_df['num_reads'] > 5) & (genes_df['cov_fraction'] > 0.90)]['gene'].values)
pickle.dump(genes_detected, open('results/figures/08-summarize_functions/genes_detected.pkl', 'wb'))
'''
genes_detected = pickle.load(open('results/figures/08-summarize_functions/genes_detected.pkl', 'rb'))

# 1b. Collect kegg, pfam and cazy annotations of all genes from DRAM output
        # at results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv
samples = [x.split('/')[-2] for x in glob.glob('results/08_gene_content/02_DRAM_annotations/*/annotations.tsv')]
# get functional information about each gene
rank_dict = {}  # referred to using ko_ind
kegg_dict = {}
dram_cazy_dict = {}
pfam_dict = {}
for i, sample in enumerate(samples):
    file = f'results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv'
    print(f'working on {sample} {i}/{len(samples)}', end = '\r')
    with open(file, 'r') as f:
        header_read = False
        gene_ind = 0
        ko_ind = 0
        kegg_ind = 0
        cazy_ind = 0
        pfam_ind = 0
        for line in f:
            line = line.strip()
            if not header_read:
                header = line
                header_read = True
                try:
                    gene_ind = header.split('\t').index('fasta')
                    ko_ind = header.split('\t').index('ko_id')
                    kegg_ind = header.split('\t').index('kegg_hit')
                    cazy_ind = header.split('\t').index('cazy_hits')
                    pfam_ind = header.split('\t').index('pfam_hits')
                except:
                    print(f'could not find one of the columns in {file}')
                continue
            gene = line.split('\t')[gene_ind]
            if f'{sample}_{sample}' in gene:
                gene = f'{sample}_' + gene.split(f'{sample}_')[2]
            if ko_ind != 0:
                ko = line.split('\t')[ko_ind]
                rank_dict[gene] = ko
            if kegg_ind != 0:
                kegg = line.split('\t')[kegg_ind]
                kegg_dict[gene] = kegg
            if cazy_ind != 0:
                cazy = line.split('\t')[cazy_ind].split('.hmm')[0]
                dram_cazy_dict[gene] = cazy
            if pfam_ind != 0:
                pfam = line.split('\t')[pfam_ind]
                pfam_dict[gene] = pfam

# 1c. Collect cazy annotations from Cayman output

with open('data/cayman_gene_db/20230313_gene_catalog_db.tsv', 'r') as f:
    cayman_info_dict = {line.split()[0]: line.split()[3] for line in f}

# 1d. Collect the MAGs that the genes are in from the MAGs collection 
        # and info about their species and genus
        # at results/09_MAGs_collection/all_scaffold_to_bin.tsv

all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')
handmade_spec_names = pd.read_csv('results/figures/handmade_species_names.csv')
handmade_spec_names['magotu'] = handmade_spec_names['cluster'].astype(str)
all_scaffolds_info = all_scaffolds_info.merge(handmade_spec_names, on='magotu', how='left')

mag_genus_dict = all_scaffolds_info.set_index('mag').to_dict()['genus']
mag_species_dict = all_scaffolds_info.set_index('mag').to_dict()['MAG_species_name_final']

binned_scaffolds = set(all_scaffolds_info['scaffold'])
scaffold_mag_dict = all_scaffolds_info.set_index('scaffold').to_dict()['mag']

phylo_metadata = pd.read_csv('results/11_phylogenies/phylo_genomes_metadata.tsv', sep='\t')

# 1e. Collect the information about the genes from the OrthoFinder output
        # at results/11_phylogenies/02_orthofinder_results/{group}/Results_{group}/Orthogroups/Orthogroups.txt

gene_og_dict = {}
for group in phylo_metadata['Phylogeny_group'].unique():
    og_file = f'results/11_phylogenies/02_orthofinder_results/{group}/Results_{group}/Orthogroups/Orthogroups.txt'
    if not os.path.isfile(og_file):
        print(f'{og_file} does not exist, skipping')
        continue
    for line in open(og_file, 'r'):
        if line.startswith('OG'):
            line = line.strip()
            og = line.split(': ')[0]
            genes = line.split(' ')[1:]
            for gene in genes:
                gene_og_dict[gene] = f'{group}--{og}'

# 1f. Collect the information about the coreness of the OGs from the motupan output
        # at results/08_gene_content/07_OG_coreness/{group}_motupan_output.tsv

og_coreness = {}
for group in phylo_metadata['Phylogeny_group'].unique():
    motupan_file = f'results/08_gene_content/07_OG_coreness/{group}_motupan_output.tsv'
    header_read = False
    if os.path.exists(motupan_file):    
        with open(motupan_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line == '':
                    continue
                if not header_read:
                    header_read = True
                    header = line
                    continue
                else:
                    og = line.split('\t')[0]
                    coreness = line.split('\t')[1]
                    og_coreness[f'{group}--{og}'] = coreness

# 1g. Collect the information about the clusters of the genes from the cd-hit output
                    
# 1g (1) Define the cluster class and read the cd-hit output pickle file

class cdhit_cluster:
    def __init__(self, name):
        self.name = name #cluster name
        self.genes = [] # list of genes in the cluster
        self.lengths = {} # list of lengths of the genes in the cluster
        self.identity = {} # list of perc identity to rep gene
        self.strandedness = {} # strandedness wrt rep
        self.rep_gene = '' # representative gene for the cluster
    
    def add_info(self, gene_to_add, length, is_rep, strand='+', perc=100):
        self.genes.append(gene_to_add)
        self.lengths[gene_to_add] = length
        self.identity[gene_to_add] = perc
        self.strandedness[gene_to_add] = strand
        if is_rep:
            self.rep_gene = gene_to_add
    
    def __repr__(self):
        return f'Cluster--{self.name}:\n  rep gene:{self.rep_gene}\n    genes: {self.genes}\n    lengths: {self.lengths}\n    identity: {self.identity}'

    def __str__(self):
        return f'Cluster--{self.name}:\n  rep gene:{self.rep_gene}\n    genes: {self.genes}\n    lengths: {self.lengths}\n    identity: {self.identity}'

clusters = pd.read_pickle('results/figures/clusters_objects_collection.pkl')

gene_cd_hit_dict = {}
for cluster in clusters:
    for gene in clusters[cluster].genes:
        gene_cd_hit_dict[gene] = cluster

# 2. Parse and write to the gene_info_tables
        # at 'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv'

# 2a. Define functions to use and create / read info required by them

'''
gene_position_names_ffn = {} # (start, end)
for sample in samples:
    with open(f'results/figures/08-summarize_functions/genes_detected/{sample}_df_detected_genes_info.csv') as f:
        for line in f:
            line = line.strip()
            start = line.split(',')[2]
            end = line.split(',')[3]
            gene = line.split(',')[9]
            ffn_name = '_'.join(gene.split('_')[:-1]) + ':' + start + '-' + end
            gene_position_names_ffn[gene] = ffn_name
pickle.dump(gene_position_names_ffn, open('results/figures/08-summarize_functions/gene_position_names_ffn.pkl', 'wb'))
'''
gene_position_names_ffn = pickle.load(open('results/figures/08-summarize_functions/gene_position_names_ffn.pkl', 'rb'))

def scaffold_from_geneid(geneid):
    return '_'.join(geneid.split('_')[:-1])

# get the mag info for the scaffold
def get_mag_of_gene(gene):
    if scaffold_from_geneid(gene) in binned_scaffolds:
        return scaffold_mag_dict[scaffold_from_geneid(gene)]
    else:
        return 'unbinned'

'''
gene_renamed_dict = {}
# the faa and ffn files from results/11_phylogenies/00_genomes/ have different IDs and a different 
# number of genes than the files that were used to name and collected "detected genes"
# so use the gene positions of detected genes from 
# 'results/figures/08-summarize_functions/genes_detected/{sample}_df_detected_genes_info.csv'
# to then get the re-name of the gene in the ffn file based on its name in the ffn_original

# use the ffn to get the old and new names!
for i, mag in enumerate(glob.glob('results/11_phylogenies/00_genomes/*')):
    mag_name = mag.split('/')[-1]
    print(f'working on mag {i}/{len(glob.glob("results/11_phylogenies/00_genomes/*"))}', end = '\r')
    gene_renamed_dict[mag_name] = {}
    with open(f'{mag}/{mag_name}.ffn', 'r') as f, open(f'{mag}/{mag_name}_original.ffn', 'r') as f_o:
        for line, line_o in zip(f, f_o):
            if line.startswith('>') and line_o.startswith('>'):
                gene = line.strip().split('>')[1]
                gene_o = line_o.strip().split('>')[1]
                gene_renamed_dict[mag_name][gene_o] = gene
pickle.dump(gene_renamed_dict, open('results/figures/08-summarize_functions/gene_renamed_dict.pkl', 'wb'))
'''
gene_renamed_dict = pickle.load(open('results/figures/08-summarize_functions/gene_renamed_dict.pkl', 'rb'))

def og_gene_id(gene):
    mag = get_mag_of_gene(gene)
    if mag == 'unbinned':
        return(None)
    else:
        if mag in gene_renamed_dict.keys(): 
            if gene_position_names_ffn[gene] in gene_renamed_dict[mag]:
                gene_id = gene_renamed_dict[mag][gene_position_names_ffn[gene]]
                return(gene_id)
    return(None)

# gene='C2-3_NODE_740_length_70557_cov_30.867252_28'
# 'C2-3_18_28' - wrong
# 'C2-3_18_1179' - right
# 'C4-3_NODE_4_length_644241_cov_51.191471_609' - was not found in gene_renamed_dict
# gene_positions_dict['C4-3_6']
# og_gene_id('C2-3_NODE_740_length_70557_cov_30.867252_28')
# og_gene_id('C4-3_NODE_4_length_644241_cov_51.191471_609')

# exception : 'C4-3_NODE_174_length_110371_cov_13.960894:94400-95152' gave key error before
# if gene_position_names_ffn[gene] in gene_renamed_dict[mag]: was added - dig in

os.makedirs('results/figures/08-summarize_functions/gene_info_tables', exist_ok=True)
for i, sample in enumerate(samples):
    if os.path.isfile(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv'):
        continue
    print(f'workin on {sample} {i}/{len(samples)}', end = '\r')
    df_detected_genes_info = pd.DataFrame()
    df_detected_genes_info['gene'] = [x for x in genes_detected[sample]]
    df_detected_genes_info['ko'] = df_detected_genes_info['gene'].map(kegg_dict)
    df_detected_genes_info['cazyme-DRAM'] = df_detected_genes_info['gene'].map(dram_cazy_dict)
    df_detected_genes_info['cazyme'] = df_detected_genes_info['gene'].map(cayman_info_dict)
    df_detected_genes_info['mag'] = pd.Series([get_mag_of_gene(gene) for gene in genes_detected[sample]])
    df_detected_genes_info['species'] = df_detected_genes_info['mag'].map(mag_species_dict)
    df_detected_genes_info['genus'] = df_detected_genes_info['mag'].map(mag_genus_dict)
    og_info_dict = {}
    for gene in genes_detected[sample]:
        if og_gene_id(gene) is not None and og_gene_id(gene) in gene_og_dict.keys():
            og_info_dict[gene] = gene_og_dict[og_gene_id(gene)]
    df_detected_genes_info['og'] = df_detected_genes_info['gene'].map(og_info_dict)
    df_detected_genes_info['coreness'] = df_detected_genes_info['og'].map(og_coreness)
    df_detected_genes_info['cluster'] = df_detected_genes_info['gene'].map(gene_cd_hit_dict)
    df_detected_genes_info.to_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')


# # delete large dicts from memory
# del mag_species_dict
# del mag_genus_dict
# del binned_scaffolds
# del scaffold_mag_dict
# del gene_renamed_dict
# del gene_og_dict
# del og_coreness