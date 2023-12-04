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
import random
# from Bio.KEGG import REST
# from Bio.KEGG.KGML import KGML_parser
# os.chdir('/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison')
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

samples = [x.split('/')[-1].split('_gene_coverage')[0] for x in glob.glob('results/08_gene_content/01_profiling_bowtie2/*_gene_coverage.txt')]
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
    # only consider a gene detected if there are >5 reads (1000/150) and at least 90% of the bases in the gene are covered
    genes_detected[sample] = set(genes_df.loc[(genes_df['num_reads'] > 5) & (genes_df['cov_fraction'] > 0.90)]['gene'].values)

genes_per_sample = {}
for sample in samples:
    genes_per_sample[sample] = len(genes_detected[sample])
with open('results/figures/08-summarize_functions/genes_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tgenes\n')
    for sample in genes_per_sample:
        out_fh.write(f'{sample}\t{genes_per_sample[sample]}\n')

# now number of clusters per sample
'''
define a class to read and handle clusters. Each cluster represents a list of genes and is represented
'''

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
# for each cluster read the names of genes in cluster 
# if any of them are among the genes detected in the sample,

# add the cluster to the set of clusters detected in that sample
# clusters_per_sample = {sample: set() for sample in samples}
# for i, sample in enumerate(samples):
#     print(f'working on {sample} {i}/{len(samples)}', end = '\r')
#     for cluster in clusters:
#         for gene in clusters[cluster].genes:
#             if gene in genes_detected[sample]:
#                 clusters_per_sample[sample].add(cluster)
#                 break
# pickle.dump(clusters_per_sample, open('results/figures/08-summarize_functions/clusters_per_sample.pkl', 'wb'))

clusters_per_sample = pd.read_pickle('results/figures/08-summarize_functions/clusters_per_sample.pkl')

with open('results/figures/08-summarize_functions/clusters_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tclusters\n')
    for sample in clusters_per_sample:
        out_fh.write(f'{sample}\t{len(clusters_per_sample[sample])}\n')


with open('results/figures/08-summarize_functions/cluster_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tgenes\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(10):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/11 and sample size {j}/{len(samples_host)}', end = '\r')
                cum_clusters = 0
                samples_iter = random.sample(samples_host, sample_size)
                clusters_detected_iter = set()
                for sample in samples_iter:
                    clusters_detected_iter.update(clusters_per_sample[sample])
                cum_clusters = len(clusters_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_clusters}\n')

# get functional information about each gene
rank_dict = {}
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

# prepare the disctionary to get kegg ortholog information
current_a = ''
current_b = ''
current_c = ''
kegg_info_dict = {}
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
            kegg_info_dict[ko] = {'A': current_a, 'B': current_b, 'C': current_c}

with open('results/figures/08-summarize_functions/kegg_info_dict.txt', 'w+') as f:
    header = 'ko\tA\tB\tC\n'
    f.write(header)
    for ko in kegg_info_dict.keys():
        line = f'{ko}\t{kegg_info_dict[ko]["A"]}\t{kegg_info_dict[ko]["B"]}\t{kegg_info_dict[ko]["C"]}\n'
        f.write(line)

# count the number of kos and under each category for each sample
kos_detected = {sample: {} for sample in samples}
ko_count = {sample: {'Total': 0, 'Annotated': 0} for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}')
    for gene in genes_detected[sample]:
        # kos_detected[sample]['NA'] = set()
        kos_detected[sample]['NA'] = 0
        ko_gene = kegg_dict[gene]
        ko_count[sample]['Total'] += 1
        if ko_gene == '':
            # kos_detected[sample]['NA'].add(gene)
            kos_detected[sample]['NA'] += 1
        else:
            ko_count[sample]['Annotated'] += 1
            if ko_gene in kos_detected[sample].keys():
                kos_detected[sample][ko_gene] += 1
            else:
                kos_detected[sample][ko_gene] = 1
                # kos_detected[sample][ko_gene] = set([gene])

# write an output to summarize Total and annotated genes from ko_count

with open('results/figures/08-summarize_functions/kos_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tkos\tunannotated\n')
    for sample in clusters_per_sample:
        ann = str(ko_count[sample]['Annotated'])
        unann = str(ko_count[sample]['Total'] - ko_count[sample]['Annotated'])
        out_fh.write(f'{sample}\t{ann}\t{unann}\n')


# summarize number of kos found per sample
with open('results/figures/08-summarize_functions/kos_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tkos\n')
    for sample in kos_detected:
        out_fh.write(f'{sample}\t{len(kos_detected[sample].keys())}\n')

# cum curve of #ko and #category
num_iters = 20
with open('results/figures/08-summarize_functions/KO_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tkos\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(num_iters):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/{num_iters} and sample size {j}/{len(samples_host)}', end = '\r')
                cum_kos = 0
                samples_iter = random.sample(samples_host, sample_size)
                kos_detected_iter = set()
                for sample in samples_iter:
                    kos_detected_iter.update([x for x in kos_detected[sample].keys()])
                cum_kos = len(kos_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_kos}\n')

# beta div compare ko - sample vs ko matrix

# pick up enriched kos and / categories

##############################################################################################

# count the number of kos and under each category for each sample
kos_A_detected = {sample: {} for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}')
    for gene in genes_detected[sample]:
        ko_cat = 'NA'
        # kos_detected[sample]['NA'] = set()
        kos_A_detected[sample]['NA'] = 0
        ko_gene = kegg_dict[gene]
        if ko_gene == '':
            kos_A_detected[sample]['NA'] += 1
        else:
            ko_cat = kegg_info_dict[ko_gene]['A']
            if 'Unclassified: ' in ko_cat:
                ko_cat = ko_cat.split('Unclassified: ')[1]
        if ko_cat in kos_A_detected[sample].keys():
            kos_A_detected[sample][ko_cat] += 1
        else:
            kos_A_detected[sample][ko_cat] = 1
            # kos_A_detected[sample][ko_gene] = set([gene])

# summarize number of kos found per sample
with open('results/figures/08-summarize_functions/kos_A_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tkos\tcount\n')
    for sample in kos_A_detected:
        for cat in kos_A_detected[sample]:
            out_fh.write(f'{sample}\t{cat}\t{kos_A_detected[sample][cat]}\n')

##############################################################################################

# count the number of kos and under each category for each sample
kos_B_detected = {sample: {} for sample in samples}
ko_B_count = {sample: {'Total': 0, 'Annotated': 0} for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}')
    for gene in genes_detected[sample]:
        ko_cat = 'NA'
        # kos_detected[sample]['NA'] = set()
        kos_B_detected[sample]['NA'] = 0
        ko_gene = kegg_dict[gene]
        if ko_gene == '':
            kos_B_detected[sample]['NA'] += 1
        else:
            ko_cat = kegg_info_dict[ko_gene]['B']
            if 'Unclassified: ' in ko_cat:
                ko_cat = ko_cat.split('Unclassified: ')[1]
            ko_B_count[sample]['Total'] += 1
            ko_B_count[sample]['Annotated'] += 1
        if ko_cat in kos_B_detected[sample].keys():
            kos_B_detected[sample][ko_cat] += 1
        else:
            kos_B_detected[sample][ko_cat] = 1
            # kos_B_detected[sample][ko_gene] = set([gene])

# summarize number of kos found per sample
with open('results/figures/08-summarize_functions/kos_B_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tkos\tcount\n')
    for sample in kos_B_detected:
        for cat in kos_B_detected[sample]:
            out_fh.write(f'{sample}\t{cat}\t{kos_B_detected[sample][cat]}\n')

# cum curve of #ko and #category
num_iters = 20
with open('results/figures/08-summarize_functions/KO_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tkos\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(num_iters):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/{num_iters} and sample size {j}/{len(samples_host)}', end = '\r')
                cum_kos = 0
                samples_iter = random.sample(samples_host, sample_size)
                kos_B_detected_iter = set()
                for sample in samples_iter:
                    kos_B_detected_iter.update([x for x in kos_B_detected[sample].keys()])
                cum_kos = len(kos_B_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_kos}\n')


##############################################################################################

# count the number of kos and under each category for each sample
kos_C_detected = {sample: {} for sample in samples}
ko_C_count = {sample: {'Total': 0, 'Annotated': 0} for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}')
    for gene in genes_detected[sample]:
        ko_cat = 'NA'
        # kos_detected[sample]['NA'] = set()
        kos_C_detected[sample]['NA'] = 0
        ko_gene = kegg_dict[gene]
        if ko_gene == '':
            kos_C_detected[sample]['NA'] += 1
        else:
            ko_cat = kegg_info_dict[ko_gene]['C']
            if 'Unclassified: ' in ko_cat:
                ko_cat = ko_cat.split('Unclassified: ')[1]
            ko_C_count[sample]['Total'] += 1
            ko_C_count[sample]['Annotated'] += 1
        if ko_cat in kos_C_detected[sample].keys():
            kos_C_detected[sample][ko_cat] += 1
        else:
            kos_C_detected[sample][ko_cat] = 1
            # kos_C_detected[sample][ko_gene] = set([gene])

# summarize number of kos found per sample
with open('results/figures/08-summarize_functions/kos_C_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tkos\tcount\n')
    for sample in kos_C_detected:
        for cat in kos_C_detected[sample]:
            out_fh.write(f'{sample}\t{cat}\t{kos_C_detected[sample][cat]}\n')

# cum curve of #ko and #category
num_iters = 50
with open('results/figures/08-summarize_functions/KO_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tkos\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(num_iters):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/{num_iters} and sample size {j}/{len(samples_host)}', end = '\r')
                cum_kos = 0
                samples_iter = random.sample(samples_host, sample_size)
                kos_C_detected_iter = set()
                for sample in samples_iter:
                    kos_C_detected_iter.update([x for x in kos_C_detected[sample].keys()])
                cum_kos = len(kos_C_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_kos}\n')

# beta div compare ko - sample vs ko matrix

# pick up enriched kos and / categories

# which taxa contributes?

# is it core?

# intermediate files in 'results/figures/08-summarize_functions'

#### KEGG category level A
# kegg_dict_A = {}
# for gene in kegg_dict.keys():
#     if kegg_dict[gene] == '':
#         kegg_dict_A[gene] = ''
#     else:
#         kegg_dict_A[gene] = kegg_info_dict[kegg_dict[gene]]['A']

################################################################################

'''
From this point on the goal is to make a table (per sample)
to summarize all the genes detected, their differen annotations
according to:
    KO (ko) - sub categories, COGs (?)
    ...
    CAzy (cazyme-dram) - from DRAM
    CAzy (cazyme) - from Cayman
    MAG (mag) - if it is found in one
    Species
    Genus
    Coreness
And then we collect the information from these table and make
summaries for further plots and then the tables themselves can
also be analyzed
'''

def scaffold_from_geneid(geneid):
    return '_'.join(geneid.split('_')[:-1])

# get the mag info for the scaffold
def get_mag_of_gene(gene):
    if scaffold_from_geneid(gene) in binned_scaffolds:
        return scaffold_mag_dict[scaffold_from_geneid(gene)]
    else:
        return 'unbinned'

def get_gene_id_num(gene):
    return gene.split('_')[-1]

def og_gene_id(gene):
    mag = get_mag_of_gene(gene)
    if mag == 'unbinned':
        return(None)
    else:
        return f'{mag}_{get_gene_id_num(gene)}'

all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')
handmade_spec_names = pd.read_csv('results/figures/handmade_species_names.csv')
handmade_spec_names['magotu'] = handmade_spec_names['cluster'].astype(str)
all_scaffolds_info = all_scaffolds_info.merge(handmade_spec_names, on='magotu', how='left')all_scaffolds_info.columns

mag_genus_dict = all_scaffolds_info.set_index('mag').to_dict()['genus']
mag_species_dict = all_scaffolds_info.set_index('mag').to_dict()['MAG_species_name_final']

binned_scaffolds = set(all_scaffolds_info['scaffold'])
scaffold_mag_dict = all_scaffolds_info.set_index('scaffold').to_dict()['mag']

# We need to find out if the gene is in a MAG and whether the mag is in the OG analysis and then what the OG is
# get the MAGs that are in the OG analysis
phylo_metadata = pd.read_csv('results/11_phylogenies/phylo_genomes_metadata.tsv', sep='\t')

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

with open('data/cayman_gene_db/20230313_gene_catalog_db.tsv', 'r') as f:
    cayman_info_dict = {line.split()[0]: line.split()[3] for line in f}

gene_cd_hit_dict = {}
for cluster in clusters:
    for gene in clusters[cluster].genes:
        gene_cd_hit_dict[gene] = cluster


# Takes a while but after it is written
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
    df_detected_genes_info['coreness'] = df_detected_genes_info['gene'].map(og_coreness)
    df_detected_genes_info['cluster'] = df_detected_genes_info['gene'].map(gene_cd_hit_dict)
    df_detected_genes_info.to_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')

del mag_species_dict
del mag_genus_dict
del binned_scaffolds
del scaffold_mag_dict

def sample_feature_matrix(feature_name='ko', path=''):
    '''
    feature can be one of gene, ko, cazyme-DRAM, cazyme,
    og, coreness, cluster i.e. one of the columns of the
    df_detected_genes_info table
    read all the sample tables and make a matrix of samples
    as rows and features as columns with the values being
    the number of genes detected in that sample for that
    feature
    '''
    samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]
    '''
        first go through all the sample and collect all the kos
        read the df with list of genes and info for each sample and count the number of genes
        for each feature and add it to the feature matrix if the feature is not present in the
        feature matrix, add it as a column for the row of this sample and set 0 for all other
        samples if the feature is present in the feature matrix, add the count for this sample
        to the column for this feature
    '''
    ko_columns = set()
    for i, sample in enumerate(samples):
        print(f'collecting {feature_name} from {sample} {i}/{len(samples)}', end = '\r')
        df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        ko_columns.update(df[feature_name].unique())
    print(f'finished collecting {feature_name} from all samples')
    feature_matrix = pd.DataFrame()
    feature_matrix['sample'] = samples
    feature_matrix = feature_matrix.set_index('sample')
    feature_names = list(ko_columns)
    feature_matrix = feature_matrix.reindex(columns = feature_names)
    feature_matrix = feature_matrix.fillna(0)
    for i, sample in enumerate(samples):
        print(f'working on {feature_name} for {sample} {i}/{len(samples)}', end = '\r')
        df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        # remove the column unnamed and set gene as index
        df = df.drop(columns = ['Unnamed: 0'])
        df = df.set_index('gene')
        # get the counts for each feature
        feature_counts = df[feature_name].value_counts()
        # add the counts to the feature matrix
        for feature in feature_counts.keys():
            feature_matrix.loc[sample, feature] = feature_counts[feature]
    print(f'finished counting {feature_name} from all samples')
    if path != '':
        # os.makedirs(os.path.dirname(path), exist_ok=True)
        feature_matrix.to_csv(path)
    return feature_matrix

sample_feature_matrix(feature_name='ko', path='results/figures/08-summarize_functions/ko_matrix.csv')
sample_feature_matrix(feature_name='cluster', path='results/figures/08-summarize_functions/cluster_matrix.csv')
sample_feature_matrix(feature_name='cazyme-DRAM', path='results/figures/08-summarize_functions/cazyme-DRAM_matrix.csv')
sample_feature_matrix(feature_name='cazyme', path='results/figures/08-summarize_functions/cazyme_matrix.csv')
sample_feature_matrix(feature_name='og', path='results/figures/08-summarize_functions/og_matrix.csv')

        

'''
make summary table of gene content information we have the following information about genes:
the gene catalog:
    results/08_gene_content/20230313_gene_catalog.ffn
profile of the gene catalog (per sample): 
    results/08_gene_content/06_cayman/{sample}.gene_counts.txt.gz:
the gene catalog with cayman annotation:
    data/cayman_gene_db/20230313_gene_catalog_db.csv (No header, columns: gene_id, start, end, gene_length, pval, family, length)
dram annotation of the gene catalog
    results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv - read using pandas!
    # check how many genes were not annotated for KEGG
name of ko annotations and their associated pathways are in data/
    data/KEGG_info/ko00001.json
we need to obtain the following information about genes:

# OG membership
gene -> contig -> mag -> id_in_og_analysis -> OG

# taxonomic information
gene -> contig -> mag -> gtdb tax from mag summary at results/09_MAGs_collection/MAGs_metadata_summary.tsv

# contig info
gene -> contig -> results/05_assembly/contig_fates/ -> kaiju/kraken/whokaryote (parse out which sample to look instide)


gene_id refers to IDs of the form: {sample}_NODE_10481_length_4970_cov_52.067955_3
    i.e. {sample}_{contig_name}_3 or {contig_id}_3
gene_id_num refers to the number at the end of the gene_id
    i.e. 3
og_gene_id refers to the gene ID in the orthogroup analysis: C3-5_4_1192
    i.e. {mag}_1192
OG refers to the orthogroup ID
    such as OG0000035
(not thinking about CD-hit for now)    
'''


'''
all_scaffold_to_bin.tsv contains the following columns
scaffold	mag	size	completeness	contamination	genus	species	magotu
compiled from across all scaffolds and includes information about the mag it belongs to
and the respective information of the mag collected using a temporary script not in this
to be recovered or rewritten later... (magotu is in the numerical form coming straight from dRep)
''' 
all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')

# cazy_dict = pd.read_csv('data/cayman_gene_db/20230313_gene_catalog_filtered_filtered_merged_CORRECT_fold_stuff.csv').set_index('sequenceID').to_dict()['family']


len(clusters)
# 710,497 in total

# show a bar plot with x axis the number of hosts that the cluster is found in and y axis the number of clusters
# that are found in that many hosts

# get the number of hosts that each cluster is found in
# singletons = 0
# cluster_size_dict = {}
# for cluster in clusters:
#     cluster_size_dict[cluster] = len(clusters[cluster].genes)
#     if len(clusters[cluster].genes) == 1:
#         singletons += 1

# plt.hist(cluster_size_dict.values(), bins=100)
# plt.savefig('results/figures/08-summarize_functions/cluster_size_hist.png')
# plt.close()

# cluster_host_dict = {}
# for cluster in clusters:
#     if len(clusters[cluster].genes) == 1:
#         continue
#     hosts_seen = Counter([x[0] for x in clusters[cluster].genes])
#     cluster_host_dict[cluster] = len(hosts_seen.keys())
# plt.hist(cluster_host_dict.values(), bins = 5)
# plt.savefig('results/figures/08-summarize_functions/cluster_host_hist.png')
# plt.close()

# do the above but using mapping result to tell if a cluster has been detected
'''
A cluster is considered detected if any gene in the cluster is found at a coverage > 0
and breadth > 0.75 in the sample
To find this, for each sample, 
'''



# # get a histogram of coverages
# plt.hist(gene_info_dfs['A1-1'].coverage, bins=10000)
# plt.xlim(0, 1)
# plt.savefig('results/figures/08-summarize_functions/coverage_hist.png')
# plt.close()


# '''
# copying some files for future use
# '''

# with open('results/export/magOTU_handmade_names.csv', 'r') as f:
#     header = f.readline()
#     for line in f:
#         line = line.strip()
#         mag = line.split(',')[0]
#         mag_name = line.split(',')[1]
#         faa_path = line.split(',')[5]
#         gff_path = line.split(',')[6]
#         print(mag)
#         out_faa=f'results/export/MAGs_annotation_faa/{mag_name}.faa'
#         out_gff=f'results/export/MAGs_annotation_gff/{mag_name}.gff'
#         out_faa_nas=f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/MAGs_collection/MAGs_faa/{mag_name}.faa'
#         out_gff_nas=f'/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/MAGs_collection/MAGs_gff/{mag_name}.gff'
#         os.system(f'cp {faa_path} {out_faa}')
#         os.system(f'cp {gff_path} {out_gff}')
#         os.system(f'cp {faa_path} {out_faa_nas}')
#         os.system(f'cp {gff_path} {out_gff_nas}')
