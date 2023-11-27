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
# from Bio.KEGG import REST
# from Bio.KEGG.KGML import KGML_parser


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

def scaffold_from_geneid(geneid):
    return '_'.join(geneid.split('_')[:-1])

# get the mag info for the scaffold
def get_mag_of_gene(gene):
    try:
        mag = all_scaffolds_info[all_scaffolds_info['scaffold'] == scaffold_from_geneid(gene)]['mag'].values[0]
    except:
        mag = 'unbinned'
    return(mag)

def get_gene_id_num(gene):
    return gene.split('_')[-1]

def og_gene_id(gene):
    mag = get_mag_of_gene(gene)
    if mag == 'unbinned':
        return(None)
    else:
        return f'{mag}_{get_gene_id_num(gene)}'

# for this we need to find out if the gene is in a MAG and whether the mag is in the OG analysis and then what the OG is
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

# combine gene_og_dict and og_coreness and write as a table
with open('results/08_gene_content/07_OG_coreness/gene_og_coreness.tsv', 'w') as out_fh:
    out_fh.write('gene\tog\tcoreness\n')
    for gene in gene_og_dict:
        if gene_og_dict[gene] not in og_coreness:
            coreness = "NA"
        else:
            coreness = og_coreness[gene_og_dict[gene]]
        out_fh.write(f'{gene}\t{gene_og_dict[gene]}\t{coreness}\n')

# samples = [x.split('/')[-2] for x in glob.glob('results/10_instrain/02_instrain_profile/*/')]
# gene_info_dfs = {}
# for sample in samples:
#     gene_info_dfs[sample] = pd.read_csv(f"results/10_instrain/02_instrain_profile/{sample}/output/{sample}_gene_info.tsv", sep = "\t")


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

# cluster_file = "results/08_gene_content/00_cdhit_clustering/20230313_gene_catalog_cdhit9590.fasta.clstr"
# clusters = {}
# cluster_num = 0
# with open(cluster_file, 'r') as f:
#     for line in f:
#         line = line.strip()
#         if line.startswith('>'):
#             cluster_num = line.split('>Cluster ')[1]
#             clust_obj = cdhit_cluster(cluster_num)
#             clusters[cluster_num] = clust_obj
#         else:
#             length = int(line.split()[1].split('nt')[0])
#             gene = line.split()[2].split('>')[1].split('...')[0]
#             rep_char = line.split()[3]
#             if rep_char == '*':
#                 clusters[cluster_num].add_info(gene, length, True)
#             else:
#                 strand = line.split()[4].split('/')[1]
#                 perc = float(line.split()[4].split('/')[2].split('%')[0])
#                 clusters[cluster_num].add_info(gene, length, False, strand, perc)
# with open('results/figures/clusters_objects_collection.pkl', "wb") as pkl_fh:
#         pickle.dump(clusters, pkl_fh)
clusters = pd.read_pickle('results/figures/clusters_objects_collection.pkl')
gene_cluster_dict = {}
for cluster in clusters:
    for gene in clusters[cluster].genes:
        gene_cluster_dict[gene] = cluster

'''
all_scaffold_to_bin.tsv contains the following columns
scaffold	mag	size	completeness	contamination	genus	species	magotu
compiled from across all scaffolds and includes information about the mag it belongs to
and the respective information of the mag collected using a temporary script not in this
to be recovered or rewritten later... (magotu is in the numerical form coming straight from dRep)
''' 
all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')

cazy_dict = pd.read_csv('data/cayman_gene_db/20230313_gene_catalog_filtered_filtered_merged_CORRECT_fold_stuff.csv').set_index('sequenceID').to_dict()['family']
# get functional information about each gene
anno_files = glob.glob(f'results/08_gene_content/02_DRAM_annotations/*/annotations.tsv')
rank_dict = {}
kegg_dict = {}
dram_cazy_dict = {}
pfam_dict = {}
for file in anno_files:
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
                    cazy_ind = header.split('\t').index('cazy_best_hit')
                    pfam_ind = header.split('\t').index('pfam_hits')
                except:
                    print(f'could not find one of the columns in {file}')
                continue
            gene = line.split('\t')[gene_ind]
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

gene_catalog = 'results/08_gene_content/20230313_gene_catalog.ffn'
total_genes = 0
genes_per_sample = Counter()
gene_per_host = Counter()
all_genes = SeqIO.parse(gene_catalog, 'fasta')
for gene in all_genes:
    total_genes += 1
    genes_per_sample[gene.id.split('_NODE')[0]] += 1
    gene_per_host[gene.id[0]] += 1
# total_genes
# 4,911,821
len([kegg_dict[x] for x in kegg_dict if kegg_dict[x] != ''])
# 2,772,118
len([dram_cazy_dict[x] for x in dram_cazy_dict if dram_cazy_dict[x] != ''])
# 122,164
len(cazy_dict) # from cayman
# 150,832
len([pfam_dict[x] for x in pfam_dict if pfam_dict[x] != ''])
# 259,590
binned_scaffolds = set(all_scaffolds_info['scaffold'])
in_mag = 0
all_genes = SeqIO.parse(gene_catalog, 'fasta')
for gene in all_genes:
    if scaffold_from_geneid(gene.id) in binned_scaffolds:
        in_mag += 1
# 3,678,762
in_og = len(gene_og_dict)
# 2,511,732

'''
At this point, multiple approaches are possible.
    Using presence-absence based on a certain breadth and coverage threshold and the counting
        the number of genes for that feature in each sample
    Making a matrix of counts (summed gene coverage) of feature across each samples
    Number of clusters for each feature per sample - this is KE's approach

Make tables, matrices etc. from that can be read and plotted in R
'''

handmade_spec_names = pd.read_csv('results/figures/handmade_species_names.csv')
handmade_spec_names['magotu'] = handmade_spec_names['cluster'].astype(str)
all_scaffolds_info = all_scaffolds_info.merge(handmade_spec_names, on='magotu', how='left')

'''
# summarize the gene clusters
Since they will be used as a unit to consider the representation of functional categories
instead of number of genes....
'''

len(clusters)
# 710497 in total

# show a bar plot with x axis the number of hosts that the cluster is found in and y axis the number of clusters
# that are found in that many hosts

# get the number of hosts that each cluster is found in
singletons = 0
cluster_size_dict = {}
for cluster in clusters:
    cluster_size_dict[cluster] = len(clusters[cluster].genes)
    if len(clusters[cluster].genes) == 1:
        singletons += 1

plt.hist(cluster_size_dict.values(), bins=100)
plt.savefig('results/figures/visualize_temp/cluster_size_hist.png')
plt.close()

cluster_host_dict = {}
for cluster in clusters:
    if len(clusters[cluster].genes) == 1:
        continue
    hosts_seen = Counter([x[0] for x in clusters[cluster].genes])
    cluster_host_dict[cluster] = len(hosts_seen.keys())
plt.hist(cluster_host_dict.values(), bins = 5)
plt.savefig('results/figures/visualize_temp/cluster_host_hist.png')
plt.close()

# do the above but using mapping result to tell if a cluster has been detected
'''
A cluster is considered detected if any gene in the cluster is found at a coverage > 0
and breadth > 0.75 in the sample
To find this, for each sample, 
'''

# RESUME HERE

samples = [x.split('/')[-1].split('_mapped')[0] for x in glob.glob('results/08_gene_content/01_profiling/*_mapped.coverage')]
for sample in samples:
    input_f = f'results/08_gene_content/01_profiling/{sample}_mapped.coverage'
    output_f = f'results/08_gene_content/01_profiling/{sample}_mapped.coverage.filtered'
    # use awk to filter the coverage file
    os.system(f"awk '$6 > 50' {input_f} > {output_f}")

coverage_outputs_filt = {}
for sample in samples:
    try:
        coverage_output_filt = pd.read_csv(f'results/08_gene_content/01_profiling/{sample}_mapped.coverage.filtered', sep='\t')
        coverage_outputs_filt[sample] = coverage_output_filt
    except:
        pass

# write a table with the number of genes per sample
genes_per_sample = {}
for sample in samples:
    genes_per_sample[sample] = len(coverage_outputs_filt[sample])
with open('results/figures/visualize_temp/genes_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tgenes\n')
    for sample in genes_per_sample:
        out_fh.write(f'{sample}\t{genes_per_sample[sample]}\n')

# now number of clusters per sample
clusters_per_sample = {}
for sample in samples:

    clusters_per_sample[sample] = len(set([gene_cluster_dict[x] for x in coverage_outputs_filt[sample]['gene'].values]))

with open('results/figures/visualize_temp/clusters_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tclusters\n')
    for sample in clusters_per_sample:
        out_fh.write(f'{sample}\t{clusters_per_sample[sample]}\n')

# intermediate files in 'results/figures/functional_analysis_data'

#### KEGG category level A
kegg_dict_A = {}
for gene in kegg_dict.keys():
    if kegg_dict[gene] == '':
        kegg_dict_A[gene] = ''
    else:
        kegg_dict_A[gene] = kegg_info_dict[kegg_dict[gene]]['A']

### summarize across samples
# 2. #clusters
#       count the number of clusters per category for each sample to
#       make a table of the form:
#           sample  category    #clusters




# 1. #genes
#       count the number of genes per category for each sample to
#       make a table of the form:
#           sample  category    #genes




pprint(clusters)
# find a good way to summarize this information as a table
pprint(gene_og_dict)
pprint(og_coreness)
# combine gene_og_dict and og_coreness and write as a table
with open('results/figures/gene_og_coreness.tsv', 'w') as out_fh:
    out_fh.write('gene\tog\tcoreness\n')
    for gene in gene_og_dict:
        if gene_og_dict[gene] not in og_coreness:
            coreness = "NA"
        else:
            coreness = og_coreness[gene_og_dict[gene]]
        out_fh.write(f'{gene}\t{gene_og_dict[gene]}\t{coreness}\n')

pprint(gene_cluster_dict)
pprint(rank_dict)
pprint(kegg_dict)
pprint(dram_cazy_dict)
pprint(pfam_dict)
pprint(kegg_info_dict)

og_mags = set(phylo_metadata['ID'])


# # get a histogram of coverages
# plt.hist(gene_info_dfs['A1-1'].coverage, bins=10000)
# plt.xlim(0, 1)
# plt.savefig('results/figures/visualize_temp/coverage_hist.png')
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


'''
more notes:
    there are 3 ways to consider coverage of a gene
        one, is directly from the back-mapping to the contig but this is 
        tedious and potentially suffers from read stealing by host contigs
        second, mapping to the MAGs this give a breadth and coverage 
        for each gene that is in the rep. MAGs but means that there is information
        missing for others so this is not the most useful here
        third, is from mapping to the set of all filtered genes (gene catalog) this is
        the approach we will use right now and consider a gene of the gene catalog
        detected if it has a coverage of > 0 and breadth > 0.75 
        (reconsider the threshold later)
'''