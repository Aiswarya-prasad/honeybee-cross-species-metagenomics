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
with open('results/figures/visualize_temp/genes_per_sample.tsv', 'w') as out_fh:
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
clusters_per_sample = {sample: set() for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}', end = '\r')
    for cluster in clusters:
        for gene in clusters[cluster].genes:
            if gene in genes_detected[sample]:
                clusters_per_sample[sample].add(cluster)
                break
pickle.dump(clusters_per_sample, open('results/figures/visualize_temp/clusters_per_sample.pkl', 'wb'))
clusters_per_sample = pd.read_pickle('results/figures/visualize_temp/clusters_per_sample.pkl')
with open('results/figures/visualize_temp/clusters_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tclusters\n')
    for sample in clusters_per_sample:
        out_fh.write(f'{sample}\t{len(clusters_per_sample[sample])}\n')


with open('results/figures/visualize_temp/cluster_cum_curve.tsv', 'w+') as out_fh:
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
for sample in samples:
    file = f'results/08_gene_content/02_DRAM_annotations/{sample}/annotations.tsv'
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

# count the number of kos and under each category

# cum curve of #ko and #category

# beta div compare ko

# pick up enriched kos and / categories

# which taxa contributes?

# is it core?

# intermediate files in 'results/figures/visualize_temp'

#### KEGG category level A
# kegg_dict_A = {}
# for gene in kegg_dict.keys():
#     if kegg_dict[gene] == '':
#         kegg_dict_A[gene] = ''
#     else:
#         kegg_dict_A[gene] = kegg_info_dict[kegg_dict[gene]]['A']



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
all_scaffold_to_bin.tsv contains the following columns
scaffold	mag	size	completeness	contamination	genus	species	magotu
compiled from across all scaffolds and includes information about the mag it belongs to
and the respective information of the mag collected using a temporary script not in this
to be recovered or rewritten later... (magotu is in the numerical form coming straight from dRep)
''' 
all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')

cazy_dict = pd.read_csv('data/cayman_gene_db/20230313_gene_catalog_filtered_filtered_merged_CORRECT_fold_stuff.csv').set_index('sequenceID').to_dict()['family']


# gene_catalog = 'results/08_gene_content/20230313_gene_catalog.ffn'
# total_genes = 0
# genes_per_sample = Counter()
# gene_per_host = Counter()
# all_genes = SeqIO.parse(gene_catalog, 'fasta')
# for gene in all_genes:
#     total_genes += 1
#     genes_per_sample[gene.id.split('_NODE')[0]] += 1
#     gene_per_host[gene.id[0]] += 1
# # total_genes
# # 4,911,821
# len([kegg_dict[x] for x in kegg_dict if kegg_dict[x] != ''])
# # 2,772,118
# len([dram_cazy_dict[x] for x in dram_cazy_dict if dram_cazy_dict[x] != ''])
# # 122,164
# len(cazy_dict) # from cayman
# # 150,832
# len([pfam_dict[x] for x in pfam_dict if pfam_dict[x] != ''])
# # 259,590
# binned_scaffolds = set(all_scaffolds_info['scaffold'])
# in_mag = 0
# all_genes = SeqIO.parse(gene_catalog, 'fasta')
# for gene in all_genes:
#     if scaffold_from_geneid(gene.id) in binned_scaffolds:
#         in_mag += 1
# # 3,678,762
# in_og = len(gene_og_dict)
# # 2,511,732

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
# plt.savefig('results/figures/visualize_temp/cluster_size_hist.png')
# plt.close()

# cluster_host_dict = {}
# for cluster in clusters:
#     if len(clusters[cluster].genes) == 1:
#         continue
#     hosts_seen = Counter([x[0] for x in clusters[cluster].genes])
#     cluster_host_dict[cluster] = len(hosts_seen.keys())
# plt.hist(cluster_host_dict.values(), bins = 5)
# plt.savefig('results/figures/visualize_temp/cluster_host_hist.png')
# plt.close()

# do the above but using mapping result to tell if a cluster has been detected
'''
A cluster is considered detected if any gene in the cluster is found at a coverage > 0
and breadth > 0.75 in the sample
To find this, for each sample, 
'''


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