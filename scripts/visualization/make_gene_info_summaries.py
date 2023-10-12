import os
import json
import gzip
from Bio import SeqIO
# from Bio.KEGG import REST
# from Bio.KEGG.KGML import KGML_parser
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

def sample_of_gene(gene_id):
    return gene_id.split('_NODE')[0]

def contig_of_gene(gene_id):
    return '_'.join(gene_id.split('_')[1:-1])

def scaffold_from_geneid(gene_id):
    return '_'.join(gene_id.split('_')[:-1])

# get the mag the gene is in
all_scaffolds_info = pd.read_csv('results/09_MAGs_collection/all_scaffold_to_bin.tsv', sep='\t')
# get the mag info for the scaffold
def get_mag_of_gene(gene):
    try:
        all_scaffolds_info[all_scaffolds_info['scaffold'] == scaffold_from_geneid(gene)]['mag'].values[0]
    except:
        'unbinned'

samples = [os.path.basename(x).split('.')[0] for x in glob.glob('results/08_gene_content/06_cayman/*cazy.txt.gz')]
# get the gene counts
gene_counts = {}
for sample in samples:
    gene_counts[sample] = pd.read_csv(f'results/08_gene_content/06_cayman/{sample}.gene_counts.txt.gz', compression='gzip', header=0, sep='\t')
def get_gene_counts(gene, sample):
    gene_counts[sample][gene_counts[sample]['gene'] == gene]['uniq_rpkm'].values[0]

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


# make an iterator for all the genes from the catalog headers
gene_catalog = 'results/08_gene_content/20230313_gene_catalog.ffn'
genes = SeqIO.parse(gene_catalog, 'fasta')

total_genes = 0
genes_per_sample = Counter()
gene_per_host = Counter()
for gene in genes:
    total_genes += 1
    genes_per_sample[gene.id.split('_NODE')[0]] += 1
    gene_per_host[gene.id[0]] += 1

print(f'Total genes: {total_genes}')
print(f'Genes per host: {gene_per_host}')
print(f'Genes per sample: {genes_per_sample}')


# for this we need to find out if the gene is in a MAG and whether the mag is in the OG analysis and then what the OG is
# get the MAGs that are in the OG analysis
phylo_metadata = pd.read_csv('results/11_phylogenies/phylo_genomes_metadata.tsv', sep='\t')
og_mags = set(phylo_metadata['ID'])

# retrive og_gene_id for each gene_id in the gene catalog
gene_og_gene_id_dict = {}
og_gene_id_dict = {}
for group in phylo_metadata['Phylogeny_group'].unique():
    gene_og_gene_id_dict[group] = {}
    og_gene_id_dict[group] = {}
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
                    coreness = line.split('\t')[1].split(';')
                    genes = line.split('\t')[5].strip().split(';')
                    for gene in genes:
                        gene_og_gene_id_dict[group][gene] = og
                        if og in og_gene_id_dict[group]:
                            og_gene_id_dict[group][og].add(gene)
                        else:
                            og_gene_id_dict[group][og] = set(gene)

# pickle.dump(gene_og_gene_id_dict, open('results/08_gene_content/07_OG_coreness/gene_og_gene_id_dict.p', 'wb'))
# pickle.dump(og_gene_id_dict, open('results/08_gene_content/07_OG_coreness/og_gene_id_dict.p', 'wb'))



# for genus in [os.path.basename('/'.join(x.split('/')[:-1])) for x in glob.glob('results/11_phylogenies/02_orthofinder_results/*/')]:
#     '''
#     IMPLEMENT more checks or a better way to get genera!
#     '''
#     if '.bak' in genus:
#         continue
#     input_orthofile = f'results/11_phylogenies/02_orthofinder_results/{genus}/Results_{genus}/Orthogroups/Orthogroups.GeneCount.tsv'
#     output_cogfile = f'results/08_gene_content/07_OG_coreness/{genus}_motupan_input.tsv'
#     if os.path.exists(output_cogfile):
#         print(f'{output_cogfile} already exists, skipping')

#     ortho_df = pd.read_csv(input_orthofile, sep = '\t')
#     ortho_df = ortho_df.set_index('Orthogroup')
#     ortho_df.index.names = [None]
#     cog_df = {}
#     for genome in ortho_df.keys():
#         if genome == 'Total':
#             continue
#         cog_df[genome] = [x for x in ortho_df.index[ortho_df[genome] > 0]]

#     with open(output_cogfile, 'w') as out_fh:
#         for genome in cog_df:
#             out_list = '\t'.join(cog_df[genome])
#             out_fh.write(f'{genome}\t{out_list}\n')
# 
'''
ran motupan using,
genus="g__Gilliamella"
genus="g__Frischella"
genus="g__Entomomonas"
genus="g__Dysgonomonas"
genus="g__Lactobacillus"
genus="g__Bombilactobacillus"
genus="g__Snodgrassella"
genus="g__Bifidobacterium"
genus="g__Q3-R57-64"
genus="g__Pectinatus"
genus="g__Bartonella_A"
genus="g__Apibacter"
genus="g__WRHT01"
genus="g__CALYQJ01"
genus="g__Enterobacter"
genus="g__Melissococcus"
genus="g__Spiroplasma"
genus="g__Saezia"
genus="g__Pluralibacter"
genus="g__CALYQQ01"
genus="g__Commensalibacter"
genus="g__Klebsiella"
genus="g__Apilactobacillus"


output_cogfile="results/08_gene_content/07_OG_coreness/${genus}_motupan_input.tsv"
checkm_info="results/09_MAGs_collection/checkm_merged.tsv"
outfile="results/08_gene_content/07_OG_coreness/${genus}_motupan_output.tsv"
log="results/08_gene_content/07_OG_coreness/${genus}_motupan.log"
mOTUpan.py --gene_clusters_file ${output_cogfile} --boots 10 -o ${outfile} --checkm ${checkm_info} | tee ${log}
'''

# plot cazy profiling output
# 'results/08_gene_content/06_cayman/C5-2.gene_counts.txt.gz'
# 'results/08_gene_content/06_cayman/C5-2.cazy.txt.gz'
# read gz file

samples = [os.path.basename(x).split('.')[0] for x in glob.glob('results/08_gene_content/06_cayman/*cazy.txt.gz')]
# get the gene counts
cazy_counts = {}
for sample in samples:
    cazy_counts[sample] = pd.read_csv(f'results/08_gene_content/06_cayman/{sample}.cazy.txt.gz', compression='gzip', header=0, sep='\t')
    # remove into a different dict the features, filtered_reads and total_reads
    cazy_counts[sample] = cazy_counts[sample][cazy_counts[sample]['feature'] != 'filtered_reads']
    cazy_counts[sample] = cazy_counts[sample][cazy_counts[sample]['feature'] != 'total_reads']

# combine cazy df into one melted df for sample, feature, uniq_rpkm
cazy_df = pd.DataFrame()
for sample in samples:
    cazy_df_sample = cazy_counts[sample].loc[:, ['feature', 'uniq_rpkm']]
    cazy_df_sample['sample'] = sample
    cazy_df = cazy_df._append(cazy_df_sample)
cazy_df.to_csv('results/figures/cazy_df.tsv', sep='\t')