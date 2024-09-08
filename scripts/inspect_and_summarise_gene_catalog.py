import os
import sys
import glob
import pickle
import re
from Bio import SeqIO
from itertools import groupby
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

'''
This script is used to inspect and summarise the gene catalog.
It can be run after kraken, kaiju and cd hit clustering on the
gene catalog are completed.
It creates a summary table and plots the results into
results/figures/....png
It does not need any input files and can be run as follows:
python3 scripts/inspect_and_summarise_gene_catalog.py
as the result paths are hard coded in here
At the moment it was run by hand, but might be integrated into the pipeline later
'''

catalog_ffn = 'results/08_gene_content/20230313_gene_catalog.ffn'
catalog_faa = 'results/08_gene_content/20230313_gene_catalog.faa'
cdhit_clustering = 'results/08_gene_content/00_cdhit_clustering/20230313_gene_catalog_cdhit9590.fasta.clstr'
# cdhit_clusters = 'results/08_gene_content/00_cdhit_clustering/cluster_host_affiliations.tsv'
kraken_out = 'results/08_gene_content/04_kraken2_on_genes/20230313_gene_catalog.kraken'
kraken_rep = 'results/08_gene_content/04_kraken2_on_genes/20230313_gene_catalog_report.txt'
kaiju_out = 'results/08_gene_content/04_kaiju_on_genes/nr/20230313_gene_catalog_taxa_full.txt'
mag_metadata = '../backups/09_MAGs_collection/MAGs_metadata_summary.tsv'
if not os.path.exists('mkdir -p results/08_gene_content/05_parse_catalog'):
    os.system('mkdir -p results/08_gene_content/05_parse_catalog')
kaiju_out_fixed = 'results/08_gene_content/05_parse_catalog/20230313_gene_catalog_taxa_full_edited.txt'

def sample_of_gene(gene_id):
    return gene_id.split('_NODE')[0]

def host_of_gene(gene_id):
    sample = sample_of_gene(gene_id)
    if sample.startswith('M'):
        return 'Apis mellifera'
    if sample.startswith('C'):
        return 'Apis cerana'
    if sample.startswith('D'):
        return 'Apis dorsata'
    if sample.startswith('F'):
        return 'Apis florea'
    if sample.startswith('A'):
        return 'Apis andreniformis'

def contig_of_gene(gene_id):
    return '_'.join(gene_id.split('_')[1:-1])

global mag_dict
if not os.path.exists('results/08_gene_content/05_parse_catalog/scaffolds_to_bin.tsv'):
    with open('results/08_gene_content/05_parse_catalog/scaffolds_to_bin.tsv', 'w+') as out_fh:
        for mag in glob.glob('results/09_MAGs_collection/MAGs/*.fa'):
            with open(mag, "r") as m:
                for line in m:
                    mag_name = os.path.basename(mag).split(".")[0]
                    if line.startswith(">"):
                        scaffold = line.strip().split(">")[1]
                        out_fh.write(f"{scaffold}\t{mag_name}\n")
else:
    mag_dict_df = pd.read_csv('results/08_gene_content/05_parse_catalog/scaffolds_to_bin.tsv', sep='\t', header=None)
    mag_dict = mag_dict_df.set_index(0).to_dict()[1]

def get_mag_of_gene(gene_id):
    sample = sample_of_gene(gene_id)
    contig = contig_of_gene(gene_id)
    try:
        mag = mag_dict[f'{sample}_{contig}']
    except KeyError:
        mag = 'unbinned'
    return mag

# global mag_matadata_df
# mag_matadata_df = pd.read_csv(mag_metadata, sep='\t')
# # mag_matadata_df.columns
# def gtdb_of_mag(mag):
#     species = str(mag_matadata_df[mag_matadata_df['ID'] == mag]['Genus'].values[0])
#     genus = str(mag_matadata_df[mag_matadata_df['ID'] == mag]['Species'].values[0])
#     if genus == 'nan':
#         return 'NA'
#     if species == 'nan':
#         return genus
#     return species

global gtdb_dict
gtdb_dict = {}
with open('results/09_MAGs_collection/gtdb_output/classify/20230313_MAGs.bac120.summary.tsv', 'r') as fh:
    header = fh.readline()
    for line in fh:
        mag = line.split()[0]
        tax_line = line.split()[1]
        if tax_line.split(';')[-1] == 's__':
            gtdb_dict[mag] = tax_line.split(';')[-2]
        else:
            gtdb_dict[mag] = tax_line.split(';')[-1]

def gtdb_of_mag(mag):
    if mag in gtdb_dict.keys():
        return gtdb_dict[mag]
    else:
        return 'unknown'


if not os.path.exists('results/08_gene_content/05_parse_catalog/kraken_df.pickle'):
    kraken_df = pd.read_csv(kraken_out, sep='\t', header=None)
    kraken_df.columns = ['status', 'gene_id', 'taxon_name', 'length', 'taxon_depth']
    kraken_df['ncbi_id'] = kraken_df['taxon_name'].str.split(' \(taxid ').str[1].str.split('\)').str[0]
    kraken_df['taxon_name'] = kraken_df['taxon_name'].str.split(' \(taxid ').str[0]
    pd.to_pickle(kraken_df, 'results/08_gene_content/05_parse_catalog/kraken_df.pickle')
else:
    kraken_df = pd.read_pickle('results/08_gene_content/05_parse_catalog/kraken_df.pickle')

if os.path.exists('results/08_gene_content/05_parse_catalog/kraken_df_dict.pickle'):
    kraken_df_dict = pd.read_pickle('results/08_gene_content/05_parse_catalog/kraken_df_dict.pickle')
else:
    kraken_df_dict = kraken_df.set_index('gene_id').to_dict()
    pd.to_pickle(kraken_df_dict, 'results/08_gene_content/05_parse_catalog/kraken_df_dict.pickle')

def get_kraken_tax_of_gene(gene_id):
    try:
        tax = kraken_df_dict['taxon_name'][gene_id]
        if tax == 'unclassified':
            return 'unknown'
        return tax
    except:
        return 'unknown'

if not os.path.exists('results/08_gene_content/05_parse_catalog/kaiju_df.pickle'):
    kaiju_df = pd.DataFrame()
    # cannot read unequal lines so take out lines starting with 'U'
    with open(kaiju_out, 'r') as kaiju_fh:
        if not os.path.exists(kaiju_out_fixed):
            with open(kaiju_out_fixed, 'w+') as fixed_fh:
                for line in kaiju_fh:
                    if line.startswith('U'):
                        line = line.strip()
                        # unclassified lines have only 3 columns
                        # extent to 8
                        line = '\t'.join([line]+['']*5)+'\n'
                    fixed_fh.write(line)
    kaiju_df = pd.read_csv(kaiju_out_fixed, sep='\t', header=None)
    kaiju_df.columns = ['status', 'gene_id', 'taxon_id', 'best_match_score', 'best_match_taxa', 'best_match_accessions', 'best_match_sequences', 'taxon_name']                
    pd.to_pickle(kaiju_df, 'results/08_gene_content/05_parse_catalog/kaiju_df.pickle')
else:
    kaiju_df = pd.read_pickle('results/08_gene_content/05_parse_catalog/kaiju_df.pickle')

if os.path.exists('results/08_gene_content/05_parse_catalog/kaiju_df_dict.pickle'):
    kaiju_df_dict = pd.read_pickle('results/08_gene_content/05_parse_catalog/kaiju_df_dict.pickle')
else:
    kaiju_df_dict = kaiju_df.set_index('gene_id').to_dict()
    pd.to_pickle(kaiju_df_dict, 'results/08_gene_content/05_parse_catalog/kaiju_df_dict.pickle')

def get_kaiju_tax_of_gene(gene_id):
    try:
        tax = kaiju_df_dict['taxon_name'][gene_id]
        if str(tax) == 'nan':
            return 'unknown'
        return tax
    except:
        return 'unknown'

pickle_file = cdhit_clustering + ".pickle"

if os.path.exists(pickle_file):
    print(f"loading pickle file from {pickle_file}")
    with open(pickle_file, "rb") as pkl_fh:
        cluster_groups = pickle.load(pkl_fh)
else:
    with open(cdhit_clustering, "r") as cfh:
            cluster_groups = [
                list(group)
                for key, group in groupby(cfh, lambda line: line.startswith(">Cluster"))
                if not key
            ]
    with open(pickle_file, "wb") as pkl_fh:
        pickle.dump(cluster_groups, pkl_fh)

global cluster_genes_dict
global gene_clusters_dict
total_clusters = len(cluster_groups)
if os.path.exists('results/08_gene_content/05_parse_catalog/cluster_gene_dict.pickle') and os.path.exists('results/08_gene_content/05_parse_catalog/gene_cluster_dict.pickle'):
    cluster_genes_dict = pd.read_pickle('results/08_gene_content/05_parse_catalog/cluster_gene_dict.pickle')
    gene_clusters_dict = pd.read_pickle('results/08_gene_content/05_parse_catalog/gene_cluster_dict.pickle')
else:
    cluster_genes_dict = {}
    gene_clusters_dict = {}
    for i, genes in enumerate(cluster_groups, start=1):
        print(f"Processing cluster {i} of {total_clusters}", end="\r")
        cluster_num = f'cluster_{i}'
        ids = [x.split("...")[0].split(">")[1] for x in [g.split()[2] for g in genes]]
        rep_id = [g.split()[2].split("...")[0].split(">")[1] for g in genes if g.strip().endswith("*")].pop()
        cluster_genes_dict[cluster_num] = ids
        for gene in ids:
            if gene in gene_clusters_dict.keys():
                print(f'{gene} was already in keys')
            else:
                gene_clusters_dict[gene] = cluster_num
        pd.to_pickle(cluster_genes_dict, 'results/08_gene_content/05_parse_catalog/cluster_gene_dict.pickle')
        pd.to_pickle(gene_clusters_dict, 'results/08_gene_content/05_parse_catalog/gene_cluster_dict.pickle')
        # outstring = "\t".join([str(x) for x in cluster_info.values()])

def kraken_contig_tax(gene_id):
    sample = sample_of_gene(gene_id)
    contig = contig_of_gene(gene_id)
    kraken_file = f'results/05_assembly/contig_fates/kraken2/{sample}.kraken'
    kraken_contigs_df = pd.read_csv(f'results/05_assembly/contig_fates/kraken2/{sample}.kraken', header=None, sep='\t')
    if list(kraken_contigs_df.loc[kraken_contigs_df[1] == contig][0])[0] == 'U':
        return 'unknown'
    else:
        tax = list(kraken_contigs_df.loc[kraken_contigs_df[1] == contig][2])[0].split(' (')[0]
        return tax

num_genes = len(gene_clusters_dict)
phage_count = 0
both_unclassified = 0
binned_not_classified = 0
one_classified = 0
only_one_classified = 0
both_classifed = 0
apis_gene = 0
homo_sapiens_gene = 0
binned_in_mag = 0
filtered_records = []
phage_records = []
total_records_read = 0
'''
    Considerations:
        Considering contig kraken result, _kraken and kaiju for each gene and the binning result + gtdb annotation
        sometimes both kraken and kaiju are unknown but it is binned with gtdb annotation (mostly matching kraken at contig level)
        so excluding those with both unknown does not make sense.
        
    Approach for filtering:
        First check if the gene is binned. 
            If binned, (binned_in_mag) include it unless the gtdb annotation is unclassified or unknown meaning
                the mag is not on the bac120 summary made by gtdb. (binned_not_classified)
            If it is not binned, check the gene-level taxonomy from kaiju and kraken.
                If both are unknown, exclude it. (both_unclassified)
                If only one of them is unknown,
                    (one_classified)
                    exclude it if the kraken calls it Apis mellifera (apis_gene)
                    otherwise include it
                    If both are classified,(both_classifed) count those but these will not be included if
                        kraken calls it Apis mellifera or Homo sapiens
                        else it means only one is classified (not unknown) include these!
'''
for i, rec in enumerate(SeqIO.parse(catalog_ffn, 'fasta')):
    include_record = False
    print(f'{i}/{num_genes}', end = '\r')
    total_records_read +=1
    gene_id = rec.id
    kaiju_taxon = get_kaiju_tax_of_gene(gene_id)
    kraken_taxon = get_kraken_tax_of_gene(gene_id)
    if kraken_taxon == 'Archaea':
        phage_count += 1
        phage_records.append(rec)
    if get_mag_of_gene(gene_id) not in ['unbinned'] and 'unbinned' not in get_mag_of_gene(gene_id):
        binned_in_mag += 1
        include_record = True
        if 'Apis' in kraken_taxon and gtdb_of_mag(get_mag_of_gene(gene_id)) in ['Unclassified', 'unknown', 'unclassified']:
            binned_not_classified += 1
            include_record = False
    else:
        if kaiju_taxon == 'unknown' and kraken_taxon == 'unknown':
            both_unclassified +=1
            include_record = False
        else:
            include_record = True
            # both kraken or kaiju are classified
            # ie at least one classified
            one_classified += 1
            if 'Apis mellifera' in kraken_taxon:
                # kraken is classified, kaiju could be unknown or classified
                # kaiju not considered here becuase this is already not a binned
                # gene and kaiju shows quite a few false positives
                apis_gene += 1
                include_record = False
            if 'Homo sapiens' in kraken_taxon:
                # kraken is classified, kaiju could be unknown or classified
                # kaiju not considered here becuase this is already not a binned
                # gene and kaiju shows quite a few false positives
                homo_sapiens_gene += 1
                include_record = False
            else:
                if kaiju_taxon != 'unknown' and kraken_taxon != 'unknown':
                    both_classifed += 1
                    include_record = True
                else:
                    # one of them is classified and kraken is not calling it Apis mellifera
                    # print(f'gtdb: {gtdb_of_mag(get_mag_of_gene(gene_id))}\nKaiju: {kaiju_taxon}\nKraken: {kraken_taxon}\n')
                    # if super strict, exclude these too as these have only one of kraken or kaiju classifying them 
                    # and the other is 'unknown'
                    # If more relaxed, include these as they are classified by one of them - This is often
                    # not necessarily the same species but still likely
                    # one of the core genera or even related taxa
                    # so being strict here may not be ideal...
                    # Not being strict:
                    only_one_classified += 1
                    include_record = True
                    # if kraken_taxon == 'unknown':
                    #     print(f'Kaiju: {kaiju_taxon}\nKraken: {kraken_taxon}\n')
    if include_record:
        filtered_records.append(rec)

# just check to see that the logic for filtering doesn't make
# any duplicates!
# seen = set()
# for rec in filtered_records:
#     if rec.id in seen:
#         print(f'{rec.id} is repeated')
#     seen.add(rec.id)


print(f'total records read = {total_records_read:,}')
n = len(filtered_records)
print(f'filtered records = {n:,}')
print(f'kept after filtering % = {len(filtered_records)/num_genes*100}')
n = binned_in_mag
print(f'genes in mags = {n:,}')
print(f'genes in mags % = {(binned_in_mag)/num_genes*100}')
print(f'genes annotated as apis mellifera (not binned and classified) = {apis_gene:,}')
print(f'genes annotated as homo sapiens (not binned and classified) = {homo_sapiens_gene:,}')
print(f'unclassified by both = {both_unclassified:,}')
print(f'classified by both (excludes Apis mellifera and Homo sapiens) = {both_classifed:,}')
print(f'classified by atleast one of them (excludes Apis mellifera and Homo sapiens) = {one_classified:,}')
print(f'classified by only one of them (excludes Apis mellifera and Homo sapiens) = {only_one_classified:,}')

print('Writing filtered records to file...')
with open('results/08_gene_content/05_parse_catalog/20230313_gene_catalog_filtered.ffn', 'w+') as fh:
    SeqIO.write(filtered_records, fh, 'fasta')

# prodigal_ffn = 'results/06_metagenomicORFs/D1-1/prodigal_out/D1-1.ffn'
# dict_name_pos = {}
# for rec in SeqIO.parse(prodigal_ffn, 'fasta'):
#         gene_name = rec.id
#         node = '_'.join(rec.id.split('_')[0:-1])
#         gene_start = rec.description.split(' # ')[1]
#         gene_end = rec.description.split(' # ')[2]
#         dict_name_pos.update({f'{node}:{gene_start}-{gene_end}': gene_name})

# dict_genename_in_orthofiles = {}
# for or_rec, re_rec in zip(SeqIO.parse(original, 'fasta'), SeqIO.parse(renamed, 'fasta')):
#     dict_genename_in_orthofiles.update({dict_name_pos[or_rec.id]: re_rec.id})

# # use these for making
# dict_name_pos = {}
# prodigal_ffn = 'results/06_metagenomicORFs/{sample_of_mag}/prodigal_out/{sample_of_mag}.ffn'
# # and
# dict_genename_in_orthofiles = {}
# f'results/11_phylogenies/00_genomes/{mag}/{mag}_original.ffn'
# f'results/11_phylogenies/00_genomes/{mag}/{mag}.ffn'
# # as seen above
# # finish with dict indexed by gene name (as in gene catalog)
# # mapping to MAG_num which can be used to parse orthofinder results
# # To get ortholog information for genes,
# # make a dict with index as gene name (format MAG_num) and the OG it corresponds to
# # using orthofinder information
# 'results/11_phylogenies/02_orthofinder_results/*/Results_*/Orthogroups/Orthogroups.tsv'
# # use this to get the OGs that as scp
# 'results/11_phylogenies/02_orthofinder_results/*/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt'
# # finally using these, get the OG and coreness of every gene!