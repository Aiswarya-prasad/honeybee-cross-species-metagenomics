
import os
import glob
import pandas as pd

header = False
for file in glob.glob("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/*"):
    with open(file, "r") as file_fh:
        if file.endswith("SNVs.tsv"):
            print(f"{file}")
            for line in file_fh:
                line = line.strip()
                if not header:
                    header = line
                else:
                    print(header)
                    print(line)
                    break
snvs_df = pd.read_csv("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/C1.4_profile.IS_SNVs.tsv", sep = "\t")
num_SNVs = {}
for i, gene in enumerate(snvs_df["gene"]):
    alleles = snvs_df["allele_count"].iloc[i]
    class_name = snvs_df["class"].iloc[i]
    if "SNV" in class_name:
        if gene in num_SNVs.keys():
            num_SNVs[gene] += 1
        else:
            num_SNVs[gene] = 1

genes_df = pd.read_csv("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/C1.4_profile.IS_SNVs.tsv", sep = "\t")

glob.glob("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/*")

'''
Gene info is mentioned for each sample the instrain results only include genes from representative MAGs that reads of the sample happened to map to
If a certain pattern is of interest (eg. dN/dS ratio), read the instrain output of the given sample to find the genes of interest and then use the gene_info.tsv file to get the gene information
'''

rep_mag_scaffolds = []
with open('results/10_instrain/00_prepare_mags/mag_rep_database.fa', 'r') as fh:
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            rep_mag_scaffolds.append(line.split('>')[1])
rep_mags = pd.read_csv('results/10_instrain/00_prepare_mags/scaffold_to_bin_file.tsv', sep='\t', header=None)[1].unique()
# get gene info for all rep mag genes
genes_of_interest = []
mag_annotations_list = glob.glob('results/09_MAGs_collection/dram_output/*/annotations.tsv')
for file in mag_annotations_list:
    with open(file, 'r') as fh:
        for line in fh:
            line = line.strip()
            line = line.split('\t')
            mag = '_'.join(line[0].split('_')[0:2])
            if mag in rep_mags:
                ko = line[3]
                if ko == 'K21572':
                    gene = line[0].split(f'{mag}_')[1]
                    genes_of_interest.append(gene)
                    print(gene)

samples_to_exclude = ["F4-5", "F5-1", "M6-2", "D9-5", "F7-5"]
samples = [x.split('/')[-2] for x in glob.glob('results/10_instrain/02_instrain_profile/*/') if x not in samples_to_exclude]
for sample in samples:
    instrain_gene_info = pd.read_csv(f'results/10_instrain/02_instrain_profile/{sample}/output/{sample}_gene_info.tsv', sep='\t')
    instrain_gene_info[instrain_gene_info['gene'].isin(genes_of_interest)]['dNdS_substitutions']

instrain_gene_info[instrain_gene_info['dNdS_substitutions'] > 1]['gene']

# identify genes in D2-1 with dN/dS ratio > 1

# # SusC
# sus_genes = [x for x in gene_info_df[gene_info_df['ko'] == 'K21573']['gene']]
# # SusD
# sus_genes += [x for x in gene_info_df[gene_info_df['ko'] == 'K21572']['gene']]

# find these strains in instrain output
instrain_gene_info[instrain_gene_info['gene'].isin(sus_genes)]