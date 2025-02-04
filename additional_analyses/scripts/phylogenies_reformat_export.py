import ete3
import numpy as np
import pandas as pd
from skbio import TreeNode
from skbio.stats.evolve import hommola_cospeciation
import skbio as sk
import pickle
import seaborn as sns
from io import StringIO
from os.path import join, basename, exists
from collections import OrderedDict
from itertools import combinations
from scipy import stats
from os import makedirs, system, getcwd
from shutil import copytree
from glob import glob
import matplotlib.pyplot as plt
from Bio import Phylo

phylo_metadata = pd.read_csv('results/11_phylogenies/phylo_genomes_metadata.tsv', sep='\t')
handmade_spec_names = pd.read_csv('results/figures/handmade_species_names.csv')
handmade_spec_names['magotu'] = handmade_spec_names['cluster'].astype(str)
phylo_metadata = phylo_metadata.merge(handmade_spec_names[['magotu', 'MAG_species_name_final', 'MAG_species_name_final_nospace']], left_on='magOTU', right_on='magotu', how='left')
phylo_metadata['MAG_species_name_final'] = phylo_metadata['MAG_species_name_final'].fillna(phylo_metadata['Species'])
phylo_metadata['MAG_species_name_final_nospace'] = phylo_metadata['MAG_species_name_final_nospace'].fillna(phylo_metadata['Species'])
mag_info = pd.read_csv('results/09_MAGs_collection/MAGs_metadata_summary.tsv', sep = '\t')

def get_gtdb_taxonomy():
    if not exists('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-bac120/bac120_taxonomy_r214.tsv'):
        print('Downloading GTDB taxonomy')
        system('wget -O results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-bac120/bac120_taxonomy_r214.tsv.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.0/bac120_taxonomy_r214.tsv.gz')
        system('gunzip results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-bac120/bac120_taxonomy_r214.tsv.gz')
        # read, parse tax string and pickle it
        gtdb_taxonomy = pd.read_csv('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-bac120/bac120_taxonomy_r214.tsv', sep='\t', header=None)
        # split the tax string in column 1 into domain, phylum, class, order, family, genus, species
        # using d__ prefix for domain, p__ for phylum, c__ for class, o__ for order, f__ for family, g__ for genus, s__ for species
        gtdb_taxonomy['domain'] = gtdb_taxonomy[1].apply(lambda x: x.split(';')[0].split('d__')[1])
        gtdb_taxonomy['phylum'] = gtdb_taxonomy[1].apply(lambda x: x.split(';')[1].split('p__')[1])
        gtdb_taxonomy['class'] = gtdb_taxonomy[1].apply(lambda x: x.split(';')[2].split('c__')[1])
        gtdb_taxonomy['order'] = gtdb_taxonomy[1].apply(lambda x: x.split(';')[3].split('o__')[1])
        gtdb_taxonomy['family'] = gtdb_taxonomy[1].apply(lambda x: x.split(';')[4].split('f__')[1])
        gtdb_taxonomy['genus'] = gtdb_taxonomy[1].apply(lambda x: x.split(';')[5].split('g__')[1])
        gtdb_taxonomy['species'] = gtdb_taxonomy[1].apply(lambda x: x.split(';')[6].split('s__')[1])
        gtdb_taxonomy.to_pickle('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-bac120/bac120_taxonomy_r214.pkl')
    if 'gtdb_taxonomy_downoaded' in globals():
        pass
    else:
        global gtdb_taxonomy_downoaded
        gtdb_taxonomy_downoaded = pd.read_pickle('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-bac120/bac120_taxonomy_r214.pkl')

def get_tax_info(ID):
    get_gtdb_taxonomy()
    if ID in set(gtdb_taxonomy_downoaded[0]):
        # return a dict
        return gtdb_taxonomy_downoaded[gtdb_taxonomy_downoaded[0] == ID].to_dict(orient='records')[0]
    else:
        return None

# Bombilactobacillus

bombi_isolate_names = {'GCF_000970795.1': 'Bombilactobacillus mellifer (strain=Bin4)',
                        'GCF_025290075.1': 'Bombilactobacillus mellifer (strain=MRS2-bin.14)',
                        'GCF_025290985.1': 'Bombilactobacillus mellifer (strain=MRS1-bin.8)',
                        'GCF_025291455.1': 'Bombilactobacillus mellifer (strain=GUT-bin.11)',
                        'GCF_042663165.1': 'Bombilactobacillus mellifer (strain=CCUG 57507)',
                        'GCF_000967245.1': 'Bombilactobacillus mellis (strain=Hon2)',
                        'GCF_013345055.1': 'Bombilactobacillus mellis (strain=ESL0295)',
                        'GCF_013346925.1': 'Bombilactobacillus mellis (strain=ESL0294)',
                        'GCF_013347085.1': 'Bombilactobacillus mellis (strain=ESL0394)',
                        'GCF_026229285.1': 'Bombilactobacillus mellis (strain=LB26)',
                        'GCF_042676405.1': 'Bombilactobacillus mellis (strain=CCUG 63289)',
                        'GCF_003515755.1': 'Bombilactobacillus bombi (strain=LV-8.1)',
                        'GCF_003515805.1': 'Bombilactobacillus bombi (strain=BI-1.1)',
                        'GCF_003522965.1': 'Bombilactobacillus bombi (strain=BI-2.5)',
                        'GCF_013607485.1': 'Bombilactobacillus bombi (strain=XV6)',
                        'GCF_042432605.1': 'Bombilactobacillus bombi (strain=CCM 8440)',
                        'GCF_013385145.1': 'Bombilactobacillus apium (strain=DCY120)',
                        'GCF_023380265.1': 'Bombilactobacillus folatiphilus (strain=SG4_D2)',
                        'GCF_023380245.1': 'Bombilactobacillus thymidiniphilus (strain=SG4_A1)'
}

recomputed_tree = ete3.Tree('additional_analyses/results/add_isolates/g__Bombilactobacillus/Results_Nov18/Species_Tree/SpeciesTree_rooted.txt')
recomputed_tree.get_leaves()
# rename tip labels by adding the species name instead of the acession ID or MAG ID
for leaf in recomputed_tree:
    print(leaf.name)
    if leaf.name in bombi_isolate_names:
        leaf.name = bombi_isolate_names[leaf.name]
    else:
        species_name = phylo_metadata[phylo_metadata['ID'] == leaf.name].iloc[0]['MAG_species_name_final']
        if 'GCA' in leaf.name:
            leaf.name = f'{species_name}'
        else:
            leaf.name = f'{leaf.name}--{species_name}'
recomputed_tree.get_leaf_names()
recomputed_tree.write(outfile='additional_analyses/results/add_isolates/g__Bombilactobacillus/SpeciesTree_renamed.txt', format=1)
# mark mags and isolates and decorate tree to save as pdf (not using render because it's not working)
recomputed_tree = ete3.Tree('additional_analyses/results/add_isolates/g__Bombilactobacillus/SpeciesTree_renamed.txt')
recomputed_tree.get_leaves()
for leaf in recomputed_tree:
    if '--' in leaf.name:
        leaf.add_feature('color', 'red')
    else:
        leaf.add_feature('color', 'blue')
ts = ete3.TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_branch_support = False
ts.show_scale = False
ts.show_scale = False
ts.show_border = False
ts.show_branch_support = False
ts.show_branch_length = False
ts.show_leaf_name = False

def layout(node):
    if node.is_leaf():
        F = ete3.TextFace(node.name, fgcolor=node.color)
        ete3.add_face_to_node(F, node, column=0, position="branch-right")
ts.layout_fn = layout
recomputed_tree.render('additional_analyses/results/add_isolates/g__Bombilactobacillus/SpeciesTree_renamed.pdf', tree_style=ts)

# Lactobacillus

lacto_isolate_names = {'GCF_002837055.1': 'Lactobacillus apis (strain=LMG 26964)',
                       'GCF_900112665.1': 'Lactobacillus bombicola (strain=R-53102)',
                       'GCF_026428255.1': 'Lactobacillus helsingborgensis (strain=IBH002)',
                       'GCF_019972815.1': 'Lactobacillus huangpiensis (strain=F306-1)',
                       'GCF_000970755.1': 'Lactobacillus kimbladii (strain=Hma2)',
                       'GCF_000967195.1': 'Lactobacillus kullabergensis (strain=Biut2)',
                       'GCF_019972835.1': 'Lactobacillus laiwuensis (strain=F551-2)',
                       'GCF_026185255.2': 'Lactobacillus melliventris (strain=IBH004)',
                       'GCF_003693045.1': 'Lactobacillus melliventris (strain=ESL0259)',
                       'GCF_019469265.1': 'Lactobacillus panisapium (strain=ESL0416)',
                       'GCF_016100975.1': 'Lactobacillus sp016100975 (strain=W8174)',
                       'GCF_014323605.1': 'Lactobacillus kimbladii (strain=Dan47)'
}

recomputed_tree = ete3.Tree('additional_analyses/results/add_isolates/g__Lactobacillus/Results_Nov23/Species_Tree/SpeciesTree_rooted.txt')
recomputed_tree.get_leaves()
# rename tip labels by adding the species name instead of the acession ID or MAG ID
for leaf in recomputed_tree:
    print(leaf.name)
    if leaf.name in lacto_isolate_names:
        leaf.name = lacto_isolate_names[leaf.name]
    else:
        species_name = phylo_metadata[phylo_metadata['ID'] == leaf.name].iloc[0]['MAG_species_name_final']
        if 'GCA' in leaf.name:
            leaf.name = f'{species_name}'
        else:
            leaf.name = f'{leaf.name}--{species_name}'
recomputed_tree.get_leaf_names()
recomputed_tree.write(outfile='additional_analyses/results/add_isolates/g__Lactobacillus/SpeciesTree_renamed.txt', format=1)


# Bifidobacterium


bifido_isolate_names = {'GCF_007559275.1': 'Bifidobacterium apousia (strain=W8102)',
                        'GCF_002715865.1': 'Bifidobacterium asteroides (strain=DSM 20089)',
                        'GCF_003202755.1': 'Bifidobacterium asteroides_F (strain=ESL0199)',
                        'GCF_003202695.1': 'Bifidobacterium asteroides_G (strain=ESL0200)',
                        'GCF_009683175.1': 'Bifidobacterium asteroides_H (strain=VRA_9sq_n)',
                        'GCF_019469425.1': 'Bifidobacterium asteroides_I (strain=ESL0447)',
                        'GCF_016102005.1': 'Bifidobacterium choladohabitans (strain=B14384H11)',
                        'GCF_000706765.1': 'Bifidobacterium indicum (strain=strain=LMG 11587)',
                        'GCF_020884755.1': 'Bifidobacterium mizhiense (strain=S053-2)',
                        'GCF_016101585.1': 'Bifidobacterium polysaccharolyticum (strain=W8117)',
                        'GCF_000499285.1': 'Bifidobacterium sp000499285 (strain=7101)'
}

recomputed_tree = ete3.Tree('additional_analyses/results/add_isolates/g__Bifidobacterium/Results_Nov21/Species_Tree/SpeciesTree_rooted.txt')
recomputed_tree.get_leaves()
# rename tip labels by adding the species name instead of the acession ID or MAG ID
for leaf in recomputed_tree:
    print(leaf.name)
    if leaf.name in bifido_isolate_names:
        leaf.name = bifido_isolate_names[leaf.name]
    else:
        species_name = phylo_metadata[phylo_metadata['ID'] == leaf.name].iloc[0]['MAG_species_name_final']
        if 'GCA' in leaf.name:
            leaf.name = f'{species_name}'
        else:
            leaf.name = f'{leaf.name}--{species_name}'
recomputed_tree.get_leaf_names()
recomputed_tree.write(outfile='additional_analyses/results/add_isolates/g__Bifidobacterium/SpeciesTree_renamed.txt', format=1)

# Gilliamella

gilli_isolate_genomes = {'GCF_000599985.1': 'Gilliamella apicola (strain=wkB1)',
                         'GCF_019469165.1': 'Gilliamella apicola_E (strain=ESL0443)',
                         'GCF_001693755.1': 'Gilliamella apicola_F (strain=wkB2)',
                         'GCF_001690185.1': 'Gilliamella apicola_I (strain=wkB308)',
                         'GCF_001690685.1': 'Gilliamella apicola_J (strain=wkB112)',
                         'GCF_001690705.1': 'Gilliamella apicola_K (strain=wkB178)',
                         'GCF_001690195.1': 'Gilliamella apicola_L (strain=wkB112)',
                         'GCF_001693435.1': 'Gilliamella apicola_N (strain=wkB7)',
                         'GCF_003202915.1': 'Gilliamella apicola_Q (strain=ESL0177)',
                         'GCF_002142155.1': 'Gilliamella apis (strain=NO3)',
                         'GCF_030758615.1': 'Gilliamella apis_A (strain=ESL0172)',
                         'GCF_002142215.1': 'Gilliamella sp002142215 (strain=N-G2)',
                         'GCF_019469185.1': 'Gilliamella sp019469185 (strain=ESL0441)',
                         'GCF_019469205.1': 'Gilliamella sp019469205 (strain=ESL0405)',
                         'GCF_028751545.1': 'Gilliamella sp945273075 (strain=B3022)',
                         'GCF_026536085.1': 'Gilliamella sp945276085 (strain=B3781)'
}

recomputed_tree = ete3.Tree('additional_analyses/results/add_isolates/g__Gilliamella/Results_Nov23/Species_Tree/SpeciesTree_rooted.txt')
recomputed_tree.get_leaves()
# rename tip labels by adding the species name instead of the acession ID or MAG ID
for leaf in recomputed_tree:
    print(leaf.name)
    if leaf.name in gilli_isolate_genomes:
        leaf.name = gilli_isolate_genomes[leaf.name]
    else:
        species_name = phylo_metadata[phylo_metadata['ID'] == leaf.name].iloc[0]['MAG_species_name_final']
        if 'GCA' in leaf.name:
            leaf.name = f'{species_name}'
        else:
            leaf.name = f'{leaf.name}--{species_name}'
recomputed_tree.get_leaf_names()
recomputed_tree.write(outfile='additional_analyses/results/add_isolates/g__Gilliamella/SpeciesTree_renamed.txt', format=1)

# Snodgrassella

snod_isolate_genomes = {'GCF_000600005.1': 'Snodgrassella alvi (strain=wkB2)',
                        'GCF_002777855.1': 'Snodgrassella alvi_E (strain=wkB298)',
                        'GCF_002777745.1': 'Snodgrassella alvi_D (strain=WF3-3)',
                        'GCF_026535915.1': 'Snodgrassella sp945267035 (strain=B3882)'
}

recomputed_tree = ete3.Tree('additional_analyses/results/add_isolates/g__Snodgrassella/Results_Nov21/Species_Tree/SpeciesTree_rooted.txt')
recomputed_tree.get_leaves()
# rename tip labels by adding the species name instead of the acession ID or MAG ID
for leaf in recomputed_tree:
    print(leaf.name)
    if leaf.name in snod_isolate_genomes:
        leaf.name = snod_isolate_genomes[leaf.name]
    else:
        species_name = phylo_metadata[phylo_metadata['ID'] == leaf.name].iloc[0]['MAG_species_name_final']
        if 'GCA' in leaf.name:
            leaf.name = f'{species_name}'
        else:
            leaf.name = f'{leaf.name}--{species_name}'
recomputed_tree.get_leaf_names()
recomputed_tree.write(outfile='additional_analyses/results/add_isolates/g__Snodgrassella/SpeciesTree_renamed.txt', format=1)