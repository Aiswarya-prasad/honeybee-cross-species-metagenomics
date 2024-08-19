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


'''
rename and prune the trees made from aa sequences at the genus level
and write them to a new directory
'''
# "F4-5", "F5-1", "M6-2", "D9-5", "F7-5"    
samples_to_remove = ['F4-5', 'F5-1', 'M6-2', 'D9-5', 'F7-5']
# following ete3 tutorial to understand how to parse the trees
for genus in [basename(x).split('.')[0] for x in glob('results/11_phylogenies/03_iqtree_trees/*/*.treefile')]:
    t = ete3.PhyloTree(f'results/11_phylogenies/03_iqtree_trees/{genus}/{genus}.treefile')

    for i, node in enumerate(t.traverse('postorder')):
        if not node.is_leaf():
            continue
        print(f'{i}/{len(t)}', end = '\r')
        # for each node, if it matches an ID in phylo_metadata, rename it to the final species name prefixed to MAG name
        if node.name in list(phylo_metadata['ID']):
            node.name = phylo_metadata[phylo_metadata['ID'] == node.name]['MAG_species_name_final_nospace'].values[0] + '--' + node.name
        # if not, append quality info from mag_info to it
        else:
            if node.name in list(mag_info['ID']):
                node.name = node.name + '--' + mag_info[mag_info['ID'] == node.name]['Quality'].values[0]
            else:
                if node.name != '':
                    node.name = get_tax_info(node.name)['species']
        # replace spaces with underscores
        node.name = node.name.replace(' ', '_')
    tips_to_keep = [x.name for x in t.iter_leaves()]
    len(tips_to_keep)
    for tip in tips_to_keep:
        for sample in samples_to_remove:
            if sample in tip:
                print(f'removing {tip}')
                tips_to_keep.remove(tip)
    len(tips_to_keep)
    t.prune(tips_to_keep)

    t.write(outfile=f'results/figures/visualize_temp/renamed_trees/{genus}.treefile', format=1)

'''
code from Sanders et al 2022 (modified!)
https://github.com/CUMoellerLab/Sanders-etal-2022-analysis/blob/main/notebooks/prepare_GWCodeML.ipynb
'''

def count_clade_hosts(tree, focal_host):
    '''
    '''
    tips = [x.name for x in tree.tips()]
    spp = [x.split('--')[1][0] for x in tips]
    focal_tips = sum([x == focal_host for x in spp])
    if focal_tips == 0:
        return(None)
    tree.assign_ids()
    clade_counts = {}
    for clade in tree.postorder():
        if clade.is_tip():
            clade_tips = [clade.name]
        else:
            clade_tips = [x.name for x in clade.tips()]
        clade_spp = [x.split('--')[1][0] for x in clade_tips]
        clade_focal_tips = sum([x == focal_host for x in clade_spp])
        clade_counts[clade.id] = {'id': clade.id,
                                  'tips': len(clade_tips),
                                  'focal_spec': focal_host,
                                  'focal_tips': clade_focal_tips,
                                  'prop_of_foci': clade_focal_tips / focal_tips,
                                  'focal_prop': clade_focal_tips / len(clade_tips)}
    clade_df = pd.DataFrame.from_dict(clade_counts, orient='index')
    return(clade_df)
    
trees_dir = 'results/figures/visualize_temp/renamed_trees/*.treefile'
tree_fps = glob(trees_dir)

trees = {}
for tree_fp in tree_fps:
    if 'g__' not in tree_fp:
        continue
    tree = TreeNode.read(tree_fp, 
                         convert_underscores=False)
    cluster = basename(tree_fp).split('.')[0].split('__')[1]
    trees[cluster] = tree

# trees = {}
# dfs = []
# focal_species_list = ['M', 'C', 'D', 'F', 'A']

# for tree_fp in tree_fps:
#     tree = TreeNode.read(tree_fp, 
#                          convert_underscores=False)
#     cluster = basename(tree_fp).split('.')[0].split('__')[1]
#     for focus in focal_species_list:
#         cluster_df = count_clade_hosts(tree, focus)
#         if cluster_df is not None:
#             trees[cluster] = tree
#             cluster_df['cluster'] = cluster
#             dfs.append(cluster_df)

# 
def is_mag(name):
    if '--' in name:
        return True
    else:
        return False


def reroot(target, source):
    target_outgroups = []
    for child in target.children:
        if child.is_tip():
            target_outgroups.append(child.name)
    for child in source.children:
        if child.is_tip():
            if child.name in target_outgroups:
                return(target)
            else:
                return(target.root_at(target.find(child.name).parent))
        tips = [x.name for x in child.tips()]
        target_node = target.lowest_common_ancestor(tips)
        if len([x for x in target_node.tips()]) == len(tips):
            return(target.root_at(target_node))
    raise ValueError

rerooted = {}

# list of bee species related genomes not to be used to root
bee_species_related = ['Apis mellifera', 'Apis cerana', 'Apis dorsata', 'honeybee',
                    'Apis florea', 'Apis andreniformis', 'Bee gut', 'bee', 'adult honey bee queen gut',
                    'Honey', 'Apis mellifera gut', 'Apis mellifera gut', 'apis mellifera']
other_genomes_df = phylo_metadata[~phylo_metadata.Host.isin(bee_species_related)]
other_genomes = other_genomes_df['MAG_species_name_final_nospace'].apply(lambda x: x.replace(' ', '_')) + '--' + other_genomes_df['ID']

for clade in trees:
    tree = trees[clade]
    root_genomes = [x.name for x in tree.tips() if x.name in list(other_genomes.values)]
    if len(root_genomes) == 0:
        print(f'No genomes from other environments in clade {clade}\n')
        continue
    super_clade = tree.lowest_common_ancestor(root_genomes)
    # super_clade = super_tree.lowest_common_ancestor([x.name for x in super_tree.tips() if clade in x.name])
    # where which node to use as root? get LCA of the genomes from 'other environments'
    try:
        rerooted[clade] = reroot(tree, super_clade)
    except:
        print('Failure on clade %s\n' % clade)
        print('Target tree:\n')
        print(tree.ascii_art())
        print('Source tree:\n')
        print(super_clade.ascii_art())


'''
cophylogeny tests using Hommola test as in Sanders et al 2022
'''

def write_node(f,
               r,
               p,
               node,
               sym_tips,
               host_tips,
               host_tree):
    '''
    Writes calculated output for a node to an open filehandle.
    
    f: open filehandle
    r: correlation coefficient from Hommola test
    p: p value from Hommola test
    node: skbio.TreeNode object for tested symbiont node
    sym_tips: list of symbiont tip names in node
    host_tips: host tips represented by symbionts from tested node
    host_tree: complete host tree
    '''
    nodeio = StringIO()
    node.write(nodeio)
    nodetree_str = nodeio.getvalue().strip()
    host_subtree = host_tree.lowest_common_ancestor(host_tips)
    hostio = StringIO()
    host_subtree.write(hostio)
    hosttree_str = hostio.getvalue().strip()
    host_span = len([x.name for x in host_subtree.tips()])
    host_depth = host_subtree.get_max_distance()[0]
    sym_depth = node.get_max_distance()[0]
    outline = '\t'.join([str(r),
                         str(p),
                         str(len(sym_tips)),
                         str(len(host_tips)),
                         str(host_span),
                         nodetree_str,
                         hosttree_str,
                         str(host_depth),
                         str(sym_depth)])
    f.write(outline + '\n')
    f.flush()


def hommola_traverse(host_tree,
                     sym_tree,
                     interact,
                     min_node_size=7,
                     max_node_depth=1.0,
                     perms=100,
                     results_fp=None,
                     signodes_fp=None,
                     sigval=0.05):
    '''
    Recursive test using Hommola test at each node of the symbiont phylogeny.
    
    Returns dictionary of results, one item per node tested.
    
    host_tree: the skbio.TreeNode object for the hosts
    sym_tree: the skbio.TreeNode object for the symbionts
    interact: pd.DataFrame object with symbiont tip names as index and host tip names as column names
    min_node_size: minimum number of symbiont tips in a node in order to execute test
    max_node_depth: maximum depth of node in order to execute test (max tip to tip distance)
    perms: permutations to run for Hommola test
    results_fp: fp at which to write results as Pickle file
    signodes_fp: fp at which to write calculated values for significant nodes
    sigval: significance threshold for writing signodes
    
    '''
    node_dict = {}
    if signodes_fp:
        signodes_f = open(signodes_fp, 'w')
        header = '\t'.join(['r',
                            'p',
                            'sym_tips_count',
                            'host_tips_count',
                            'host_tip_span',
                            'sym_subtree',
                            'host_subtree',
                            'host_depth',
                            'sym_depth'])
        signodes_f.write('{0}\n'.format(header))
    
    host_dists = host_tree.tip_tip_distances()
    nodes_tested = 0
    nodes_skipped = 0
    for node in sym_tree.postorder():
        nodes_tested += 1
        print(nodes_tested, end = '\r')
        sym_tips = [x.name for x in node.tips()]
        host_tips = interact.loc[sym_tips, interact.loc[sym_tips, ].sum() > 0].columns
        node_depth = node.get_max_distance()[0]
        if len(sym_tips) < min_node_size or len(host_tips) < 3 or node_depth > max_node_depth:
            nodes_skipped += 1
            continue
        sym_dists = node.tip_tip_distances()
        subhost_dists =  host_dists.filter(host_tips)
        subint = interactions.loc[sym_tips, host_tips]
        r, p, _ = hommola_cospeciation(subhost_dists, sym_dists, subint)
        node_dict[node.id] = [r, p, node]
        if signodes_fp and p <= sigval:
#             try:
            write_node(signodes_f,
                       r,
                       p,
                       node,
                       sym_tips,
                       host_tips,
                       host_tree)
#             except:
#                 break
    print(f'{nodes_tested} nodes tested, {nodes_skipped} nodes skipped')
    if signodes_fp:
        signodes_f.close()
    if results_fp:
        pickle.dump(node_dict,
                    open(results_fp, 'wb'))
    return(node_dict)

for clade in rerooted:
    if exists(f'results/figures/12-cophylogeny_test/{clade}'):
        system(f'mv results/figures/12-cophylogeny_test/{clade} results/figures/12-cophylogeny_test/{clade}.old')
    makedirs(f'results/figures/12-cophylogeny_test/{clade}')
    makedirs(f'results/figures/12-cophylogeny_test/{clade}/intermediate_files')
    host_tree = TreeNode.read('results/figures/visualize_temp/host_tree.treefile', convert_underscores=False)
    for tip in host_tree.tips():
        tip.name = tip.name.replace('A_', 'Apis ')
    print(host_tree.ascii_art())

    bact_tree = rerooted[clade]
    bact_tree.get_max_distance()[0]
    bact_tips = [x.name for x in bact_tree.tips()]
    host_tips = [x.name for x in host_tree.tips()]

    max_node_depth = bact_tree.get_max_distance()[0] #/4

    incidence_dict = {}
    species_list = {x:0 for x in ['M', 'C', 'D', 'F', 'A']}
    for genome in [x.name for x in bact_tree.tips()]:
        species = genome.split('--')[1][0]
        genome_host = species_list.copy()
        genome_host[species] = 1

        incidence_dict[genome] = genome_host

    incidence_table = pd.DataFrame.from_dict(incidence_dict,
                                            orient='index')
    interactions = incidence_table.copy()
    interactions.columns = ['Apis mellifera', 'Apis cerana', 'Apis dorsata', 'Apis florea', 'Apis andreniformis', 'Other']
    # drop the column Other
    interactions = interactions.drop(columns=['Other'])

    set(host_tips) == set(interactions.columns)
    set(bact_tips) == set(interactions.index)
    
    for n in range(100):
        host_tips_permuted = np.random.permutation(host_tips)
        host_tree_permuted = host_tree.copy()
        for i, orig in enumerate(host_tree_permuted.tips()):
            orig.name = host_tips_permuted[i]
        nodes_permuted = hommola_traverse(host_tree_permuted,
                                bact_tree, 
                                interactions,
                                min_node_size=7,
                                max_node_depth=max_node_depth,
                                signodes_fp=f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % n,
                                results_fp=f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_nodes.permuted.%s.pickle' % n)
        
    nodes = hommola_traverse(host_tree,
                         bact_tree, 
                         interactions,
                         min_node_size=7,
                         max_node_depth=max_node_depth,
                         signodes_fp=f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt',
                         results_fp=f'results/figures/12-cophylogeny_test/{clade}/host_nodes.pickle')
    
    sizes = []
    for i in range(100):
        
        signodes_permuted = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % i, sep='\t')
        size = signodes_permuted.loc[(signodes_permuted['r'] > 0.75) &
                        (signodes_permuted['p'] < 0.01)].shape
        print(size)
        sizes.append(size[0])

    permuted = sns.histplot(sizes)
    permuted.axvline(206, color='darkred')

    fig = permuted.get_figure()
    fig.savefig(f'results/figures/12-cophylogeny_test/{clade}/permuted_nodes.pdf')

def print_tree_inline(x):
    try:
        return(str(ete3.Tree(x)))
    except:
        return(str(x))

clade = 'Lactobacillus'
clade = 'Frischella'
clade = 'Gilliamella'
clade = 'Snodgrassella'
clade = 'Dysgonomonas'
clade = 'Bombilactobacillus'
clade = 'Bifidobacterium'
df_nodes_tested_info = pd.DataFrame()
host_tree = TreeNode.read('results/figures/visualize_temp/host_tree.treefile', convert_underscores=False)
for tip in host_tree.tips():
    tip.name = tip.name.replace('A_', 'Apis ')
# print(host_tree.ascii_art())
for clade in rerooted:
    bact_tree = rerooted[clade]
    bact_tree.get_max_distance()[0]
    bact_tips = [x.name for x in bact_tree.tips()]
    max_node_depth = bact_tree.get_max_distance()[0] #/4
    incidence_dict = {}
    species_list = {x:0 for x in ['M', 'C', 'D', 'F', 'A']}
    for genome in [x.name for x in bact_tree.tips()]:
        species = genome.split('--')[1][0]
        genome_host = species_list.copy()
        genome_host[species] = 1

        incidence_dict[genome] = genome_host

    incidence_table = pd.DataFrame.from_dict(incidence_dict,
                                            orient='index')
    interactions = incidence_table.copy()
    interactions.columns = ['Apis mellifera', 'Apis cerana', 'Apis dorsata', 'Apis florea', 'Apis andreniformis', 'Other']
    # drop the column Other
    interactions = interactions.drop(columns=['Other'])
    host_dists = host_tree.tip_tip_distances()
    min_node_size = 7
    interact = interactions
    node_dict = {}
    host_dists = host_tree.tip_tip_distances()
    max_node_depth = bact_tree.get_max_distance()[0]
    nodes_skipped = 0
    nodes_present = 0
    host_specific_node = 0
    sig_nodes = 0
    for node in bact_tree.postorder():
        nodes_present += 1
        # print(nodes_present, end = '\r')
        sym_tips = [x.name for x in node.tips()]
        host_tips = interact.loc[sym_tips, interact.loc[sym_tips, ].sum() > 0].columns
        node_depth = node.get_max_distance()[0]
        if len(sym_tips) < min_node_size or len(host_tips) < 3 or node_depth > max_node_depth:
            nodes_skipped += 1
            if len(host_tips) == 1:
                host_specific_node += 1
            continue
        sym_dists = node.tip_tip_distances()
        subhost_dists =  host_dists.filter(host_tips)
        subint = interactions.loc[sym_tips, host_tips]
        r, p, _ = hommola_cospeciation(subhost_dists, sym_dists, subint)
        if p < 0.01 and r > 0.75:
            sig_nodes += 1
        node_dict[node.id] = [r, p, node]
    # add info abour number of nodes of diff kinds to dataframe
    df_nodes_tested_info = df_nodes_tested_info._append({'clade': clade,
                                                        'nodes_present': nodes_present,
                                                        'nodes_skipped': nodes_skipped,
                                                        'nodes_tested': nodes_present - nodes_skipped,
                                                        'host_specific_node': host_specific_node,
                                                        'sig_nodes': sig_nodes}, ignore_index=True)
    print(f'{clade}: {nodes_present} nodes tested, {nodes_skipped} nodes skipped ({nodes_skipped/nodes_present*100}), {host_specific_node} host specific nodes ({host_specific_node/nodes_present*100}), {sig_nodes} significant nodes ({sig_nodes/nodes_present*100})')
df_nodes_tested_info.to_csv('results/figures/12-cophylogeny_test/nodes_tested_info.tsv', sep='\t', index=False)
'''
Snodgrassella: 274 nodes tested, 260 nodes skipped (94.89), 100 host specific nodes (36.49), 2 significant nodes (0.72)
Bifidobacterium: 492 nodes tested, 458 nodes skipped (93.08), 185 host specific nodes (37.60), 3 significant nodes (0.60)
Bartonella_A: 182 nodes tested, 157 nodes skipped (86.26), 44 host specific nodes (24.17), 1 significant nodes (0.54)
Frischella: 152 nodes tested, 144 nodes skipped (94.73), 65 host specific nodes (42.76), 1 significant nodes (0.65)
Lactobacillus: 680 nodes tested, 641 nodes skipped (94.26), 246 host specific nodes (36.17), 1 significant nodes (0.14)
Commensalibacter: 112 nodes tested, 109 nodes skipped (97.32), 43 host specific nodes (38.39), 3 significant nodes (2.67)
Pectinatus: 72 nodes tested, 68 nodes skipped (94.44), 17 host specific nodes (23.61), 0 significant nodes (0.0)
Gilliamella: 484 nodes tested, 453 nodes skipped (93.59), 158 host specific nodes (32.64), 0 significant nodes (0.0)
Bombilactobacillus: 504 nodes tested, 471 nodes skipped (93.45), 162 host specific nodes (32.14), 2 significant nodes (0.39)
Apilactobacillus: 60 nodes tested, 60 nodes skipped (100.0), 12 host specific nodes (20.0), 0 significant nodes (0.0)
Apibacter: 144 nodes tested, 144 nodes skipped (100.0), 68 host specific nodes (47.22), 0 significant nodes (0.0)
Dysgonomonas: 260 nodes tested, 255 nodes skipped (98.07), 93 host specific nodes (35.76), 1 significant nodes (0.38)
'''

clade = 'Lactobacillus'
clade = 'Frischella'
clade = 'Gilliamella'
clade = 'Snodgrassella'
clade = 'Dysgonomonas'
clade = 'Bombilactobacillus'
clade = 'Bifidobacterium'
sizes = []
for i in range(100):
        
    signodes_permuted = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % i, sep='\t')
    size = signodes_permuted.loc[(signodes_permuted['r'] > 0.75) &
                    (signodes_permuted['p'] < 0.01)].shape
    # print(size)
    sizes.append(size[0])
signodes = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep='\t')
size_main = signodes.loc[(signodes['r'] > 0.75) &
                        (signodes['p'] < 0.01)].shape
# sizes
np.median(sizes)
np.std(sizes)
size[0]

# the paper says r > 0.75, non-parametric P < 0.01 is co-diversifying

clade = 'Lactobacillus'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
df_res['sym_subtree_drawn'] = df_res['sym_subtree'].apply(print_tree_inline)
df_res['host_subtree_drawn'] = df_res['host_subtree'].apply(print_tree_inline)
tree = ete3.Tree(df_res['sym_subtree'][0])
# write dataframe to file
df_res.to_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes_drawn.txt', sep = '\t', index = False)
type(print(tree))
tree = ete3.Tree(df_res['sym_subtree'][12])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][14])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][7])
print(tree)


clade = 'Frischella'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
df_res['sym_subtree_drawn'] = df_res['sym_subtree'].apply(print_tree_inline)
df_res['host_subtree_drawn'] = df_res['host_subtree'].apply(print_tree_inline)
df_res.to_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes_drawn.txt', sep = '\t', index = False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][3])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][1])
print(tree)


clade = 'Gilliamella'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
df_res['sym_subtree_drawn'] = df_res['sym_subtree'].apply(print_tree_inline)
df_res['host_subtree_drawn'] = df_res['host_subtree'].apply(print_tree_inline)
df_res.to_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes_drawn.txt', sep = '\t', index = False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][4])
print(tree)


clade = 'Snodgrassella'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
df_res['sym_subtree_drawn'] = df_res['sym_subtree'].apply(print_tree_inline)
df_res['host_subtree_drawn'] = df_res['host_subtree'].apply(print_tree_inline)
df_res.to_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes_drawn.txt', sep = '\t', index = False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)


clade = 'Dysgonomonas'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
df_res['sym_subtree_drawn'] = df_res['sym_subtree'].apply(print_tree_inline)
df_res['host_subtree_drawn'] = df_res['host_subtree'].apply(print_tree_inline)
df_res.to_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes_drawn.txt', sep = '\t', index = False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)


clade = 'Bombilactobacillus'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
df_res['sym_subtree_drawn'] = df_res['sym_subtree'].apply(print_tree_inline)
df_res['host_subtree_drawn'] = df_res['host_subtree'].apply(print_tree_inline)
df_res.to_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes_drawn.txt', sep = '\t', index = False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][17])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][18])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][20])
print(tree)


clade = 'Bifidobacterium'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
df_res['sym_subtree_drawn'] = df_res['sym_subtree'].apply(print_tree_inline)
df_res['host_subtree_drawn'] = df_res['host_subtree'].apply(print_tree_inline)
df_res.to_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes_drawn.txt', sep = '\t', index = False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][1])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][19])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][20])
print(tree)

def get_species_from_tree(x):
    try:
        species_names = ete3.Tree(x).get_leaf_names()
        species_names = [x.split('--')[0] for x in species_names]
        return(' '.join(list(set(species_names))))
    except:
        return('')
    
df_res_all = pd.DataFrame()
files = [x for x in glob('results/figures/12-cophylogeny_test/*/host_signodes.txt') if 'old' not in x]
for file in files:
    df_read = pd.read_csv(file, sep = '\t')
    df_res_all = pd.concat([df_res_all, df_read], axis=0)
# only keep p < 0.01
df_res_all = df_res_all[df_res_all['p'] < 0.01]
df_res_all = df_res_all.sort_values(by='r', ascending=False)
df_res_all['species_list'] = df_res_all['sym_subtree'].apply(lambda x: get_species_from_tree(x))
# df_res_all['sym_subtree_drawn'] = df_res_all['sym_subtree'].apply(print_tree_inline)
# df_res_all['host_subtree_drawn'] = df_res_all['host_subtree'].apply(print_tree_inline)
# drop columns
df_res_all  = df_res_all.drop(columns=['sym_subtree', 'host_subtree'])
df_res_all.to_csv(f'results/figures/12-cophylogeny_test/sig_nodes_all.txt', sep = ',', index = False)

clades = ['Apibacter', 'Apilactobacillus', 'Bartonella_A', 'Commensalibacter', 'Lactobacillus', 'Frischella', 'Gilliamella', 'Snodgrassella', 'Dysgonomonas', 'Bombilactobacillus', 'Bifidobacterium']
df_out = pd.DataFrame()
for clade in clades:
    sig_nodes_perm = []
    sig_nodes_perm_05 = []
    sig_nodes_perm_075 = []
    bact_tree = rerooted[clade]
    max_node_depth = bact_tree.get_max_distance()[0] #/4
    min_node_size = 7
    nodes_tested = 0
    host_tree = TreeNode.read('results/figures/visualize_temp/host_tree.treefile', convert_underscores=False)
    bact_tips = [x.name for x in bact_tree.tips()]
    host_tips = [x.name for x in host_tree.tips()]
    for node in bact_tree.postorder():
            node_depth = node.get_max_distance()[0]
            if len(bact_tips) < min_node_size or len(host_tips) < 3 or node_depth > max_node_depth:
                continue
            nodes_tested += 1
    for i in range(100):
        df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_signodes.permuted.{i}.txt', sep = '\t')
        sig_nodes_perm.append(df_res.loc[(df_res['r'] > 0.0) & (df_res['p'] < 0.05)].shape[0])
        sig_nodes_perm_05.append(df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.01)].shape[0])
        sig_nodes_perm_075.append(df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.05)].shape[0])
    df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
    sig_nodes = df_res.loc[(df_res['r'] > 0) & (df_res['p'] < 0.05)].shape[0]
    median_sig_nodes_perm = np.median(sig_nodes_perm)
    std_sig_nodes_perm = np.std(sig_nodes_perm)
    sig_nodes_05 = df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.01)].shape[0]
    median_sig_nodes_perm_05 = np.median(sig_nodes_perm_05)
    std_sig_nodes_perm_05 = np.std(sig_nodes_perm_05)
    sig_nodes_075 = df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.05)].shape[0]
    median_sig_nodes_perm_075 = np.median(sig_nodes_perm_075)
    std_sig_nodes_perm_075 = np.std(sig_nodes_perm_075)
    df_out = df_out._append({'clade': clade,
                            'total_tested_nodes': nodes_tested,
                            'sig_rel': sig_nodes,
                            'sig_perm_median_re,': median_sig_nodes_perm,
                            'sig_perm_sd_rel': std_sig_nodes_perm,
                            'sig_med': sig_nodes_05,
                            'sig_perm_median_med': median_sig_nodes_perm_05,
                            'sig_perm_sd_med': std_sig_nodes_perm_05,
                            'sig_strict': sig_nodes_075,
                            'sig_perm_median_strict': median_sig_nodes_perm_075,
                            'sig_perm_sd_strict': std_sig_nodes_perm_075}, ignore_index=True)
df_out.to_csv('results/figures/12-cophylogeny_test/sig_nodes_summary.txt', sep = '\t', index = False)

for clade in clades:
    # get the number of total nodes and the number of significant nodes
    # in df res and average across each of the intermediate files
    sig_nodes_perm = []
    bact_tree = rerooted[clade]
    max_node_depth = bact_tree.get_max_distance()[0] #/4
    min_node_size = 7
    nodes_tested = 0
    host_tree = TreeNode.read('results/figures/visualize_temp/host_tree.treefile', convert_underscores=False)
    bact_tips = [x.name for x in bact_tree.tips()]
    host_tips = [x.name for x in host_tree.tips()]
    for node in bact_tree.postorder():
            node_depth = node.get_max_distance()[0]
            if len(bact_tips) < min_node_size or len(host_tips) < 3 or node_depth > max_node_depth:
                continue
            nodes_tested += 1
    for i in range(100):
        df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_signodes.permuted.{i}.txt', sep = '\t')
        sig_nodes_perm.append(df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.05)].shape[0])
        # sig_nodes_perm.append(df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.01)].shape[0])
    df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
    sig_nodes = df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.05)].shape[0]
    # sig_nodes = df_res.loc[(df_res['r'] > 0.75) & (df_res['p'] < 0.05)].shape[0]
    median_sig_nodes_perm = np.median(sig_nodes_perm)
    print(f'{clade}: {sig_nodes}/{nodes_tested} and in permutations: {median_sig_nodes_perm}/{nodes_tested}')

'''
using p < 0.05 and r > 0
Apibacter: 0/144 and in permutations: 0/144
Apilactobacillus: 0/60 and in permutations: 0/60
Bartonella_A: 1/182 and in permutations: 1/182
Commensalibacter: 3/112 and in permutations: 3/112
Lactobacillus: 2/680 and in permutations: 0/680
Frischella: 2/152 and in permutations: 1/152
Gilliamella: 0/484 and in permutations: 0/484
Snodgrassella: 2/274 and in permutations: 2/274
Dysgonomonas: 1/260 and in permutations: 1/260
Bombilactobacillus: 2/504 and in permutations: 1/504
Bifidobacterium: 3/492 and in permutations: 4/492

using p < 0.05 and r > 0
Apibacter: 0/144 and in permutations: 0/144
Apilactobacillus: 0/60 and in permutations: 0/60
Bartonella_A: 16/182 and in permutations: 9/182
Commensalibacter: 3/112 and in permutations: 3/112
Lactobacillus: 24/680 and in permutations: 22/680
Frischella: 7/152 and in permutations: 6/152
Gilliamella: 12/484 and in permutations: 11/484
Snodgrassella: 5/274 and in permutations: 5/274
Dysgonomonas: 3/260 and in permutations: 3/260
Bombilactobacillus: 29/504 and in permutations: 24/504
Bifidobacterium: 33/492 and in permutations: 33/492
'''

'''
identify species (clades) to consider for estimate substitution rates
per time for as the one that are monophyletic in the bact tree
for at least 3 hosts based on the species tree made using the
bac120 marker genes amino acid sequences

and then calculate the substitution rate per time for each of those
using the codon-aware alignment of nucleotide sequences for these 
bac120 marker genes
'''

# trees = {}
# dfs = []
# focal_species_list = ['M', 'C', 'D', 'F', 'A']

# for tree_fp in tree_fps:
#     tree = TreeNode.read(tree_fp, 
#                          convert_underscores=False)
#     cluster = basename(tree_fp).split('.')[0].split('__')[1]
#     for focus in focal_species_list:
#         cluster_df = count_clade_hosts(tree, focus)
#         if cluster_df is not None:
#             trees[cluster] = tree
#             cluster_df['cluster'] = cluster
#             dfs.append(cluster_df)

'''
Rates of sequence evolution were estimated using a linear regression of the genetic distance calculated between each pair 
of MAGs and the evolutionary divergence time of their respective hosts
'''

def host_from_tip(tipname):
    letter = tipname.split('--')[1][0]
    if letter == 'M':
        return 'Apis mellifera'
    elif letter == 'C':
        return 'Apis cerana'
    elif letter == 'D':
        return 'Apis dorsata'
    elif letter == 'F':
        return 'Apis florea'
    elif letter == 'A':
        return 'Apis andreniformis'
    else:
        return 'Other'

def species_from_tip(tipename):
    return tipename.split('--')[0]


# # for each genus, (only makes sense if found in > 3 hosts)
# # calculate the genetic distance between each pair of MAGs
# # calculate the evolutionary divergence time of their respective hosts
# # iterate through each pair of MAGs hat are tips of the bact tree using combinations
# # and calculate the slope of the linear regression
def host_divergence(a, b, type = 'dist'):
    # Carr, S. M. Multiple mitogenomes indicate Things Fall Apart with Out of Africa or Asia hypotheses for the phylogeographic evolution of Honey Bees (Apis mellifera). Sci Rep 13, 9386 (2023).
    host_tree = TreeNode.read(
        StringIO('(Bombus_ignitus:1938, ((A_florea:553, A_andreniformis:593):421, ((A_dorsata:608, A_laboriosa:579):346, (A_mellifera:651, (A_koschevnikovi:634, (A_nuluensis:419, (A_cerana:338, A_nigrocincta:379):161):381):262):335):312):107);'),
        convert_underscores=False
        )
    for tip in host_tree.tips():
        tip.name = tip.name.replace('A_', 'Apis ')
    host_dists = host_tree.tip_tip_distances().to_data_frame()
    if type == 'dist':
        return host_dists.loc[a, b]
    if type == 'time':
        # from literature https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-023-35937-4/MediaObjects/41598_2023_35937_MOESM6_ESM.jpg
        # lengths measured using https://eleif.net/photomeasure#howto for nodes that are not marked
        # 0.127/0.0115 = 11.04 Myr for flo - mel
        # check
        if set([a, b]) == set(['Apis mellifera', 'Apis cerana']):
            return 7.15
        if set([a, b]) == set(['Apis mellifera', 'Apis dorsata']):
            return 9.92
        if set([a, b]) == set(['Apis mellifera', 'Apis florea']):
            return 11.04
        if set([a, b]) == set(['Apis mellifera', 'Apis andreniformis']):
            return 11.04
        if set([a, b]) == set(['Apis dorsata', 'Apis cerana']):
            return 9.92
        if set([a, b]) == set(['Apis dorsata', 'Apis florea']):
            return 11.04
        if set([a, b]) == set(['Apis dorsata', 'Apis andreniformis']):
            return 11.04
        if set([a, b]) == set(['Apis cerana', 'Apis florea']):
            return 11.04
        if set([a, b]) == set(['Apis cerana', 'Apis andreniformis']):
            return 11.04
        if set([a, b]) == set(['Apis florea', 'Apis andreniformis']):
            return 6.42

'''
rename and prune the tree made from dna sequences of all MAGs
and write to a new file
'''

samples_to_remove = ['F4-5', 'F5-1', 'M6-2', 'D9-5', 'F7-5']
t = ete3.PhyloTree(f'results/11_phylogenies/05_MAG_bac120_nucleotide_trees/MAGs_bac120_nuc/MAGs_bac120_nuc.treefile')
for i, node in enumerate(t.traverse('postorder')):
        if not node.is_leaf():
            continue
        print(f'{i}', end = '\r')
        # for each node, if it matches an ID in phylo_metadata, rename it to the final species name prefixed to MAG name
        if node.name in list(phylo_metadata['ID']):
            node.name = phylo_metadata[phylo_metadata['ID'] == node.name]['MAG_species_name_final_nospace'].values[0] + '--' + node.name
        # if not, append quality info from mag_info to it
        else:
            if node.name in list(mag_info['ID']):
                node.name = node.name + '--' + mag_info[mag_info['ID'] == node.name]['Quality'].values[0]
            else:
                if node.name != '':
                    node.name = get_tax_info(node.name)['species']
        # replace spaces with underscores
        node.name = node.name.replace(' ', '_')
        tips_to_keep = [x.name for x in t.iter_leaves()]
        for tip in tips_to_keep:
            for sample in samples_to_remove:
                if sample in tip:
                    print(f'removing {tip}')
                    tips_to_keep.remove(tip)
        t.prune(tips_to_keep)
print(t)
t.write(outfile=f'results/figures/visualize_temp/renamed_trees/MAGs_bac120_nuc.treefile', format=1)

def get_mag_and_id(gene):
    dict_return = {}
    dict_return["id"] = gene.split("_")[-1]
    dict_return["mag"] = "_".join(gene.split("_")[:-1])
    return dict_return



'''
calculate the genetic distance between each pair of MAGs
and the evolutionary divergence time of their respective hosts
in results/figures/12-phylo_distances/bact_host_dists.tsv
'''

t_filt = TreeNode.read('results/figures/visualize_temp/renamed_trees/MAGs_bac120_nuc.treefile', convert_underscores=False)
t = TreeNode.read('results/11_phylogenies/05_MAG_bac120_nucleotide_trees/MAGs_bac120_nuc/MAGs_bac120_nuc.treefile', convert_underscores=False)
mags_in_tree = set()
for tip in t.tips():
    print(tip.name)
    mags_in_tree.add(tip.name)
    continue
    spec = get_tax_info(tip.name)
    if spec is None:
        continue
    print(spec['species'])
tip_tip_distances = t_filt.tip_tip_distances().to_data_frame()
header = f'MAG1\tMAG2\tHost1\tHost2\tSpecies1\tSpecies2\tHost_dist\tBact_dist\tHost_time\n'
with open(f'results/figures/12-phylo_distances/bact_host_dists.tsv', 'w+') as f:
    f.write(header)
for pair in combinations(t_filt.tips(), 2):
    if host_from_tip(pair[0].name) == 'Other' or host_from_tip(pair[1].name) == 'Other':
        continue
    else:
        host_dist = host_divergence(host_from_tip(pair[0].name), host_from_tip(pair[1].name), 'dist')
        # https://www.nature.com/articles/s41598-023-35937-4#Sec2
        # host_time = host_dist/11043 * 0.0115
        host_time = host_divergence(host_from_tip(pair[0].name), host_from_tip(pair[1].name), 'time')
        bact_dist = tip_tip_distances.loc[pair[0].name, pair[1].name]
        with open(f'results/figures/12-phylo_distances/bact_host_dists.tsv', 'a') as f:
            success = f.write(f'{pair[0].name}\t{pair[1].name}\t{host_from_tip(pair[0].name)}\t{host_from_tip(pair[1].name)}\t{species_from_tip(pair[0].name)}\t{species_from_tip(pair[1].name)}\t{host_dist}\t{bact_dist}\t{host_time}\n')

tip_dist_out = pd.read_csv('results/figures/12-phylo_distances/bact_host_dists.tsv', sep = '\t')
tip_dist_out['Genus1'] = tip_dist_out['Species1'].apply(lambda x: x.split('_')[0])
tip_dist_out['Genus2'] = tip_dist_out['Species2'].apply(lambda x: x.split('_')[0])
tip_dist_out['diff_genus'] = tip_dist_out['Genus1'] != tip_dist_out['Genus2']
np.median(tip_dist_out['Bact_dist'])
np.median(tip_dist_out[tip_dist_out['diff_genus'] == True]['Bact_dist'])
# for each genus calculate the minimum and maximum distance
# and the median and standard deviation
df_info_out = pd.DataFrame()
for genus in tip_dist_out['Genus1'].unique():
    df_genus = tip_dist_out[(tip_dist_out['Genus1'] == genus) & (tip_dist_out['Genus2'] == genus)]
    if genus == 's':
        print(f'{df_genus}')
    print(f'{genus} median: {np.median(df_genus["Bact_dist"])}, std: {np.std(df_genus["Bact_dist"])}, min: {np.min(df_genus["Bact_dist"])}, max: {np.max(df_genus["Bact_dist"])}')
    df_info_out = df_info_out._append({'Genus': genus,
                                        'median': np.median(df_genus["Bact_dist"]),
                                        'std': np.std(df_genus["Bact_dist"]),
                                        'min': np.min(df_genus["Bact_dist"]),
                                        'max': np.max(df_genus["Bact_dist"])}, ignore_index=True)
df_info_out.to_csv('results/figures/12-phylo_distances/bact_dist_genus_info.tsv', sep = '\t', index = False)

# make a set of all pairs of species and collect the list of bac-dist for each
species_pairs = {}
for row in tip_dist_out.iterrows():
    species_pair = tuple(sorted([row[1]['Species1'], row[1]['Species2']]))
    if row[1]['Genus1'] != row[1]['Genus2']:
        continue
    if species_pair in species_pairs:
        species_pairs[species_pair].append(row[1]['Bact_dist'])
    else:
        species_pairs[species_pair] = [row[1]['Bact_dist']]
df_spec_pairs_out = pd.DataFrame()
for pair in species_pairs:
    print(f'{pair}: {np.median(species_pairs[pair])}, {np.std(species_pairs[pair])}, {np.min(species_pairs[pair])}, {np.max(species_pairs[pair])}')
    df_spec_pairs_out = df_spec_pairs_out._append({'Species1': pair[0],
                                                    'Species2': pair[1],
                                                    'same_genus': pair[0].split('_')[0] == pair[1].split('_')[0],
                                                    'same_species': pair[0] == pair[1],
                                                    'median': np.median(species_pairs[pair]),
                                                    'std': np.std(species_pairs[pair]),
                                                    'min': np.min(species_pairs[pair]),
                                                    'max': np.max(species_pairs[pair])}, ignore_index=True)
df_spec_pairs_out.to_csv('results/figures/12-phylo_distances/bact_dist_species_pairs_info.tsv', sep = '\t', index = False)

species_pairs_all = {}
for row in tip_dist_out.iterrows():
    species_pair = tuple(sorted([row[1]['Species1'], row[1]['Species2']]))
    # if row[1]['Genus1'] != row[1]['Genus2']:
    #     continue
    if species_pair in species_pairs_all:
        species_pairs_all[species_pair].append(row[1]['Bact_dist'])
    else:
        species_pairs_all[species_pair] = [row[1]['Bact_dist']]
df_spec_pairs_out = pd.DataFrame()
for pair in species_pairs_all:
    print(f'{pair}: {np.median(species_pairs_all[pair])}, {np.std(species_pairs_all[pair])}, {np.min(species_pairs_all[pair])}, {np.max(species_pairs_all[pair])}')
    df_spec_pairs_out = df_spec_pairs_out._append({'Species1': pair[0],
                                                    'Species2': pair[1],
                                                    'same_genus': pair[0].split('_')[0] == pair[1].split('_')[0],
                                                    'same_species': pair[0] == pair[1],
                                                    'median': np.median(species_pairs_all[pair]),
                                                    'std': np.std(species_pairs_all[pair]),
                                                    'min': np.min(species_pairs_all[pair]),
                                                    'max': np.max(species_pairs_all[pair])}, ignore_index=True)
df_spec_pairs_out.to_csv('results/figures/12-phylo_distances/bact_dist_species_pairs_all_info.tsv', sep = '\t', index = False)

# box plot of same spec, diff spec, same genus and diff genus

fig, ax = plt.subplots()
sns.boxplot(x = 'same_genus', y = 'median', data = df_spec_pairs_out, ax = ax)
fig.savefig('results/figures/12-phylo_distances/bact_dist_same_diff_genus.png')
fig, ax = plt.subplots()
sns.boxplot(x = 'same_species', y = 'median', data = df_spec_pairs_out, ax = ax)
fig.savefig('results/figures/12-phylo_distances/bact_dist_same_diff_species.png')



np.std(tip_dist_out[tip_dist_out['diff_genus'] == True]['Bact_dist'])
np.min(tip_dist_out[tip_dist_out['diff_genus'] == True]['Bact_dist'])
np.max(tip_dist_out[tip_dist_out['diff_genus'] == True]['Bact_dist'])
# # show the row that has the min distance (print the whole thing)
# tip_dist_out[tip_dist_out['Bact_dist'] == np.min(tip_dist_out[tip_dist_out['diff_genus'] == True]['Bact_dist'])]
np.median(tip_dist_out[tip_dist_out['diff_genus'] == False]['Bact_dist'])
np.std(tip_dist_out[tip_dist_out['diff_genus'] == False]['Bact_dist'])
np.min(tip_dist_out[tip_dist_out['diff_genus'] == False]['Bact_dist'])
np.max(tip_dist_out[tip_dist_out['diff_genus'] == False]['Bact_dist'])
# make a boxplot of bact dist for same and diff genus
fig, ax = plt.subplots()
sns.boxplot(x = 'diff_genus', y = 'Bact_dist', data = tip_dist_out, ax = ax)
fig.savefig('results/figures/12-phylo_distances/bact_dist_same_diff_genus.png')
# for clade in rerooted:
#     makedirs(f'results/figures/12-phylo_distances/{clade}', exist_ok=True)
#     print(clade)
#     bact_tree = rerooted[clade]
#     bact_tips = [x.name for x in bact_tree.tips()]
#     bact_dists = bact_tree.tip_tip_distances().to_data_frame()
#     # get pairwise the distance between each pair of MAGs
#     # from bact_dist
#     # and the distance between their respective hosts
#     # from host_dists
#     # also write a tsv with mag pairs in columns 1 and 2 and host
#     # dist and bact dist in columns 3 and 4
#     # if one of the hosts is Other, put the host dist down as
#     # twice the max host dist
# # 0.0115 substitutions/site/Myr (Brower 1994 in Ref.31)
# # Papadopoulou, A., Anastasiou, I. & Vogler, A. P. Revisiting the insect mitochondrial molecular clock: The Mid-Aegean trench calibration. Mol. Biol. Evol. 27, 1659–1672 (2010).    
# # confirm that the length of mitogenome 
#     header = f'MAG1\tMAG2\tHost1\tHost2\tSpecies1\tSpecies2\tHost_dist\tBact_dist\tHost_time\n'
#     with open(f'results/figures/12-phylo_distances/{clade}/bact_host_dists.tsv', 'w+') as f:
#         f.write(header)
#     for pair in combinations(bact_tips, 2):
#         if host_from_tip(pair[0]) == 'Other' or host_from_tip(pair[1]) == 'Other':
#             continue
#         else:
#             host_dist = host_divergence(host_from_tip(pair[0]), host_from_tip(pair[1]), 'dist')
#         # https://www.nature.com/articles/s41598-023-35937-4#Sec2
#         # host_time = host_dist/11043 * 0.0115
#         host_time = host_divergence(host_from_tip(pair[0]), host_from_tip(pair[1]), 'time')
#         bact_dist = bact_dists.loc[pair[0], pair[1]]
#         with open(f'results/figures/12-phylo_distances/{clade}/bact_host_dists.tsv', 'a') as f:
#             success = f.write(f'{pair[0]}\t{pair[1]}\t{host_from_tip(pair[0])}\t{host_from_tip(pair[1])}\t{species_from_tip(pair[0])}\t{species_from_tip(pair[1])}\t{host_dist}\t{bact_dist}\t{host_time}\n')