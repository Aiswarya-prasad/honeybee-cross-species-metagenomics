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
    for node in bact_tree.postorder():
    # for node in sym_tree.postorder():
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

with open(f'additional_analyses/results/12-cophylogeny_test_repeat/total_and_actual_nodes.csv', 'w') as f:
    f.write(f'clade,type,nodes\n')

for clade in rerooted:
    if exists(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}'):
        system(f'mv additional_analyses/results/12-cophylogeny_test_repeat/{clade} additional_analyses/results/12-cophylogeny_test_repeat/{clade}.old')
    makedirs(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}')
    makedirs(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/intermediate_files')
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
    # write incidence table to file
    interactions.to_csv(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/intermediate_files/interactions.tsv', sep='\t')

    set(host_tips) == set(interactions.columns)
    set(bact_tips) == set(interactions.index)
    
    for n in range(100):
        host_tips_permuted = np.random.permutation(host_tips)
        host_tree_permuted = host_tree.copy()
        for i, orig in enumerate(host_tree_permuted.tips()):
            orig.name = host_tips_permuted[i]
        # print(host_tree_permuted.ascii_art())
        nodes_permuted = hommola_traverse(host_tree_permuted,
                                bact_tree, 
                                interactions,
                                min_node_size=7,
                                max_node_depth=max_node_depth,
                                signodes_fp=f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % n,
                                results_fp=f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/intermediate_files/host_nodes.permuted.%s.pickle' % n)
        
    nodes = hommola_traverse(host_tree,
                         bact_tree, 
                         interactions,
                         min_node_size=7,
                         max_node_depth=max_node_depth,
                         signodes_fp=f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/host_signodes.txt',
                         results_fp=f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/host_nodes.pickle')
    
    min_node_size = 7
    nodes_tested = 0
    nodes_skipped = 0
    for node in bact_tree.postorder():
        nodes_tested += 1
        print(nodes_tested, end = '\r')
        sym_tips = [x.name for x in node.tips()]
        # host_tips = interact.loc[sym_tips, interact.loc[sym_tips, ].sum() > 0].columns
        # get number of hosts in the subtree
        host_tips = list(set([x.split('--')[1][0] for x in sym_tips]))
        node_depth = node.get_max_distance()[0]
        if len(sym_tips) < min_node_size or len(host_tips) < 3 or node_depth > max_node_depth:
            nodes_skipped += 1
            continue

    sizes = []
    for i in range(100):
        
        signodes_permuted = pd.read_csv(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % i, sep='\t')
        size = signodes_permuted.loc[(signodes_permuted['r'] > 0.75) &
                        (signodes_permuted['p'] < 0.01)].shape
        print(size)
        sizes.append(size[0])

    signodes_true = pd.read_csv(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/host_signodes.txt', sep='\t')
    size_main = signodes_true.loc[(signodes_true['r'] > 0.75) &
                        (signodes_true['p'] < 0.01)].shape

    permuted = sns.histplot(sizes)
    permuted.axvline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axvline(nodes_tested-nodes_skipped, color='black')

    # mark what the vertical lines mean
    plt.text(size_main[0], 0, f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(nodes_tested-nodes_skipped, 0, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')

    # add legend and x and y labels
    plt.xlabel('Number of Significant Nodes')
    plt.ylabel('Frequency')

    fig = permuted.get_figure()
    # add extra space on the right
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes.pdf', bbox_inches='tight')
    plt.close()
    # repeat but with kde plot
    permuted = sns.kdeplot(sizes)
    permuted.axvline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axvline(nodes_tested-nodes_skipped, color='black')

    # mark what the vertical lines mean
    plt.text(size_main[0], 0, f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(nodes_tested-nodes_skipped, 0, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')

    # add legend and x and y labels
    plt.xlabel('Number of Significant Nodes')
    plt.ylabel('Density')

    fig = permuted.get_figure()
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes_kde.pdf', bbox_inches='tight')
    plt.close()
    # also add a scatter plot with round in x axis and actual marked in red as the last point
    permuted = sns.scatterplot(x=range(100), y=sizes)
    # add a point for the actual
    permuted.scatter(x=100, y=size_main[0], color='darkred', label='actual')
    permuted.axhline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axhline(nodes_tested-nodes_skipped, color='black')

    # # mark what the vertical lines mean
    plt.text(100, size_main[0], f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(100, nodes_tested-nodes_skipped, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')

    # add legend and x and y labels
    plt.xlabel('Permutation round #')
    plt.ylabel('Number of Significant Nodes')

    fig = permuted.get_figure()
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes_scatter.pdf', bbox_inches='tight')

    # make a plot where the y axis is the number of significant nodes x axis is one point for the genus/clade and the size of the point 
    # is the frequency y value with a red line for the actual number of significant nodes
    # make a frequency table of the number of significant nodes
    df_plot_points = pd.DataFrame(sizes, columns=['size'])
    df_plot_points = df_plot_points['size'].value_counts().reset_index()
    df_plot_points.columns = ['size', 'frequency']
    df_plot_points = df_plot_points.sort_values(by='size')
    df_plot_points['clade'] = clade
    # write the info
    df_plot_points.to_csv(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/df_size_freq.csv')
    with open(f'additional_analyses/results/12-cophylogeny_test_repeat/total_and_actual_nodes.csv', 'a') as f:
        f.write(f'{clade},total,{nodes_tested-nodes_skipped}\n')
        f.write(f'{clade},actual,{size_main[0]}\n')
    # ensure the same size range across all loops by adding a fized size scale 
    # add size scale 100 =2, 50 = 1, 1 =0.1
    permuted = sns.scatterplot(x='clade', y='size', size='frequency', data=df_plot_points, color='#000000')
    permuted.axhline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axhline(nodes_tested-nodes_skipped, color='black')

    # mark what the vertical lines mean
    plt.text(clade, size_main[0], f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(clade, nodes_tested-nodes_skipped, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')
    
    fig = permuted.get_figure()

    # labels
    plt.xlabel('')
    plt.ylabel('Number of Significant Nodes')
    

    
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes_scatter_size.pdf', bbox_inches='tight')
    plt.close()
    
    # ensure the same size range across all loops by adding a fized size scale 
    # add size scale 100 =2, 50 = 1, 1 =0.1
    permuted = sns.scatterplot(x='clade', y='size', size='frequency', data=df_plot_points, color='#000000')
    permuted.axhline(size_main[0], color='darkred', linewidth=1, label='actual')

    # mark what the vertical lines mean
    plt.text(clade, size_main[0], f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    fig = permuted.get_figure()

    # labels
    plt.xlabel('')
    plt.ylabel('Number of Significant Nodes')
    

    
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes_scatter_size_wo_total.pdf', bbox_inches='tight')
    plt.close()

    plt.close('all')

# collect all the additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/intermediate_files/interactions.tsv into one folder
# with all the tsv files
makedirs('additional_analyses/results/12-cophylogeny_test_repeat/interactions')
for clade in rerooted:
    system(f'cp additional_analyses/results/12-cophylogeny_test_repeat/{clade}/intermediate_files/interactions.tsv additional_analyses/results/12-cophylogeny_test_repeat/interactions/{clade}_interactions.tsv')

# using these, The mapping is a text file that ends with the extension .mapping and specifies the association of the tips of the parasite tree to the tips of the host tree.  Each line in the file is of the form:
# parasiteTipName : hostTipName

# for each clade, write a mapping file
makedirs('additional_analyses/results/12-cophylogeny_test_repeat/mappings')
for clade in rerooted:
    bact_tree = rerooted[clade]
    with open(f'additional_analyses/results/12-cophylogeny_test_repeat/mappings/{clade}.mapping', 'w') as f:
        for tip in bact_tree.tips():
            tip.name = tip.name.replace(' ', '_')
            species = tip.name.split('--')[1][0]
            # print(f'{tip.name}:{species}')
            if species == 'M':
                species_name = 'Apis_mellifera'
            elif species == 'C':
                species_name = 'Apis_cerana'
            elif species == 'D':
                species_name = 'Apis_dorsata'
            elif species == 'F':
                species_name = 'Apis_florea'
            elif species == 'A':
                species_name = 'Apis_andreniformis'
            else:
                species_name = 'Outgroup'
                # continue
            f.write(f'{tip.name}:{species_name}\n')
        
        

# each is a dict of tree nodes now write them in nwk format into a file
rooted_mod = rerooted.copy()
makedirs(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees')
for clade in rooted_mod:
    # replace whitespace with underscores
    for tip in rooted_mod[clade].tips():
        tip.name = tip.name.replace(' ', '_')
    rooted_mod[clade].bifurcate()
    rooted_mod[clade].write(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/{clade}.nwk')


# only keep the species representative MAGs
makedirs(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees_mod')
rerooted_mod = rerooted.copy()
for clade in rerooted_mod:
    # delete mags that are not species representative
    tips_to_keep = []
    for tip in rerooted_mod[clade].tips():
        mag = tip.name.split('--')[1]
        if mag.startswith('G'):
            continue
        else:
            rep = mag_info[mag_info['ID'] == mag]['Representative'].values[0]
        if rep == 1:
            print(f'keeping {tip.name}')
            tips_to_keep.append(tip.name)
    rerooted_mod[clade]=rerooted_mod[clade].shear(tips_to_keep)
    # replace whitespace with underscores
    for tip in rerooted_mod[clade].tips():
        tip.name = tip.name.replace(' ', '_')
    rerooted_mod[clade].bifurcate()
    rerooted_mod[clade].write(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees_mod/{clade}.nwk')

magotu_prevs_df = pd.read_csv('results/figures/species_abundance_table.csv')
# only select MAG_species_name_final and prevalence_in_host
magotu_prevs_df = magotu_prevs_df[['MAG_species_name_final', 'prevalence_in_host', 'Host']]
# get unique rows
magotu_prevs_df = magotu_prevs_df.drop_duplicates()

makedirs(f'additional_analyses/results/12-cophylogeny_test_repeat/mappings_mod')
for clade in rerooted_mod:
    bact_tree = TreeNode.read(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees_mod/{clade}.nwk', convert_underscores=False)
    with open(f'additional_analyses/results/12-cophylogeny_test_repeat/mappings_mod/{clade}.mapping', 'w') as f:
        tips = list(bact_tree.tips())  # Static copy of tips to avoid traversal issues
        for tip in tips:
            print(tip.name)
            # tip.name = tip.name.replace(' ', '_')
            mag = tip.name.split('--')[1]
            species = phylo_metadata[phylo_metadata['ID'] == mag]['MAG_species_name_final'].values[0]
            prevs_sub_df = magotu_prevs_df[magotu_prevs_df['MAG_species_name_final'] == species]
            # replace NaN with 0
            prevs_sub_df['prevalence_in_host'] = prevs_sub_df['prevalence_in_host'].fillna(0)
            # if there is only one non-zero entry, take that
            if prevs_sub_df.shape[0] == 1:
                suffix = prevs_sub_df['Host'].values[0].split(' ')[1]
                f.write(f'{tip.name}:Apis_{suffix}\n')
            else:
                if prevs_sub_df.shape[0] == 0 or max(prevs_sub_df['prevalence_in_host']) == 0:
                    # remove from the tree
                    bact_tree.remove_deleted(lambda x: x.name == tip.name)
                else:
                    # Handle multiple non-zero entries
                    # Filter for prevalence > 0.1
                    filtered_df = prevs_sub_df[prevs_sub_df['prevalence_in_host'] > 0.1]

                    if filtered_df.empty:
                        bact_tree.remove_deleted(lambda x: x.name == tip.name)
                    elif filtered_df.shape[0] == 1:
                        suffix = filtered_df['Host'].values[0].split(' ')[1]
                        f.write(f'{tip.name}:Apis_{suffix}\n')
                    else:
                        # Convert tip to internal node
                        internal_node = TreeNode(name='')  # Create internal node
                        parent_node = tip.parent
                        
                        # Add the internal node as a child to the parent
                        parent_node.append(internal_node)
                        parent_node.remove(tip)  # Remove the original tip from the parent
                        
                        for _, row in filtered_df.iterrows():
                                suffix = row['Host'].split(' ')[1]
                                new_tip_name = f"{tip.name}_Apis_{suffix}"
                                new_tip = TreeNode(name=new_tip_name)
                                internal_node.append(new_tip)  # Add new tip as a child of the internal node
                                f.write(f'{new_tip.name}:Apis_{suffix}\n')
                        # remove the original tip
                        internal_node.remove(tip)
        print(bact_tree.ascii_art())
        bact_tree.bifurcate()
        bact_tree.prune()
        print(bact_tree.ascii_art())
        # # print all tips 
        # for tip in bact_tree.tips():
        #     print(tip.name)
        bact_tree.write(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees_mod/{clade}_mod.nwk')


def is_bifurcating(tree):
    """Check if all internal nodes of the tree are bifurcating."""
    for node in tree.traverse():
        if not node.is_tip() and len(node.children) != 2:
            return False
    return True


    # # remove all the quotes within the file
    # with open(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/{clade}.nwk', 'r') as f:
    #     lines = f.readlines()
    # with open(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/{clade}.nwk', 'w') as f:
    #     for line in lines:
    #         f.write(line.replace("'", ""))

    # with open(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/{clade}.nwk', 'r') as f:
    #     lines = f.readlines()
    # with open(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/{clade}.nwk', 'w') as f:
    #     for line in lines:
    #         f.write(line.replace(",", ", "))

    # with open(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/{clade}.nwk', 'r') as f:
    #     lines = f.readlines()
    # with open(f'additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/{clade}.nwk', 'w') as f:
    #     for line in lines:
    #         f.write(line.replace("s0", " s0"))

# make trww look lik:
# ((O_atricolis, A_subflava) h1, (((G_ianthinogaster, G_granatia) h64, ((E_dybowskii, (H_niveoguttatus, H_margaritatus) h71) h68, ((C_monteiri, ((L_larvata, (L_rara, ((L_r_cingica, L_r_rubricata) h100, (L_virata, (L_rhodopareia, L_sanguinodorsalis) h105) h101) h99) h97) h88, ((L_rufopicta, L_nitidula) h90, (L_senegala_rendalii, L_senegala_rhodopsis) h91) h89) h87) h74, ((P_melba_grotei, P_melba_citerior) h76, (P_afra, (P_lineata, (P_phoenicoptera, P_hypogrammica) h83) h81) h77) h75) h69) h65) h48, ((C_quartinia, C_melanotis) h50, ((E_astrild, (E_paludicola, (E_rhodophyga, (E_troglodytes, E_melpoda) h61) h59) h57) h54, E_erythronotos) h51) h49) h2) h0;


host_tree_carr = TreeNode.read('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree_Carr2023.nwk', convert_underscores=False)
host_tree_carr = host_tree_carr.shear(['Apis_mellifera', 'Apis_cerana', 'Apis_dorsata', 'Apis_florea', 'Apis_andreniformis', 'Outgroup'])
host_tree_carr.write('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree_Carr2023_sheared.nwk')

# replace whitespace with underscores

# add node labels to host tree
host_tree = TreeNode.read('results/figures/visualize_temp/host_tree.treefile', convert_underscores=False)
host_tree_mod = host_tree
for tip in host_tree_mod.tips():
    tip.name = tip.name.replace(' ', '_')
    if tip.name == 'A_mellifera':
        tip.name = 'Apis_mellifera'
    elif tip.name == 'A_cerana':
        tip.name = 'Apis_cerana'
    elif tip.name == 'A_dorsata':
        tip.name = 'Apis_dorsata'
    elif tip.name == 'A_florea':
        tip.name = 'Apis_florea'
    elif tip.name == 'A_andreniformis':
        tip.name = 'Apis_andreniformis'
    else:
        continue
# n = 1
# for node in host_tree_mod.traverse():
#     node.length = None
#     if node.is_tip():
#         node.name = node.name
#     # if node is non empty
#     else:
#         # write node number
#         node.name = f'h0{n}'
#         n += 1

host_tree_mod.write('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree.nwk')

# remove all the quotes within the file
with open('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree.nwk', 'r') as f:
    lines = f.readlines()
with open('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree.nwk', 'w') as f:
    for line in lines:
        f.write(line.replace("'", ""))


with open('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree.nwk', 'r') as f:
    lines = f.readlines()
with open('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree.nwk', 'w') as f:
    for line in lines:
        f.write(line.replace(",", ", "))

# with open('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree.nwk', 'r') as f:
#     lines = f.readlines()
# with open('additional_analyses/results/12-cophylogeny_test_repeat/nwk_trees/host_tree.nwk', 'w') as f:
#     for line in lines:
#         f.write(line.replace("h0", " h0"))



# collect kde and histogram plots into one folder
makedirs('additional_analyses/results/12-cophylogeny_test_repeat/histograms')
for clade in rerooted:
    system(f'cp additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes_kde.pdf additional_analyses/results/12-cophylogeny_test_repeat/histograms/{clade}_permuted_nodes_kde.pdf')
    system(f'cp additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes.pdf additional_analyses/results/12-cophylogeny_test_repeat/histograms/{clade}_permuted_nodes.pdf')

makedirs('additional_analyses/results/12-cophylogeny_test_repeat/scatter_size')
for clade in rerooted:
    system(f'cp additional_analyses/results/12-cophylogeny_test_repeat/{clade}/permuted_nodes_scatter_size.pdf additional_analyses/results/12-cophylogeny_test_repeat/scatter_size/{clade}_permuted_nodes_scatter_size.pdf')

# move all this to a new script within additional_analyses and keep this as was on github...

# for clade in rerooted:
#     if exists(f'results/figures/12-cophylogeny_test/{clade}'):
#         system(f'mv results/figures/12-cophylogeny_test/{clade} results/figures/12-cophylogeny_test/{clade}.old')
#     makedirs(f'results/figures/12-cophylogeny_test/{clade}')
#     makedirs(f'results/figures/12-cophylogeny_test/{clade}/intermediate_files')
#     host_tree = TreeNode.read('results/figures/visualize_temp/host_tree.treefile', convert_underscores=False)
#     for tip in host_tree.tips():
#         tip.name = tip.name.replace('A_', 'Apis ')
#     print(host_tree.ascii_art())

#     bact_tree = rerooted[clade]
#     bact_tree.get_max_distance()[0]
#     bact_tips = [x.name for x in bact_tree.tips()]
#     host_tips = [x.name for x in host_tree.tips()]

#     max_node_depth = bact_tree.get_max_distance()[0] #/4

#     incidence_dict = {}
#     species_list = {x:0 for x in ['M', 'C', 'D', 'F', 'A']}
#     for genome in [x.name for x in bact_tree.tips()]:
#         species = genome.split('--')[1][0]
#         genome_host = species_list.copy()
#         genome_host[species] = 1

#         incidence_dict[genome] = genome_host

#     incidence_table = pd.DataFrame.from_dict(incidence_dict,
#                                             orient='index')
#     interactions = incidence_table.copy()
#     interactions.columns = ['Apis mellifera', 'Apis cerana', 'Apis dorsata', 'Apis florea', 'Apis andreniformis', 'Other']
#     # drop the column Other
#     interactions = interactions.drop(columns=['Other'])

#     set(host_tips) == set(interactions.columns)
#     set(bact_tips) == set(interactions.index)
    
#     for n in range(100):
#         host_tips_permuted = np.random.permutation(host_tips)
#         host_tree_permuted = host_tree.copy()
#         for i, orig in enumerate(host_tree_permuted.tips()):
#             orig.name = host_tips_permuted[i]
#         nodes_permuted = hommola_traverse(host_tree_permuted,
#                                 bact_tree, 
#                                 interactions,
#                                 min_node_size=7,
#                                 max_node_depth=max_node_depth,
#                                 signodes_fp=f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % n,
#                                 results_fp=f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_nodes.permuted.%s.pickle' % n)
        
#     nodes = hommola_traverse(host_tree,
#                          bact_tree, 
#                          interactions,
#                          min_node_size=7,
#                          max_node_depth=max_node_depth,
#                          signodes_fp=f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt',
#                          results_fp=f'results/figures/12-cophylogeny_test/{clade}/host_nodes.pickle')
    
#     sizes = []
#     for i in range(100):
        
#         signodes_permuted = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % i, sep='\t')
#         size = signodes_permuted.loc[(signodes_permuted['r'] > 0.75) &
#                         (signodes_permuted['p'] < 0.01)].shape
#         print(size)
#         sizes.append(size[0])

#     permuted = sns.histplot(sizes)
#     permuted.axvline(206, color='darkred')

#     fig = permuted.get_figure()
#     fig.savefig(f'results/figures/12-cophylogeny_test/{clade}/permuted_nodes.pdf')



def print_tree_inline(x):
    try:
        return(str(ete3.Tree(x)))
    except:
        return(str(x))

clade = 'Lactobacillus'
clade = 'Frischella'
clade = 'Gilliamella'
clade = 'Snodgrassella'
clade = 'Bombilactobacillus'
clade = 'Bifidobacterium'
clade = 'Dysgonomonas'
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
        # drop 'root'
        species_names = [x.split('--')[0] for x in species_names if x != 'root']
        return(','.join(list(set(species_names))))
    except:
        species_names = ete3.Tree(x, quoted_node_names=True, format=1).get_leaf_names()
        # drop 'root'
        species_names = [x.split('--')[0] for x in species_names if x != 'root']
        return(','.join(list(set(species_names))))


df_res_all = pd.DataFrame()
files = [x for x in glob('results/figures/12-cophylogeny_test/*/host_signodes.txt') if 'old' not in x]
for file in files:
    df_read = pd.read_csv(file, sep = '\t')
    df_res_all = pd.concat([df_res_all, df_read], axis=0)
# only keep p < 0.01
df_res_all = df_res_all[df_res_all['p'] < 0.01]
df_res_all = df_res_all.sort_values(by='r', ascending=False)
df_res_all['species_list'] = df_res_all['sym_subtree'].apply(lambda x: get_species_from_tree(x))
# rearrange columns
df_res_all = df_res_all[['r','p','sym_tips_count','host_tips_count','host_tip_span','species_list','sym_subtree']]
# df_res_all['sym_subtree_drawn'] = df_res_all['sym_subtree'].apply(print_tree_inline)
# df_res_all['host_subtree_drawn'] = df_res_all['host_subtree'].apply(print_tree_inline)
# drop columns
# df_res_all  = df_res_all.drop(columns=['sym_subtree', 'host_subtree'])
df_res_all.to_csv(f'additional_analyses/results/12-cophylogeny_test_repeat/sig_nodes_all_with_trees.txt', sep = '\t', index = False)
# df_res_all.to_csv(f'results/figures/12-cophylogeny_test/sig_nodes_all.txt', sep = ',', index = False)

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

def get_species_from_tree_ext(x):
    try:
        species_names = ete3.Tree(x).get_leaf_names()
        # drop 'root'
        species_names = [x.split('--')[0]+'--'+x.split('--')[1][0] for x in species_names if x != 'root']
        return(','.join(list(set(species_names))))
    except:
        species_names = ete3.Tree(x, quoted_node_names=True, format=1).get_leaf_names()
        # drop 'root'
        species_names = [x.split('--')[0]+'--'+x.split('--')[1][0] for x in species_names if x != 'root']
        return(','.join(list(set(species_names))))

df_sig_nodes = pd.read_csv('additional_analyses/results/12-cophylogeny_test_repeat/sig_nodes_all_with_trees.txt', sep = '\t')
# read sym_subtree as tree and get tip label (split by -- and take first part)
df_sig_nodes['species_list'] = df_sig_nodes['sym_subtree'].apply(lambda x: get_species_from_tree(x))

df_specificity=pd.read_csv('additional_analyses/results/12-cophylogeny_test_repeat/species_specificity.csv')
df_specificity['species_name'] = df_specificity['Species'].apply(lambda x: x.replace(' ', '_'))
df_specificity['species_name_w_mag'] = df_specificity['species_name'] + '--' + df_specificity['Reference MAG']
specificity_dict = df_specificity.set_index('species_name').to_dict()['Host specificity']

def get_spec_list(x):
    spec_list = x.split(',')
    spec_list = [x.split('--')[0] for x in spec_list]
    spec_list_final = []
    for x in spec_list:
        if x.replace('\'','') not in specificity_dict:
            spec_list_final.append(f'unknown')
        else:
            spec_list_final.append(specificity_dict[x.replace('\'','')])
    return(','.join(spec_list_final))

df_sig_nodes['specificity'] = df_sig_nodes['species_list'].apply(lambda x: get_spec_list(x))

df_sig_nodes['species_list'] = df_sig_nodes['species_list'].apply(lambda x: x.replace('\'', ''))

df_sig_nodes = df_sig_nodes[['r', 'p', 'sym_tips_count', 'host_tips_count', 'host_tip_span', 'species_list', 'specificity', 'sym_subtree']]

# sort by species name
df_sig_nodes = df_sig_nodes.sort_values(by='species_list')

df_sig_nodes.to_csv('additional_analyses/results/12-cophylogeny_test_repeat/sig_nodes_all_with_trees_clean.txt', sep = '\t', index = False)

header = True
with open('additional_analyses/results/12-cophylogeny_test_repeat/sig_nodes_all_with_trees_clean.txt', 'r') as f:
    with open('additional_analyses/results/12-cophylogeny_test_repeat/sig_nodes_all_with_trees_ext_all.txt', 'w') as f_out:
        for line in f:
            if header:
                f_out.write(line)
                header = False
                continue
            line_split = line.strip().split('\t')
            with open('additional_analyses/results/12-cophylogeny_test_repeat/temp_tree.txt', 'w') as f_temp:
                f_temp.write(line_split[7])
            f_out.write('\t'.join(line_split[:7])+'\n')
            my_tree = TreeNode.read('additional_analyses/results/12-cophylogeny_test_repeat/temp_tree.txt', convert_underscores=False)
            f_out.write(my_tree.ascii_art())
            f_out.write('\n')

header = True
with open('additional_analyses/results/12-cophylogeny_test_repeat/sig_nodes_all_with_trees_clean.txt', 'r') as f:
    with open('additional_analyses/results/12-cophylogeny_test_repeat/sig_nodes_all_with_trees_ext.txt', 'w') as f_out:
        for line in f:
            if header:
                f_out.write(line)
                header = False
                continue
            line_split = line.strip().split('\t')
            if float(line_split[0]) < 0.75:
                continue
            with open('additional_analyses/results/12-cophylogeny_test_repeat/temp_tree.txt', 'w') as f_temp:
                f_temp.write(line_split[7])
            f_out.write('\t'.join(line_split[:7])+'\n')
            my_tree = TreeNode.read('additional_analyses/results/12-cophylogeny_test_repeat/temp_tree.txt', convert_underscores=False)
            f_out.write(my_tree.ascii_art())
            f_out.write('\n')

magotu_prevs_df = pd.read_csv('results/figures/species_abundance_table_fixed_manual.csv')
'''
During this analysis it became apprent that the following species names were mislabelled.
Probably did not have an impact in earlier analysis so the fix is done here
NA_NA-F5-1_7 to Bifidobacterium_polysaccharolyticum-F5-1_7
NA_NA-F4-5_12 to Gilliamella_apicola_K-F4-5_12
NA_NA-F4-5_6 to Apilactobacillus_zhangqiuensis-F4-5_6
NA_NA-F4-5_22 to Snodgrassella_SGB-58_1-F4-5_22
'''
# only select MAG_species_name_final and prevalence_in_host
magotu_prevs_df = magotu_prevs_df[['MAG_species_name_final', 'prevalence_in_host', 'Host']]
# get unique rows
magotu_prevs_df = magotu_prevs_df.drop_duplicates()
magotu_prevs_df['MAG_species_to_compare'] = magotu_prevs_df['MAG_species_name_final'].apply(lambda x: str(x).replace(' ', '_'))
magotu_prevs_df['Host_to_compare'] = magotu_prevs_df['Host'].apply(lambda x: x.split(' ')[1][0].upper())

for clade in rerooted:
    if exists(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}'):
        system(f'mv additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade} additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}.old')
    makedirs(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}')
    makedirs(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/intermediate_files')
    host_tree = TreeNode.read('results/figures/visualize_temp/host_tree.treefile', convert_underscores=False)
    for tip in host_tree.tips():
        tip.name = tip.name.replace('A_', 'Apis ')
    # print(host_tree.ascii_art())

    bact_tree = rerooted[clade]
    bact_tree.get_max_distance()[0]
    bact_tips = [x.name for x in bact_tree.tips()]
    host_tips = [x.name for x in host_tree.tips()]

    max_node_depth = bact_tree.get_max_distance()[0] #/4

    incidence_dict = {}
    species_list = {x:0 for x in ['M', 'C', 'D', 'F', 'A']}
    for genome in [x.name for x in bact_tree.tips()]:
        mag_species = genome.split('--')[0]
        if mag_species not in magotu_prevs_df['MAG_species_to_compare'].values:
            species = genome.split('--')[1][0]
            if species == 'G':
                genome_host = species_list
            else:
                print(f'{mag_species} not in species list')
                genome_host = species_list.copy()
                genome_host[species] = 1
        else:
            iter_df = magotu_prevs_df[magotu_prevs_df['MAG_species_to_compare'] == mag_species]
            # only keep prevalence in host > 0.1
            iter_df = iter_df[iter_df['prevalence_in_host'] > 0.1]
            genome_host = species_list.copy()
            for host_iter in iter_df['Host_to_compare'].values:
                genome_host[host_iter] = 1

        incidence_dict[genome] = genome_host

    incidence_table = pd.DataFrame.from_dict(incidence_dict,
                                            orient='index')
    interactions = incidence_table.copy()
    interactions.columns = ['Apis mellifera', 'Apis cerana', 'Apis dorsata', 'Apis florea', 'Apis andreniformis']
    # # drop the column Other
    # interactions = interactions.drop(columns=['Other'])
    # write incidence table to file
    interactions.to_csv(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/intermediate_files/interactions.tsv', sep='\t')

    set(host_tips) == set(interactions.columns)
    set(bact_tips) == set(interactions.index)
    
    for n in range(100):
        host_tips_permuted = np.random.permutation(host_tips)
        host_tree_permuted = host_tree.copy()
        for i, orig in enumerate(host_tree_permuted.tips()):
            orig.name = host_tips_permuted[i]
        # print(host_tree_permuted.ascii_art())
        nodes_permuted = hommola_traverse(host_tree_permuted,
                                bact_tree, 
                                interactions,
                                min_node_size=7,
                                max_node_depth=max_node_depth,
                                signodes_fp=f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % n,
                                results_fp=f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/intermediate_files/host_nodes.permuted.%s.pickle' % n)
        
    nodes = hommola_traverse(host_tree,
                         bact_tree, 
                         interactions,
                         min_node_size=7,
                         max_node_depth=max_node_depth,
                         signodes_fp=f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/host_signodes.txt',
                         results_fp=f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/host_nodes.pickle')
    
    min_node_size = 7
    nodes_tested = 0
    nodes_skipped = 0
    for node in bact_tree.postorder():
        nodes_tested += 1
        print(nodes_tested, end = '\r')
        sym_tips = [x.name for x in node.tips()]
        # host_tips = interact.loc[sym_tips, interact.loc[sym_tips, ].sum() > 0].columns
        # get number of hosts in the subtree
        host_tips = list(set([x.split('--')[1][0] for x in sym_tips]))
        node_depth = node.get_max_distance()[0]
        if len(sym_tips) < min_node_size or len(host_tips) < 3 or node_depth > max_node_depth:
            nodes_skipped += 1
            continue

    sizes = []
    for i in range(100):
        
        signodes_permuted = pd.read_csv(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/intermediate_files/host_signodes.permuted.%s.txt' % i, sep='\t')
        size = signodes_permuted.loc[(signodes_permuted['r'] > 0.75) &
                        (signodes_permuted['p'] < 0.01)].shape
        print(size)
        sizes.append(size[0])

    signodes_true = pd.read_csv(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/host_signodes.txt', sep='\t')
    size_main = signodes_true.loc[(signodes_true['r'] > 0.75) &
                        (signodes_true['p'] < 0.01)].shape

    permuted = sns.histplot(sizes)
    permuted.axvline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axvline(nodes_tested-nodes_skipped, color='black')

    # mark what the vertical lines mean
    plt.text(size_main[0], 0, f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(nodes_tested-nodes_skipped, 0, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')

    # add legend and x and y labels
    plt.xlabel('Number of Significant Nodes')
    plt.ylabel('Frequency')

    fig = permuted.get_figure()
    # add extra space on the right
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/permuted_nodes.pdf', bbox_inches='tight')
    plt.close()
    # repeat but with kde plot
    permuted = sns.kdeplot(sizes)
    permuted.axvline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axvline(nodes_tested-nodes_skipped, color='black')

    # mark what the vertical lines mean
    plt.text(size_main[0], 0, f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(nodes_tested-nodes_skipped, 0, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')

    # add legend and x and y labels
    plt.xlabel('Number of Significant Nodes')
    plt.ylabel('Density')

    fig = permuted.get_figure()
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/permuted_nodes_kde.pdf', bbox_inches='tight')
    plt.close()
    # also add a scatter plot with round in x axis and actual marked in red as the last point
    permuted = sns.scatterplot(x=range(100), y=sizes)
    # add a point for the actual
    permuted.scatter(x=100, y=size_main[0], color='darkred', label='actual')
    permuted.axhline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axhline(nodes_tested-nodes_skipped, color='black')

    # # mark what the vertical lines mean
    plt.text(100, size_main[0], f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(100, nodes_tested-nodes_skipped, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')

    # add legend and x and y labels
    plt.xlabel('Permutation round #')
    plt.ylabel('Number of Significant Nodes')

    fig = permuted.get_figure()
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/permuted_nodes_scatter.pdf', bbox_inches='tight')

    # make a plot where the y axis is the number of significant nodes x axis is one point for the genus/clade and the size of the point 
    # is the frequency y value with a red line for the actual number of significant nodes
    # make a frequency table of the number of significant nodes
    df_plot_points = pd.DataFrame(sizes, columns=['size'])
    df_plot_points = df_plot_points['size'].value_counts().reset_index()
    df_plot_points.columns = ['size', 'frequency']
    df_plot_points = df_plot_points.sort_values(by='size')
    df_plot_points['clade'] = clade
    # write the info
    df_plot_points.to_csv(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/df_size_freq.csv')
    with open(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/total_and_actual_nodes.csv', 'a') as f:
        f.write(f'{clade},total,{nodes_tested-nodes_skipped}\n')
        f.write(f'{clade},actual,{size_main[0]}\n')
    # ensure the same size range across all loops by adding a fized size scale 
    # add size scale 100 =2, 50 = 1, 1 =0.1
    permuted = sns.scatterplot(x='clade', y='size', size='frequency', data=df_plot_points, color='#000000')
    permuted.axhline(size_main[0], color='darkred', linewidth=1, label='actual')
    permuted.axhline(nodes_tested-nodes_skipped, color='black')

    # mark what the vertical lines mean
    plt.text(clade, size_main[0], f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    plt.text(clade, nodes_tested-nodes_skipped, f'#Nodes tested (after filtering)', rotation=0, verticalalignment='bottom')
    
    fig = permuted.get_figure()

    # labels
    plt.xlabel('')
    plt.ylabel('Number of Significant Nodes')
    

    
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/permuted_nodes_scatter_size.pdf', bbox_inches='tight')
    plt.close()
    
    # ensure the same size range across all loops by adding a fized size scale 
    # add size scale 100 =2, 50 = 1, 1 =0.1
    permuted = sns.scatterplot(x='clade', y='size', size='frequency', data=df_plot_points, color='#000000')
    permuted.axhline(size_main[0], color='darkred', linewidth=1, label='actual')

    # mark what the vertical lines mean
    plt.text(clade, size_main[0], f'#significant in actual comparison', rotation=0, verticalalignment='bottom')
    fig = permuted.get_figure()

    # labels
    plt.xlabel('')
    plt.ylabel('Number of Significant Nodes')
    

    
    fig.savefig(f'additional_analyses/results/12-cophylogeny_test_repeat_int_by_mapping/{clade}/permuted_nodes_scatter_size_wo_total.pdf', bbox_inches='tight')
    plt.close()

    plt.close('all')