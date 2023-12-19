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
    
samples_to_remove = ['F4-5', 'F5-1', 'M6-2']
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
dfs = []
focal_species_list = ['M', 'C', 'D', 'F', 'A']

for tree_fp in tree_fps:
    tree = TreeNode.read(tree_fp, 
                         convert_underscores=False)
    cluster = basename(tree_fp).split('.')[0].split('__')[1]
    for focus in focal_species_list:
        cluster_df = count_clade_hosts(tree, focus)
        if cluster_df is not None:
            trees[cluster] = tree
            cluster_df['cluster'] = cluster
            dfs.append(cluster_df)

# 
def is_mag(name):
    if '--' in name:
        return True
    else:
        return False

# t = ete3.PhyloTree('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-bac120/MAGs_refs.treefile')
# for node in t.traverse('postorder'):
#     if '--low' in node.name:
#         node.delete()
#     else:
#         if not any([is_mag(x) for x in node.get_leaf_names()]):
#             pass
#         else:
#             node.delete()

# for i, node in enumerate(t.traverse('postorder')):
#     if not node.is_leaf():
#         continue
#     print(f'{i}/{len(t)}', end = '\r')
#         # for each node, if it matches an ID in phylo_metadata, rename it to the final species name prefixed to MAG name
#     if node.name in list(phylo_metadata['ID']):
#         node.name = phylo_metadata[phylo_metadata['ID'] == node.name]['MAG_species_name_final_nospace'].values[0] + '--' + node.name
#     # if not, append quality info from mag_info to it
#     else:
#         if node.name in list(mag_info['ID']):
#             node.name = node.name + '--' + mag_info[mag_info['ID'] == node.name]['Quality'].values[0]
#         else:
#             if node.name != '':
#                 node.name = get_tax_info(node.name)['species']
#     # replace spaces with underscores
#     node.name = node.name.replace(' ', '_')

# t.write(outfile='results/figures/visualize_temp/MAGs_refs_renamed_nospace.treefile', format=1)
# print(t)

# super_tree_fp = 'results/figures/visualize_temp/MAGs_refs_renamed_nospace.treefile'

# super_tree = TreeNode.read(super_tree_fp,
#                            convert_underscores=False)

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

# host_dfs = {}
# target_hosts = ['M', 'C', 'D', 'F', 'A']

# for clade in rerooted:
#     tree = rerooted[clade]
#     tips = [x.name for x in tree.tips()]
#     host = []
#     for tip in tips:
#         host.append(tip.split('--')[1][0])
#     clade_df = pd.DataFrame.from_dict({'Tip': tips,
#                                        'Host': host})
#     clade_df['Clade'] = clade
#     clade_df['Branch'] = 0
#     host_dfs[clade] = clade_df

#     for i, target in enumerate(target_hosts):
#         target_tips = list(clade_df.loc[clade_df['Host'] == target, 'Tip'])
#         if len(target_tips) < 1:
#             continue
#         if len(target_tips) == 1:
#             clade_df.loc[clade_df['Tip'] == target_tips[0], 'Branch'] = i + 1
#             continue
#         target_node = tree.lowest_common_ancestor(target_tips)
#         lca_tips = [x.name for x in target_node.tips()]
#         if len(lca_tips) == len(target_tips):
#             clade_df.loc[clade_df['Tip'].isin(target_tips), 'Branch'] = i + 1
    

# manually transcribed tree from 
# # Carr, S. M. Multiple mitogenomes indicate Things Fall Apart with Out of Africa or Asia hypotheses for the phylogeographic evolution of Honey Bees (Apis mellifera). Sci Rep 13, 9386 (2023).
# host_tree = ete3.Tree(
#     '''(Bombus_ignitus:1938, ((A_florea:553, A_andreniformis:593):421, ((A_dorsata:608, A_laboriosa:579):346, (A_mellifera:651, (A_koschevnikovi:634, (A_nuluensis:419, (A_cerana:338, A_nigrocincta:379):161):381):262):335):312):107);'''
# )

# manually deleted other species on itol and not the tree is,
host_tree = ete3.Tree(
    '''(Bombus_ignitus:1938,((A_florea:553,A_andreniformis:593):421,((A_dorsata:608):346,(A_mellifera:651,(((A_cerana:338):161):381):262):335):312):107);'''
)
# write host tree into file
host_tree.prune(['A_mellifera', 'A_cerana', 'A_dorsata', 'A_florea', 'A_andreniformis'])
print(host_tree)
host_tree.write(outfile='results/figures/visualize_temp/host_tree.treefile', format=1)


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
    for node in sym_tree.postorder():
        nodes_tested += 1
        print(nodes_tested)
        sym_tips = [x.name for x in node.tips()]
        host_tips = interact.loc[sym_tips, interact.loc[sym_tips, ].sum() > 0].columns
        node_depth = node.get_max_distance()[0]
        if len(sym_tips) < min_node_size or len(host_tips) < 3 or node_depth > max_node_depth:
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

clade = 'Lactobacillus'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][12])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][14])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][7])
print(tree)

clade = 'Frischella'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
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
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][4])
print(tree)


clade = 'Snodgrassella'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)


clade = 'Dysgonomonas'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)


clade = 'Bombilactobacillus'
df_res = pd.read_csv(f'results/figures/12-cophylogeny_test/{clade}/host_signodes.txt', sep = '\t')
df_res = df_res.sort_values(by='r', ascending=False)
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
print(df_res)
tree = ete3.Tree(df_res['sym_subtree'][0])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][1])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][19])
print(tree)
tree = ete3.Tree(df_res['sym_subtree'][20])
print(tree)


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
        StringIO('(Bombus_ignitus:1938, ((A_florea:553, A_andreniformis:593):421, ((A_dorsata:608, A_laboriosa:579):346, (A_mellifera:651, (A_koschevnikovi:634, (A_nuluensis:419, (A_cerana:338, A_nigrocincta:379):161):381):262):335):312):107);')
        )
    for tip in host_tree.tips():
        tip.name = tip.name.replace('A_', 'Apis ')
        tip.name = tip.name.replace('A ', 'Apis ')
    host_dists = host_tree.tip_tip_distances().to_data_frame()
    if type == 'dist':
        return host_dists.loc[a, b]
    if type == 'time':
        return (host_dists.loc[a, b]/11043)/0.0115
    # from literature https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-023-35937-4/MediaObjects/41598_2023_35937_MOESM6_ESM.jpg
    # 0.0115 substitutions/site/Myr (Brower 1994 in Ref.31)
    # 11043 sites in mitogenome
    


for clade in rerooted:
    makedirs(f'results/figures/12-phylo_distances/{clade}', exist_ok=True)
    print(clade)
    bact_tree = rerooted[clade]
    bact_tips = [x.name for x in bact_tree.tips()]
    bact_dists = bact_tree.tip_tip_distances().to_data_frame()
    # get pairwise the distance between each pair of MAGs
    # from bact_dist
    # and the distance between their respective hosts
    # from host_dists
    # also write a tsv with mag pairs in columns 1 and 2 and host
    # dist and bact dist in columns 3 and 4
    # if one of the hosts is Other, put the host dist down as
    # twice the max host dist
# 0.0115 substitutions/site/Myr (Brower 1994 in Ref.31)
# Papadopoulou, A., Anastasiou, I. & Vogler, A. P. Revisiting the insect mitochondrial molecular clock: The Mid-Aegean trench calibration. Mol. Biol. Evol. 27, 1659–1672 (2010).    
# confirm that the length of mitogenome 
    header = f'MAG1\tMAG2\tHost1\tHost2\tSpecies1\tSpecies2\tHost_dist\tBact_dist\tHost_time\n'
    with open(f'results/figures/12-phylo_distances/{clade}/bact_host_dists.tsv', 'w+') as f:
        f.write(header)
    for pair in combinations(bact_tips, 2):
        if host_from_tip(pair[0]) == 'Other' or host_from_tip(pair[1]) == 'Other':
            continue
        else:
            host_dist = host_divergence(host_from_tip(pair[0]), host_from_tip(pair[1]), 'dist')
        # https://www.nature.com/articles/s41598-023-35937-4#Sec2
        # host_time = host_dist/11043 * 0.0115
        host_time = host_divergence(host_from_tip(pair[0]), host_from_tip(pair[1]), 'time')
        bact_dist = bact_dists.loc[pair[0], pair[1]]
        with open(f'results/figures/12-phylo_distances/{clade}/bact_host_dists.tsv', 'a') as f:
            success = f.write(f'{pair[0]}\t{pair[1]}\t{host_from_tip(pair[0])}\t{host_from_tip(pair[1])}\t{species_from_tip(pair[0])}\t{species_from_tip(pair[1])}\t{host_dist}\t{bact_dist}\t{host_time}\n')
