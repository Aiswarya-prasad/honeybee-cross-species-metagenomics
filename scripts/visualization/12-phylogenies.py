import ete3

phylo_metadata = pd.read_csv('results/11_phylogenies/phylo_genomes_metadata.tsv', sep='\t')

# following ete3 tutorial to understand how to parse the trees
t = ete3.PhyloTree('results/11_phylogenies/03_iqtree_trees/g__Gilliamella/g__Gilliamella.treefile')
print(t)
