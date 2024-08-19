# 3. Look at alpha diversity (number of units per sample, cumulative curve)
    # the output of these can then be parsed and plotted in R

host_species = ['Apis mellifera', 'Apis cerana',
                'Apis dorsata', 'Apis florea', 'Apis andreniformis']

def host_of_sample(sample_name):
    '''
    function to get host name to add to summary files below
    '''
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

# 3a. genes

genes_per_sample = {}
for sample in samples:
    genes_per_sample[sample] = len(genes_detected[sample])
with open('results/figures/08-summarize_functions/genes_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tgenes\n')
    for sample in genes_per_sample:
        success = out_fh.write(f'{sample}\t{genes_per_sample[sample]}\n')

# 3b. Clusters
'''
clusters_per_sample = {sample: set() for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}', end = '\r')
    for cluster in clusters:
        for gene in clusters[cluster].genes:
            if gene in genes_detected[sample]:
                clusters_per_sample[sample].add(cluster)
                break
pickle.dump(clusters_per_sample, open('results/figures/08-summarize_functions/clusters_per_sample.pkl', 'wb'))
'''
clusters_per_sample = pd.read_pickle('results/figures/08-summarize_functions/clusters_per_sample.pkl')

samples_to_exclude = ["F4-5", "F5-1", "M6-2", "D9-5", "F7-5"]

with open('results/figures/08-summarize_functions/clusters_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tclusters\n')
    for sample in clusters_per_sample:
        if sample in samples_to_exclude:
            continue
        success = out_fh.write(f'{sample}\t{len(clusters_per_sample[sample])}\n')

with open('results/figures/08-summarize_functions/cluster_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tgenes\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(10):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/11 and sample size {j}/{len(samples_host)}', end = '\r')
                cum_clusters = 0
                samples_iter = r.sample(samples_host, sample_size)
                clusters_detected_iter = set()
                for sample in samples_iter:
                    clusters_detected_iter.update(clusters_per_sample[sample])
                cum_clusters = len(clusters_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_clusters}\n')

# 3c. KOs

# count the number of kos and under each category for each sample
kos_detected = {sample: {} for sample in samples}
ko_count = {sample: {'Total': 0, 'Annotated': 0} for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}')
    for gene in genes_detected[sample]:
        # kos_detected[sample]['NA'] = set()
        kos_detected[sample]['NA'] = 0
        ko_gene = kegg_dict[gene]
        ko_count[sample]['Total'] += 1
        if ko_gene == '':
            # kos_detected[sample]['NA'].add(gene)
            kos_detected[sample]['NA'] += 1
        else:
            ko_count[sample]['Annotated'] += 1
            if ko_gene in kos_detected[sample].keys():
                kos_detected[sample][ko_gene] += 1
            else:
                kos_detected[sample][ko_gene] = 1
                # kos_detected[sample][ko_gene] = set([gene])

# write an output to summarize Total and annotated genes from ko_count
with open('results/figures/08-summarize_functions/kos_per_sample_anno.tsv', 'w') as out_fh:
    out_fh.write('sample\tkos\tunannotated\tperc_anno\n')
    for sample in genes_per_sample:
        ann = str(ko_count[sample]['Annotated'])
        unann = str(ko_count[sample]['Total'] - ko_count[sample]['Annotated'])
        if ko_count[sample]['Total'] == 0:
            perc_anno = 'NA'
        else:
            perc_anno = str(round(ko_count[sample]['Annotated']/ko_count[sample]['Total']*100, 2))
        out_fh.write(f'{sample}\t{ann}\t{unann}\t{perc_anno}\n')

# median perc anno
median_kos = []
median_genes = []
for sample in samples:
    kos_in_sample = sum(kos_detected[sample].values())
    genes_in_sample = len(genes_detected[sample])
    median_kos.append(kos_in_sample)
    median_genes.append(genes_in_sample)
np.median(median_kos)
np.median(median_genes)

# total kos
total_kos = 0
for sample in samples:
    total_kos += sum(kos_detected[sample].values())

# summarize number of kos found per sample
with open('results/figures/08-summarize_functions/kos_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tkos\n')
    for sample in kos_detected:
        out_fh.write(f'{sample}\t{len(kos_detected[sample].keys())}\n')

# cum curve of #ko and #category
num_iters = 20
with open('results/figures/08-summarize_functions/KO_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tkos\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(num_iters):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/{num_iters} and sample size {j}/{len(samples_host)}', end = '\r')
                cum_kos = 0
                samples_iter = r.sample(samples_host, sample_size)
                kos_detected_iter = set()
                for sample in samples_iter:
                    kos_detected_iter.update([x for x in kos_detected[sample].keys()])
                cum_kos = len(kos_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_kos}\n')

##############################################################################################
# prepare the dictionary to get kegg ortholog information
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

with open('results/figures/08-summarize_functions/kegg_info_dict.txt', 'w+') as f:
    header = 'ko\tA\tB\tC\n'
    success = f.write(header)
    for ko in kegg_info_dict.keys():
        line = f'{ko}\t{kegg_info_dict[ko]["A"]}\t{kegg_info_dict[ko]["B"]}\t{kegg_info_dict[ko]["C"]}\n'
        success = f.write(line)


# 3d. Cazymes
                
cazymes_detected = {sample: {} for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}')
    for gene in genes_detected[sample]:
        if gene not in cayman_info_dict.keys():
            continue
        cazyme = cayman_info_dict[gene]
        if cazyme in cazymes_detected[sample].keys():
            cazymes_detected[sample][cazyme] += 1
        else:
            cazymes_detected[sample][cazyme] = 1

total_cazyme_genes = 0
for sample in samples:
    total_cazyme_genes += len(cazymes_detected[sample])

cazymes_per_host_median = {host: [] for host in host_species}
for sample in cazymes_detected.keys():
    host = host_of_sample(sample)
    cazymes_per_host_median[host].append(len(cazymes_detected[sample].keys()))

for host in cazymes_per_host_median.keys():
    print(f'{host} median cazymes: {np.median(cazymes_per_host_median[host])}')

GH_per_host_median = {host: [] for host in host_species}
for sample in cazymes_detected.keys():
    host = host_of_sample(sample)
    GH_per_host_median[host].append(len([x for x in cazymes_detected[sample].keys() if x.startswith('GH') or x.startswith('PL')]))

for host in GH_per_host_median.keys():
    print(f'{host} median cazymes: {np.median(GH_per_host_median[host])}')

# summarize number of cazymes found per sample
with open('results/figures/08-summarize_functions/cazymes_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\tcazymes\n')
    for sample in cazymes_detected:
        out_fh.write(f'{sample}\t{len(cazymes_detected[sample].keys())}\n')

# summarize number of cazymes (not starting with AA or GT) found per sample
with open('results/figures/08-summarize_functions/cazymes_per_sample_noAA_GT.tsv', 'w') as out_fh:
    out_fh.write('sample\tcazymes\n')
    for sample in cazymes_detected:
        out_fh.write(f'{sample}\t{len([x for x in cazymes_detected[sample].keys() if not x.startswith("AA") and not x.startswith("GT")])}\n')

# cum curve of #cazymes
num_iters = 20
with open('results/figures/08-summarize_functions/cazyme_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\tcazymes\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(num_iters):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/{num_iters} and sample size {j}/{len(samples_host)}', end = '\r')
                cum_cazymes = 0
                samples_iter = r.sample(samples_host, sample_size)
                cazymes_detected_iter = set()
                for sample in samples_iter:
                    cazymes_detected_iter.update([x for x in cazymes_detected[sample].keys()])
                cum_cazymes = len(cazymes_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_cazymes}\n')

cazymes_per_host_by_genus = {}
# contains keys as cazyme name and within that a dict with keys as genus
# and values as dicts with number of samples having that cazyme in that genus
# for each host species
genera_list = [os.path.basename(x).split('_cazyme')[0].split('g__')[1] for x in glob.glob('results/figures/08-summarize_functions/by_taxa/genus/*cazyme_matrix.csv')]
for genus in genera_list:
    cazyme_matrix_sub = pd.read_csv(f'results/figures/08-summarize_functions/cazyme_counts/{sample}_cazyme.count', sep='\t')
    # remove all columns with only zeroes
    cazyme_matrix_sub = cazyme_matrix_sub.loc[:, (cazyme_matrix_sub != 0).any(axis=0)]
    cazyme_matrix_sub = cazyme_matrix_sub.set_index('sample')
    for cazy in cazyme_matrix_sub.columns:
        if cazy not in cazymes_per_host_by_genus.keys():
            cazymes_per_host_by_genus[cazy] = {genus: {}}
        else:
            cazymes_per_host_by_genus[cazy][genus] = {}
        samples_detected = [x for x in cazyme_matrix_sub.index if cazyme_matrix_sub.loc[x][cazy] > 0]
        hosts_seen = Counter([host_of_sample(x) for x in samples_detected])
        for host in hosts_seen:
            cazymes_per_host_by_genus[cazy][genus][host] = hosts_seen[host]
    
with open('results/figures/08-summarize_functions/cazymes_per_host_by_genus.tsv', 'w') as out_fh:
    out_fh.write('cazyme\tgenus\thost\tnumber\n')
    for cazyme in cazymes_per_host_by_genus.keys():
        for genus in cazymes_per_host_by_genus[cazyme].keys():
            for host in cazymes_per_host_by_genus[cazyme][genus].keys():
                success = out_fh.write(f'{cazyme}\t{genus}\t{host}\t{cazymes_per_host_by_genus[cazyme][genus][host]}\n')


# 3e. OGs
                
ogs_detected = {sample: {} for sample in samples}
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}')
    for gene in genes_detected[sample]:
        if og_gene_id(gene) is not None and og_gene_id(gene) in gene_og_dict.keys():
            og = gene_og_dict[og_gene_id(gene)]
            if og in ogs_detected[sample].keys():
                ogs_detected[sample][og] += 1
            else:
                ogs_detected[sample][og] = 1

with open('results/figures/08-summarize_functions/ogs_per_sample.tsv', 'w') as out_fh:
    out_fh.write('sample\togs\tcore\n')
    for sample in ogs_detected:
        string_out = f'{sample}\t{len(ogs_detected[sample].keys())}\t'
        core = 0
        for og in ogs_detected[sample]:
            if og in og_coreness:
                if og_coreness[og] == 'core':
                    core += 1
        string_out += f'{core}\n'
        out_fh.write(string_out)

# cum curve of #ogs
num_iters = 20
with open('results/figures/08-summarize_functions/og_cum_curve.tsv', 'w+') as out_fh:
    out_fh.write(f'host\titeration\tsize\togs\n')
    for i, host in enumerate(host_species):
        print(f'working on {host} {i}/5 ')
        samples_host = [x for x in samples if host_of_sample(x) == host and not x.startswith('Am') and not x.startswith('Dr') and not x.startswith('Gr') and not x.startswith('Ac')]
        for iteration in range(num_iters):
            for j, sample_size in enumerate([x for x in range(len(samples_host)) if x > 0]):
                print(f'working om iteration {iteration}/{num_iters} and sample size {j}/{len(samples_host)}', end = '\r')
                cum_ogs = 0
                samples_iter = r.sample(samples_host, sample_size)
                ogs_detected_iter = set()
                for sample in samples_iter:
                    ogs_detected_iter.update([x for x in ogs_detected[sample].keys()])
                cum_ogs = len(ogs_detected_iter)
                n = out_fh.write(f'{host}\t{iteration}\t{sample_size}\t{cum_ogs}\n')

# summarize number of ogs found per sample
# try:
#     og_ko_dict = pickle.load(open('results/figures/08-summarize_functions/og_ko_dict.pkl', 'rb'))
# except FileNotFoundError:
#     pass
# all_og_matrix = pd.read_csv('results/figures/08-summarize_functions/og_matrix.csv')
# ogs_per_host = {}
# for i, og in enumerate(all_og_matrix.columns[1:]):
#     print(f'working on {i}/{len(all_og_matrix.columns[1:])}', end = '\r')
#     found_in = set([host_of_sample(x) for x in all_og_matrix[all_og_matrix[og] > 0]['sample']])
#     if len(found_in) == 0:
#         continue
#     if len(found_in) > 1:
#         status = 'shared'
#     else:
#         status = list(found_in)[0]
#     try:
#         koC = kegg_info_dict[tuple(og_ko_dict[og])[0]]['C']
#     except KeyError:
#         koC = 'Unknown'
#     genus = og.split('--')[0].split('g__')[1]
#     if genus not in ogs_per_host.keys():
#         ogs_per_host[genus] = {}
#     if koC in ogs_per_host[genus].keys():
#         ogs_per_host[genus][koC][status] += 1
#     else:
#         ogs_per_host[genus][koC] = {'shared': 0, 'Apis mellifera': 0, 'Apis cerana': 0, 'Apis dorsata': 0, 'Apis florea': 0, 'Apis andreniformis': 0}
# pickle.dump(ogs_per_host, open('results/figures/08-summarize_functions/ogs_per_host.pkl', 'wb'))
ogs_per_host = pickle.load(open('results/figures/08-summarize_functions/ogs_per_host.pkl', 'rb'))

with open('results/figures/08-summarize_functions/ogs_per_host.tsv', 'w') as out_fh:
    out_fh.write('kos_C\tgenus\tnumber\thost\n')
    for genus in ogs_per_host.keys():
        for koC in ogs_per_host[genus].keys():
            for host in ogs_per_host[genus][koC].keys():
                success = out_fh.write(f'{koC}\t{genus}\t{ogs_per_host[genus][koC][host]}\t{host}\n')

# 4. Make matrices to look at beta diversity

def sample_feature_matrix(feature_name='ko', path=''):
    '''
    feature can be one of gene, ko, cazyme-DRAM, cazyme,
    og, coreness, cluster i.e. one of the columns of the
    df_detected_genes_info table
    read all the sample tables and make a matrix of samples
    as rows and features as columns with the values being
    the number of genes detected in that sample for that
    feature
    '''
    samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]
    '''
        first go through all the sample and collect all the kos
        read the df with list of genes and info for each sample and count the number of genes
        for each feature and add it to the feature matrix if the feature is not present in the
        feature matrix, add it as a column for the row of this sample and set 0 for all other
        samples if the feature is present in the feature matrix, add the count for this sample
        to the column for this feature
    '''
    ko_columns = set()
    for i, sample in enumerate(samples):
        print(f'collecting {feature_name} from {sample} {i}/{len(samples)}', end = '\r')
        df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        ko_columns.update(df[feature_name].unique())
    print(f'finished collecting {feature_name} from all samples')
    feature_matrix = pd.DataFrame()
    feature_matrix['sample'] = samples
    feature_matrix = feature_matrix.set_index('sample')
    feature_names = list(ko_columns)
    feature_matrix = feature_matrix.reindex(columns = feature_names)
    feature_matrix = feature_matrix.fillna(0)
    for i, sample in enumerate(samples):
        print(f'working on {feature_name} for {sample} {i}/{len(samples)}', end = '\r')
        df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        # remove the column unnamed and set gene as index
        df = df.drop(columns = ['Unnamed: 0'])
        df = df.set_index('gene')
        # get the counts for each feature
        feature_counts = df[feature_name].value_counts()
        # add the counts to the feature matrix
        for feature in feature_counts.keys():
            feature_matrix.loc[sample, feature] = feature_counts[feature]
    print(f'finished counting {feature_name} from all samples')
    if path != '':
        # os.makedirs(os.path.dirname(path), exist_ok=True)
        feature_matrix.to_csv(path)
    return feature_matrix

sample_feature_matrix(feature_name='ko', path='results/figures/08-summarize_functions/ko_matrix.csv')
sample_feature_matrix(feature_name='cluster', path='results/figures/08-summarize_functions/cluster_matrix.csv')
sample_feature_matrix(feature_name='cazyme-DRAM', path='results/figures/08-summarize_functions/cazyme-DRAM_matrix.csv')
sample_feature_matrix(feature_name='cazyme', path='results/figures/08-summarize_functions/cazyme_matrix.csv')
sample_feature_matrix(feature_name='og', path='results/figures/08-summarize_functions/og_matrix.csv')

# Next these matrices are analysed in the R script, 11-gene_content.Rmd
'''
It makes the following outputs to be annotated using python:
    from maaslin2:
        results/figures/08-summarize_functions/maaslin2_results_ko/significant_results.tsv
        results/figures/08-summarize_functions/maaslin2_results_ko_pairwise.csv
        results/figures/08-summarize_functions/maaslin2_results_cazy/significant_results.tsv
        results/figures/08-summarize_functions/maaslin2_results_cazy_pairwise.csv
        results/figures/08-summarize_functions/maaslin2_results_og/significant_results.tsv
        results/figures/08-summarize_functions/maaslin2_results_og_pairwise.csv
    # parse the significant results output from the differential abundance analysis
    # add information about the KO to new column(s) and save it adjacent to it
    from the random forest approach:
        results/figures/08-summarize_functions/ko_random_forest.csv
'''
'''
og_species = {}
for sample in samples:
    with open(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv', 'r') as f:
        header = True
        for line in f:
            if header:
                header = False
            else:
                og = line.split(',')[8]
                species = line.split(',')[6]
                if og in og_species.keys():
                    og_species[og].add(species)
                else:
                    og_species[og] = set([species])
for key in og_species.keys():
    og_species[key] = ';'.join(og_species[key])
pickle.dump(og_species, open('results/figures/08-summarize_functions/og_species_dict.pkl', 'wb'))
'''
og_species = pickle.load(open('results/figures/08-summarize_functions/og_species_dict.pkl', 'rb'))
### sanity check
lengths = []
for key in og_species.keys():
    lengths.append(len(og_species[key].split(';')))
Counter(lengths)

'''
og_ko_dict = {}
samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')] 
for i, sample in enumerate(samples):
    print(f'working on {sample} {i}/{len(samples)}', end = '\r')
    with open(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv', 'r') as f:
        for line in f:
            if line.startswith('gene'):
                continue
            else:
                og = line.split(',')[8]
                ko = line.split(',')[2]
                if ko == '':
                    continue
                if og in og_ko_dict.keys():
                    og_ko_dict[og].add(ko)
                else:
                    og_ko_dict[og] = set([ko])
pickle.dump(og_ko_dict, open('results/figures/08-summarize_functions/og_ko_dict.pkl', 'wb'))
'''
og_ko_dict = pickle.load(open('results/figures/08-summarize_functions/og_ko_dict.pkl', 'rb'))
### sanity check
lengths = []
for key in og_ko_dict:
    lengths.append(len(og_ko_dict[key]))
from collections import Counter
Counter(lengths)

def kegg_cat_list_set(og, cat):
    kos = []
    kos_cat = set()
    if og in og_ko_dict.keys():
        kos = og_ko_dict[og]
        for ko in kos:
            kos_cat.add(kegg_info_dict[ko][cat])
    string_of_annos = ';'.join(kos_cat)
    return string_of_annos


# in case needed in the future
'''
gene_positions_dict = {} # contains (contig, start, end) for each gene with mag as keys
# function that can get og of a gene id by parsing or reading pre-parsed information about
# renamed genes used in IG inference by OrthoFinder
# In the orignial files containg genes (.faa), gene ids are a serial number 1 to n per contig
# and in the ffn file the gene ids are the position in the contig
# renaming was done by the rule rename_faa_and_ffn in the mag_phylogenies.smk
# script going from 1 to N simultaneously for the ffn and faa files
# mind that you should only use faa files for this because the ffn files are
# ffn files have genes named not by serial number but by position in the contig
# in the original file
# (waning if the number was different - never came up)
# the input is "results/11_phylogenies/00_genomes/{mag}/{mag_original.ffn"
# the output is "results/11_phylogenies/00_genomes/{mag}/{mag}.ffn"
# make a dict for each mag and save as pickle to be used by og_gene_id function!
for i, mag in enumerate(glob.glob('results/11_phylogenies/00_genomes/*')):
    mag_name = mag.split('/')[-1]
    print(f'working on mag {i}/{len(glob.glob("results/11_phylogenies/00_genomes/*"))}', end = '\r')
    gene_positions_dict[mag_name] = {}
    previous_contig = ''
    i = 1
    with open(f'{mag}/{mag_name}_original.ffn', 'r') as f:
        for line in f:
            if line.startswith('>') and line_o.startswith('>'):
                contig = line.strip().split('>')[1].split(':')[0]
                if previous_contig == '':
                    previous_contig = contig
                if contig != previous_contig:
                    i = 1
                    previous_contig = contig
                else:
                    i += 1
                gene = f'{mag_name}_{contig}_{i}'
                start = line.strip().split('>')[1].split(':')[1].split('-')[0]
                end = line.strip().split('>')[1].split(':')[1].split('-')[1]
                # not adding contig here because it is already in the gene name!
                gene_positions_dict[mag_name][gene] = (start, end)
pickle.dump(gene_positions_dict, open('results/figures/08-summarize_functions/gene_positions_dict.pkl', 'wb'))
'''
mag_gene_positions_dict = pickle.load(open('results/figures/08-summarize_functions/gene_positions_dict.pkl', 'rb'))



# genes_df
# gene_positions_dict['C4-3_6']
# 'C4-3_NODE_4_length_644241_cov_51.191471_609' - was not found in gene_renamed_dict


# bowtie 2 coverage file (from assembly) - using "results/08_gene_content/01_profiling_bowtie2/{sample}_filt_genes.bed" <- 
# from "results/06_metagenomicORFs/{sample}/filt_orfs/{sample}.gff"

# has fewer genes listed than the genes in the scaffolds

# of MAGs which comes from results/07_MAG_binng_QC/03_checkm_results/${sample}/bins/${bin}/genes.faa
# and results/07_MAG_binng_QC/03_checkm_results/${sample}/bins/${bin}/genes.gff

# this is because when checkm runs prodigal (?) only on MAGs it detects some more genes than
# when it was runs on the whole assembly to obtain the filtered ORFs...

# match them using positions and not IDs

# DRAM annotations are done on results/06_metagenomicORFs/${sample}/filt_orfs/${sample}.faa


def sample_feature_matrix_filtered(feature_name='ko', filter_type = '', filter_name = '', path=''):
    '''
    feature can be one of gene, ko, cazyme-DRAM, cazyme,
    og, coreness, cluster i.e. one of the columns of the
    df_detected_genes_info table
    read all the sample tables and make a matrix of samples
    as rows and features as columns with the values being
    the number of genes detected in that sample for that
    feature
    
    This function is useful to collect the KOs/caxymes etc.
    but from a specific species or genus - can be modified to also
    consider core/non-core

    '''
    samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]
    '''
        first go through all the sample and collect all the kos
        read the df with list of genes and info for each sample and count the number of genes
        for each feature and add it to the feature matrix if the feature is not present in the
        feature matrix, add it as a column for the row of this sample and set 0 for all other
        samples if the feature is present in the feature matrix, add the count for this sample
        to the column for this feature
        filter_type = 'genus' or 'species'
    '''
    ko_columns = set()
    for i, sample in enumerate(samples):
        print(f'collecting {feature_name} from {sample} {i}/{len(samples)}', end = '\r')
        df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        df = df[df[filter_type] == filter_name]
        ko_columns.update(df[feature_name].unique())
    print(f'finished collecting {feature_name} from all samples')
    feature_matrix = pd.DataFrame()
    feature_matrix['sample'] = samples
    feature_matrix = feature_matrix.set_index('sample')
    feature_names = list(ko_columns)
    feature_matrix = feature_matrix.reindex(columns = feature_names)
    feature_matrix = feature_matrix.fillna(0)
    for i, sample in enumerate(samples):
        print(f'working on {feature_name} for {sample} {i}/{len(samples)}', end = '\r')
        df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
        df = df[df[filter_type] == filter_name]
        # remove the column unnamed and set gene as index
        df = df.drop(columns = ['Unnamed: 0'])
        df = df.set_index('gene')
        # get the counts for each feature
        feature_counts = df[feature_name].value_counts()
        # add the counts to the feature matrix
        for feature in feature_counts.keys():
            feature_matrix.loc[sample, feature] = feature_counts[feature]
    print(f'finished counting {feature_name} from all samples')
    if path != '':
        # os.makedirs(os.path.dirname(path), exist_ok=True)
        feature_matrix.to_csv(path)
    return feature_matrix

'''
# examples of what can be made not being used at the moment
samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]

genera = set()
for sample in samples:
    df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # skip nans
    df = df.dropna(subset=['genus'])
    genera.update(df['genus'])

species = set()
for sample in samples:
    df = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    # skip nans
    df = df.dropna(subset=['species'])
    species.update(df['species'])

for genus in genera:
    print(f'working on {genus}')
    sample_feature_matrix_filtered(feature_name='ko', filter_type = 'genus', filter_name = f'{genus}', path=f'results/figures/08-summarize_functions/by_taxa/genus/{genus}_ko_matrix.csv')
    sample_feature_matrix_filtered(feature_name='cluster', filter_type = 'genus', filter_name = f'{genus}', path=f'results/figures/08-summarize_functions/by_taxa/genus/{genus}_cluster_matrix.csv')
    sample_feature_matrix_filtered(feature_name='cazyme-DRAM', filter_type = 'genus', filter_name = f'{genus}', path=f'results/figures/08-summarize_functions/by_taxa/genus/{genus}_cazyme-DRAM_matrix.csv')
    sample_feature_matrix_filtered(feature_name='cazyme', filter_type = 'genus', filter_name = f'{genus}', path=f'results/figures/08-summarize_functions/by_taxa/genus/{genus}_cazyme_matrix.csv')
    sample_feature_matrix_filtered(feature_name='og', filter_type = 'genus', filter_name = f'{genus}', path=f'results/figures/08-summarize_functions/by_taxa/genus/{genus}_og_matrix.csv')

for spec in species:
    spec_mod = spec.replace(' ', '--')
    print(f'working on {spec}')
    if not os.path.isfile(f'results/figures/08-summarize_functions/by_taxa/species/{spec_mod}_ko_matrix.csv'):
        sample_feature_matrix_filtered(feature_name='ko', filter_type = 'species', filter_name = f'{spec}', path=f'results/figures/08-summarize_functions/by_taxa/species/{spec_mod}_ko_matrix.csv')
    sample_feature_matrix_filtered(feature_name='cluster', filter_type = 'species', filter_name = f'{spec}', path=f'results/figures/08-summarize_functions/by_taxa/species/{spec_mod}_cluster_matrix.csv')
    sample_feature_matrix_filtered(feature_name='cazyme-DRAM', filter_type = 'species', filter_name = f'{spec}', path=f'results/figures/08-summarize_functions/by_taxa/species/{spec_mod}_cazyme-DRAM_matrix.csv')
    sample_feature_matrix_filtered(feature_name='cazyme', filter_type = 'species', filter_name = f'{spec}', path=f'results/figures/08-summarize_functions/by_taxa/species/{spec_mod}_cazyme_matrix.csv')
    sample_feature_matrix_filtered(feature_name='og', filter_type = 'species', filter_name = f'{spec}', path=f'results/figures/08-summarize_functions/by_taxa/species/{spec_mod}_og_matrix.csv')
'''


'''
# Notes:
    Cayman profiling and bwa results in results/08_gene_content/01_profiling are from mapping against the gene catalog
    but results/08_gene_content/01_profiling_bowtie2 are from mapping against the whole assembly!
'''

# code saved from the other script

samples = [x.split('/')[-1].split('_df_detected_genes_info')[0] for x in glob.glob('results/figures/08-summarize_functions/gene_info_tables/*_df_detected_genes_info.csv')]

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




# kos
df_sig = pd.read_csv(f'results/figures/08-summarize_functions/maaslin2_results_ko/significant_results.tsv', sep='\t')
df_sig['kos_A'] = df_sig['feature'].apply(lambda x: kegg_info_dict[x]['A'] if x in kegg_info_dict.keys() else 'NA')
df_sig['kos_B'] = df_sig['feature'].apply(lambda x: kegg_info_dict[x]['B'] if x in kegg_info_dict.keys() else 'NA')
df_sig['kos_C'] = df_sig['feature'].apply(lambda x: kegg_info_dict[x]['C'] if x in kegg_info_dict.keys() else 'NA')
df_sig.to_csv(f'results/figures/08-summarize_functions/maaslin2_results_ko/significant_results_annotated.tsv', sep='\t')

# 'results/figures/08-summarize_functions/maaslin2_results_ko_pairwise.csv'
df_pairwise = pd.read_csv(f'results/figures/08-summarize_functions/maaslin2_results_ko_pairwise.csv')
df_pairwise['kos_A'] = df_pairwise['feature'].apply(lambda x: kegg_info_dict[x]['A'] if x in kegg_info_dict.keys() else 'NA')
df_pairwise['kos_B'] = df_pairwise['feature'].apply(lambda x: kegg_info_dict[x]['B'] if x in kegg_info_dict.keys() else 'NA')
df_pairwise['kos_C'] = df_pairwise['feature'].apply(lambda x: kegg_info_dict[x]['C'] if x in kegg_info_dict.keys() else 'NA')
df_pairwise.to_csv(f'results/figures/08-summarize_functions/maaslin2_results_ko_pairwise_annotated.tsv', sep='\t')

# cazymes
cazy_anno_info = pd.read_csv('data/cayman_gene_db/20230607_glycan_annotations_cleaned_manually.tsv', sep='\t')
# ORIGIN	FUNCTION_IN_ORIGIN	FUNCTION_AT_DESTINATION_1	FUNCTION_AT_DESTINATION_2	FUNCTION_AT_DESTINATION_3	Family	Subfamily	Activity description	Glycan_annotation	Substrate annotation taken from characterized enzyme description CAZy.org
# feature matches subfamily
df_sig = pd.read_csv(f'results/figures/08-summarize_functions/maaslin2_results_cazyme/significant_results.tsv', sep='\t')
df_sig['Origin'] = df_sig['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['ORIGIN'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_sig['Function_in_origin'] = df_sig['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_IN_ORIGIN'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_sig['Function_at_destination_1'] = df_sig['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_AT_DESTINATION_1'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_sig['Function_at_destination_2'] = df_sig['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_AT_DESTINATION_2'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_sig['Function_at_destination_3'] = df_sig['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_AT_DESTINATION_3'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_sig['Activity_description'] = df_sig['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['Activity description'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_sig['Glycan_annotation'] = df_sig['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['Glycan_annotation'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_sig.to_csv(f'results/figures/08-summarize_functions/maaslin2_results_cazyme/significant_results_annotated.tsv', sep='\t')

# 'results/figures/08-summarize_functions/maaslin2_results_cazyme_pairwise.csv'
cazy_anno_info = pd.read_csv('data/cayman_gene_db/20230607_glycan_annotations_cleaned_manually.tsv', sep='\t')
df_pairwise = pd.read_csv(f'results/figures/08-summarize_functions/maaslin2_results_cazyme_pairwise.csv')
df_pairwise['Origin'] = df_pairwise['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['ORIGIN'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_pairwise['Function_in_origin'] = df_pairwise['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_IN_ORIGIN'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_pairwise['Function_at_destination_1'] = df_pairwise['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_AT_DESTINATION_1'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_pairwise['Function_at_destination_2'] = df_pairwise['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_AT_DESTINATION_2'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_pairwise['Function_at_destination_3'] = df_pairwise['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['FUNCTION_AT_DESTINATION_3'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_pairwise['Activity_description'] = df_pairwise['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['Activity description'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_pairwise['Glycan_annotation'] = df_pairwise['feature'].apply(lambda x: cazy_anno_info[cazy_anno_info['Subfamily'] == x]['Glycan_annotation'].values[0] if x in cazy_anno_info['Subfamily'].values else 'NA')
df_pairwise.to_csv(f'results/figures/08-summarize_functions/maaslin2_results_cazyme_pairwise_annotated.csv')


# OGs
df_sig = pd.read_csv(f'results/figures/08-summarize_functions/maaslin2_results_og/significant_results.tsv', sep='\t')
    # replace double dot in feature with double --
df_sig['feature'] = df_sig['feature'].apply(lambda x: x.replace('..', '--'))
df_sig = df_sig.set_index('feature')
df_sig['coreness'] = df_sig.index.map(og_coreness)
df_sig['og_species'] = df_sig.index.map(og_species)
# df_sig['kos'] = df_sig.index.map(og_ko_dict)
df_sig['kos_A'] = df_sig.index.map(lambda x: kegg_cat_list_set(x, 'A'))
df_sig['kos_B'] = df_sig.index.map(lambda x: kegg_cat_list_set(x, 'B'))
df_sig['kos_C'] = df_sig.index.map(lambda x: kegg_cat_list_set(x, 'C'))
df_sig.to_csv(f'results/figures/08-summarize_functions/maaslin2_results_og/significant_results_annotated.tsv', sep='\t')

df_pairwise = pd.read_csv(f'results/figures/08-summarize_functions/maaslin2_results_og_pairwise.csv')
df_pairwise['feature'] = df_pairwise['feature'].apply(lambda x: x.replace('..', '--'))
df_pairwise = df_pairwise.set_index('feature')
df_pairwise['coreness'] = df_pairwise.index.map(og_coreness)
df_pairwise['og_species'] = df_pairwise.index.map(og_species)
# df_pairwise['kos'] = df_pairwise.index.map(og_ko_dict)
df_pairwise['kos_A'] = df_pairwise.index.map(lambda x: kegg_cat_list_set(x, 'A'))
df_pairwise['kos_B'] = df_pairwise.index.map(lambda x: kegg_cat_list_set(x, 'B'))
df_pairwise['kos_C'] = df_pairwise.index.map(lambda x: kegg_cat_list_set(x, 'C'))
df_pairwise.to_csv(f'results/figures/08-summarize_functions/maaslin2_results_og_pairwise_annotated.csv')

# random forest approach results annotation

df_kruskal = pd.read_csv('results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_table.csv')
df_kruskal = df_kruskal.drop(columns = ['Unnamed: 0'])
# df_kruskal = df_kruskal.set_index('id')
df_kruskal['ko_A'] = df_kruskal['id'].apply(lambda x: kegg_info_dict[x]['A'] if x in kegg_info_dict.keys() else 'NA')
df_kruskal['ko_B'] = df_kruskal['id'].apply(lambda x: kegg_info_dict[x]['B'] if x in kegg_info_dict.keys() else 'NA')
df_kruskal['ko_C'] = df_kruskal['id'].apply(lambda x: kegg_info_dict[x]['C'] if x in kegg_info_dict.keys() else 'NA')
df_kruskal.to_csv('results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_table_annotated.csv')

df_kruskal_imp = pd.read_csv('results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance.csv')
df_kruskal_imp = df_kruskal_imp.drop(columns = ['Unnamed: 0'])
df_kruskal_imp['ko_A'] = df_kruskal_imp['Sample'].apply(lambda x: kegg_info_dict[x]['A'] if x in kegg_info_dict.keys() else 'NA')
df_kruskal_imp['ko_B'] = df_kruskal_imp['Sample'].apply(lambda x: kegg_info_dict[x]['B'] if x in kegg_info_dict.keys() else 'NA')
df_kruskal_imp['ko_C'] = df_kruskal_imp['Sample'].apply(lambda x: kegg_info_dict[x]['C'] if x in kegg_info_dict.keys() else 'NA')
df_kruskal_imp.to_csv('results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv')

df_kruskal_cazy = pd.read_csv('results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_table.csv')
df_kruskal_cazy = df_kruskal_cazy.set_index('Unnamed: 0')
cazy_anno_info = cazy_anno_info.set_index('Subfamily')
df_kruskal_cazy['Origin'] = df_kruskal_cazy.index.map(cazy_anno_info['ORIGIN'])
df_kruskal_cazy['Function_in_origin'] = df_kruskal_cazy.index.map(cazy_anno_info['FUNCTION_IN_ORIGIN'])
df_kruskal_cazy['Function_at_destination_1'] = df_kruskal_cazy.index.map(cazy_anno_info['FUNCTION_AT_DESTINATION_1'])
df_kruskal_cazy['Function_at_destination_2'] = df_kruskal_cazy.index.map(cazy_anno_info['FUNCTION_AT_DESTINATION_2'])
df_kruskal_cazy['Function_at_destination_3'] = df_kruskal_cazy.index.map(cazy_anno_info['FUNCTION_AT_DESTINATION_3'])
df_kruskal_cazy['Family'] = df_kruskal_cazy.index.map(cazy_anno_info['Family'])
df_kruskal_cazy['Activity_description'] = df_kruskal_cazy.index.map(cazy_anno_info['Activity description'])
df_kruskal_cazy['Glycan_annotation'] = df_kruskal_cazy.index.map(cazy_anno_info['Glycan_annotation'])
df_kruskal_cazy['Substrate_annotation'] = df_kruskal_cazy.index.map(cazy_anno_info['Substrate annotation taken from characterized enzyme description CAZy.org'])
df_kruskal_cazy.to_csv('results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_table_annotated.csv')

df_kruskal_cazy_imp = pd.read_csv('results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_importance.csv')
df_kruskal_cazy_imp = df_kruskal_cazy_imp.set_index('Sample')
df_kruskal_cazy_imp['Origin'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['ORIGIN'])
df_kruskal_cazy_imp['Function_in_origin'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['FUNCTION_IN_ORIGIN'])
df_kruskal_cazy_imp['Function_at_destination_1'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['FUNCTION_AT_DESTINATION_1'])
df_kruskal_cazy_imp['Function_at_destination_2'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['FUNCTION_AT_DESTINATION_2'])
df_kruskal_cazy_imp['Function_at_destination_3'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['FUNCTION_AT_DESTINATION_3'])
df_kruskal_cazy_imp['Family'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['Family'])
df_kruskal_cazy_imp['Activity_description'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['Activity description'])
df_kruskal_cazy_imp['Glycan_annotation'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['Glycan_annotation'])
df_kruskal_cazy_imp['Substrate_annotation'] = df_kruskal_cazy_imp.index.map(cazy_anno_info['Substrate annotation taken from characterized enzyme description CAZy.org'])
df_kruskal_cazy_imp.to_csv('results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_importance_annotated.csv')

df_kruskal_og = pd.read_csv('results/figures/08-summarize_functions/og_compared_kruskal_wallis_table.csv')
df_kruskal_og['id'] = df_kruskal_og['id'].apply(lambda x: x.replace('..', '--'))
df_kruskal_og['coreness'] = df_kruskal_og['id'].map(og_coreness)
df_kruskal_og['og_species'] = df_kruskal_og['id'].map(og_species)
df_kruskal_og['kos'] = df_kruskal_og['id'].map(og_ko_dict)
df_kruskal_og['kos_A'] = df_kruskal_og['id'].apply(lambda x: kegg_cat_list_set(x, 'A'))
df_kruskal_og['kos_B'] = df_kruskal_og['id'].apply(lambda x: kegg_cat_list_set(x, 'B'))
df_kruskal_og['kos_C'] = df_kruskal_og['id'].apply(lambda x: kegg_cat_list_set(x, 'C'))
df_kruskal_og.to_csv('results/figures/08-summarize_functions/og_compared_kruskal_wallis_table_annotated.csv')

df_kruskal_og_imp = pd.read_csv('results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance.csv')
df_kruskal_og_imp['Sample'] = df_kruskal_og_imp['Sample'].apply(lambda x: x.replace('..', '--'))
df_kruskal_og_imp['coreness'] = df_kruskal_og_imp['Sample'].map(og_coreness)
df_kruskal_og_imp['og_species'] = df_kruskal_og_imp['Sample'].map(og_species)
df_kruskal_og_imp['kos'] = df_kruskal_og_imp['Sample'].map(og_ko_dict)
df_kruskal_og_imp['kos_A'] = df_kruskal_og_imp['Sample'].apply(lambda x: kegg_cat_list_set(x, 'A'))
df_kruskal_og_imp['kos_B'] = df_kruskal_og_imp['Sample'].apply(lambda x: kegg_cat_list_set(x, 'B'))
df_kruskal_og_imp['kos_C'] = df_kruskal_og_imp['Sample'].apply(lambda x: kegg_cat_list_set(x, 'C'))
df_kruskal_og_imp.to_csv('results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance_annotated.csv')


os.makedirs('results/figures/visualize_temp/KO_collections', exist_ok=True)

'''
community-level: across samples, not subset by taxonomy
'''
for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    os.makedirs(f'results/figures/visualize_temp/KO_collections/by_sample/', exist_ok=True)
    # only write the list of kos withouth the header
    df_sample['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_sample/{sample}_kos.csv', index=False, header=False)
    # run annotator
    print(f'printed {sample} kos')

os.makedirs(f'results/figures/visualize_temp/microbeannotator_out/', exist_ok=True)
input_list = ' '.join([f'results/figures/visualize_temp/KO_collections/by_sample/{sample}_kos.csv' for sample in samples])
print(f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_sample')

'''
community-level: across samples, subset microbial genus
'''

for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    for genus in df_sample['genus'].unique():
        df_genus = df_sample[df_sample['genus'] == genus]
        os.makedirs(f'results/figures/visualize_temp/KO_collections/by_genus/{genus}', exist_ok=True)
        # only write the list of kos withouth the header
        df_genus['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_genus/{genus}/{sample}_kos.csv', index=False, header=False)
        # run annotator
        print(f'printed {sample} kos')
    print(f'printed {sample} kos')

# this tells you the module completeness of each pathway
"python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i results/figures/visualize_temp/KOs.csv -p test_C3-2_snod"
# actual code used for this step got lost unfortunately


'''
community-level: across samples, subset by microbial species
'''


for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    for species in df_sample['species'].unique():
        df_species = df_sample[df_sample['species'] == species]
        os.makedirs(f'results/figures/visualize_temp/KO_collections/by_species/{species}', exist_ok=True)
        # only write the list of kos withouth the header
        df_species['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_species/{species}/{sample}_kos.csv', index=False, header=False)
        # run annotator
        print(f'printed {sample} kos for {species}', end='\r')
    print(f'done for {sample}')

os.makedirs(f'results/figures/visualize_temp/microbeannotator_out/by_species', exist_ok=True)


# 
"python ../MinPath.py -ko ../../../results/figures/visualize_temp/KOs.csv -report KOs_test_gil.ko.minpath -details KOs_test_gil.ko.minpath.details"

'''
 4.3 How to read the MinPath report file?
    e.g, demo.ko.minpath
    ...
    path 00030 kegg n/a  naive 1  minpath 1  fam0  42  fam-found  18  name  Pentose phosphate pathway
    path 00031 kegg n/a  naive 1  minpath 0  fam0  12  fam-found  2  name  Inositol metabolism
    ...
    1) path 00030, path 00031 are the KEGG pathway IDs (if your input file has fig families, the pathways will then be SEED subsystems)
    2) kegg n/a, indicates pathway reconstruction of your input dataset if not available (note: this information is only available for the genomes annotated 
	 in KEGG database); for input with fig families, KEGG is replaced by SEED 
    3) naive 1 or 0: the pathway is reconstructed, or not, by the naive mapping approach
    4) minpath 1 or 0: the pathway is kept, or removed by MinPath
    5) fam0: the total number of families involved in the corresponding pathway
    6) fam-found: the total number of involved families that are annotated
    7) name: the description of the corresponding pathway (subsystem)
  
 4.4 How to read the MinPath detailed report file?
    This report file lists all the pathways found MinPath, and a list of families each pathway includes
    e.g. demo.ko.minpath.details
    ...
    path 00010 fam0 56 fam-found 27 # Glycolysis / Gluconeogenesis
       K00001 hits 6 # E1.1.1.1, adh
       K00002 hits 1 # E1.1.1.4, adh
    ...
    1) path 00010 is the KEGG pathway ID, and fam0 and fam-found are the same as in 4.3
    2) this pathway includes families, K00001 and K00002, and so on 
       (here K numbers are used for KEGG families, and FIG ids for FIG families)
    3) family K00001 has 6 hits (i.e., 6 proteins/reads are annotated as this family)
'''

# do this each time for species that are to be compared across samples

os.makedirs('results/figures/visualize_temp/KO_collections/by_species/combos/bombi_mellis_mellifer/', exist_ok=True)
species_selected = ['Bombilactobacillus mellifer', 'Bombilactobacillus mellis']
for species in species_selected:
    for file in glob.glob(f'results/figures/visualize_temp/KO_collections/by_species/{species}/*'):
        species_mod = species.replace(' ', '__')
        shutil.copy(file, f'results/figures/visualize_temp/KO_collections/by_species/combos/bombi_mellis_mellifer/{species_mod}_{os.path.basename(file)}')
    print(f'done for {species}')
input_list = ' '.join([x for x in glob.glob(f'results/figures/visualize_temp/KO_collections/by_species/combos/bombi_mellis_mellifer/*')])
command = f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_species/combos/bombi_mellis_mellifer'
cluster_header = '''#!/bin/bash\n
######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name bombi_mellis_mellifer
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 150G
#SBATCH --time 05:00:00 
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.err
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.out\n

source ~/.bashrc
conda activate 20230313_scripts_env\n

cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison\n
'''

with open('results/figures/visualize_temp/KO_collections/by_species_bombi_mellis_mellifer_command.sh', 'w+') as out_fh:
    success = out_fh.write(cluster_header)
    success = out_fh.write(command)


'''
summarize functions of MAGs from functions in all samples
    all_kos is a file containing two columns, first is the MAG name and second is the KO
    this comes from the code mentioned in mag_phylogenies.smk - collects this info from concat dram output
'''

# first get the list of KOs and name the file to have the handmade species name and MAG name (do not add any file extensions)

os.makedirs(f'results/figures/visualize_temp/KO_collections/by_mag/', exist_ok=True)

mags = set()
with open('results/09_MAGs_collection/functions_list/all_kos.txt', 'r') as in_fh:
    for line in in_fh:
        mag = line.split('\t')[0]
        ko = line.split('\t')[1].strip()
        if mag in mags:
            continue
        else:
            mags.add(mag)

for mag in mags:
    print(mag)
    with open(f'results/figures/visualize_temp/KO_collections/by_mag/{mag}', 'w+') as out_fh:
        with open('results/09_MAGs_collection/functions_list/all_kos.txt', 'r') as in_fh:
            for line in in_fh:
                mag = line.split('\t')[0]
                if mag == mag:
                    ko = line.split('\t')[1].strip()
                    success = out_fh.write(f'{ko}\n')
                else:
                    continue

os.makedirs(f'results/figures/visualize_temp/microbeannotator_out/', exist_ok=True)
input_list = ' '.join([f'results/figures/visualize_temp/KO_collections/by_mag/{mag}' for mag in mags])
command = f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_mag/by_mag_all'
# print(f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_mag')
cluster_header = '''#!/bin/bash\n
######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name all_mags_by_mag
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 250G
#SBATCH --time 15:00:00 
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.err
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.out\n

source ~/.bashrc
conda activate 20230313_scripts_env\n

cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison\n
'''

with open('results/figures/visualize_temp/KO_collections/by_mag_command.sh', 'w+') as out_fh:
    success = out_fh.write(cluster_header)
    success = out_fh.write(command)


# only bombilacto
input_list = ' '.join([f'results/figures/visualize_temp/KO_collections/by_mag/{mag}' for mag in mags if "Bombilacto" in mag])
command = f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_mag_bombi'
cluster_header = '''#!/bin/bash\n
######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name cdhit-clustering
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 150G
#SBATCH --time 05:00:00 
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotator_by_mag_bombi.err
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotator_by_mag_bombi.out\n

source ~/.bashrc
conda activate 20230313_scripts_env\n

cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison\n
'''

with open('results/figures/visualize_temp/KO_collections/by_mag_command_bombi.sh', 'w+') as out_fh:
    success = out_fh.write(cluster_header)
    success = out_fh.write(command)


'''
quick processing of the dram file to map gene id to function name when found
'''

df_functions = pd.read_csv('data/dram_ function_heatmap_form.tsv', sep ='\t')
dict_func = {}
for func_list, name in zip(df_functions['function_ids'], df_functions['function_name']):
    func_list = func_list.split(',')
    for func in func_list:
        if func not in dict_func.keys():
            dict_func[func] = name
with open('data/dram_function_map.tsv', 'w+') as out_fh:
    for func, name in dict_func.items():
        success = out_fh.write(f'{func}\t{name}\n')

# # summarizing gene info with positions and rpkms
# for sample in samples:
#     df_info = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')

# make a list of KOs with the title as the "genome name"
# run with cluster both
# this genome name can be any kind of grouping that we want to compare
os.makedirs('results/figures/visualize_temp/KO_collections', exist_ok=True)

'''
community-level: across samples, not subset by taxonomy
'''
for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    os.makedirs(f'results/figures/visualize_temp/KO_collections/by_sample/', exist_ok=True)
    # only write the list of kos withouth the header
    df_sample['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_sample/{sample}_kos.csv', index=False, header=False)
    # run annotator
    print(f'printed {sample} kos')

os.makedirs(f'results/figures/visualize_temp/microbeannotator_out/', exist_ok=True)
input_list = ' '.join([f'results/figures/visualize_temp/KO_collections/by_sample/{sample}_kos.csv' for sample in samples])
print(f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_sample')

'''
community-level: across samples, subset microbial genus
'''

for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    for genus in df_sample['genus'].unique():
        df_genus = df_sample[df_sample['genus'] == genus]
        os.makedirs(f'results/figures/visualize_temp/KO_collections/by_genus/{genus}', exist_ok=True)
        # only write the list of kos withouth the header
        df_genus['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_genus/{genus}/{sample}_kos.csv', index=False, header=False)
        # run annotator
        print(f'printed {sample} kos')
    print(f'printed {sample} kos')

# this tells you the module completeness of each pathway
"python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i results/figures/visualize_temp/KOs.csv -p test_C3-2_snod"
# actual code used for this step got lost unfortunately


'''
community-level: across samples, subset by microbial species
'''


for i, sample in enumerate(samples):
    df_sample = pd.read_csv(f'results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv')
    df_sample = df_sample[df_sample['ko'].notnull()]
    for species in df_sample['species'].unique():
        df_species = df_sample[df_sample['species'] == species]
        os.makedirs(f'results/figures/visualize_temp/KO_collections/by_species/{species}', exist_ok=True)
        # only write the list of kos withouth the header
        df_species['ko'].to_csv(f'results/figures/visualize_temp/KO_collections/by_species/{species}/{sample}_kos.csv', index=False, header=False)
        # run annotator
        print(f'printed {sample} kos for {species}', end='\r')
    print(f'done for {sample}')

os.makedirs(f'results/figures/visualize_temp/microbeannotator_out/by_species', exist_ok=True)


# 
"python ../MinPath.py -ko ../../../results/figures/visualize_temp/KOs.csv -report KOs_test_gil.ko.minpath -details KOs_test_gil.ko.minpath.details"

'''
 4.3 How to read the MinPath report file?
    e.g, demo.ko.minpath
    ...
    path 00030 kegg n/a  naive 1  minpath 1  fam0  42  fam-found  18  name  Pentose phosphate pathway
    path 00031 kegg n/a  naive 1  minpath 0  fam0  12  fam-found  2  name  Inositol metabolism
    ...
    1) path 00030, path 00031 are the KEGG pathway IDs (if your input file has fig families, the pathways will then be SEED subsystems)
    2) kegg n/a, indicates pathway reconstruction of your input dataset if not available (note: this information is only available for the genomes annotated 
	 in KEGG database); for input with fig families, KEGG is replaced by SEED 
    3) naive 1 or 0: the pathway is reconstructed, or not, by the naive mapping approach
    4) minpath 1 or 0: the pathway is kept, or removed by MinPath
    5) fam0: the total number of families involved in the corresponding pathway
    6) fam-found: the total number of involved families that are annotated
    7) name: the description of the corresponding pathway (subsystem)
  
 4.4 How to read the MinPath detailed report file?
    This report file lists all the pathways found MinPath, and a list of families each pathway includes
    e.g. demo.ko.minpath.details
    ...
    path 00010 fam0 56 fam-found 27 # Glycolysis / Gluconeogenesis
       K00001 hits 6 # E1.1.1.1, adh
       K00002 hits 1 # E1.1.1.4, adh
    ...
    1) path 00010 is the KEGG pathway ID, and fam0 and fam-found are the same as in 4.3
    2) this pathway includes families, K00001 and K00002, and so on 
       (here K numbers are used for KEGG families, and FIG ids for FIG families)
    3) family K00001 has 6 hits (i.e., 6 proteins/reads are annotated as this family)
'''

# run minpath for each species and overall across all but per sample


# do this each time for species that are to be compared across samples

os.makedirs('results/figures/visualize_temp/KO_collections/by_species/combos/bombi_mellis_mellifer/', exist_ok=True)
species_selected = ['Bombilactobacillus mellifer', 'Bombilactobacillus mellis']
for species in species_selected:
    for file in glob.glob(f'results/figures/visualize_temp/KO_collections/by_species/{species}/*'):
        species_mod = species.replace(' ', '__')
        shutil.copy(file, f'results/figures/visualize_temp/KO_collections/by_species/combos/bombi_mellis_mellifer/{species_mod}_{os.path.basename(file)}')
    print(f'done for {species}')
input_list = ' '.join([x for x in glob.glob(f'results/figures/visualize_temp/KO_collections/by_species/combos/bombi_mellis_mellifer/*')])
command = f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_species/combos/bombi_mellis_mellifer'
cluster_header = '''#!/bin/bash\n
######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name bombi_mellis_mellifer
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 150G
#SBATCH --time 05:00:00 
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.err
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.out\n

source ~/.bashrc
conda activate 20230313_scripts_env\n

cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison\n
'''

with open('results/figures/visualize_temp/KO_collections/by_species_bombi_mellis_mellifer_command.sh', 'w+') as out_fh:
    success = out_fh.write(cluster_header)
    success = out_fh.write(command)


'''
summarize functions of MAGs from functions in all samples
    all_kos is a file containing two columns, first is the MAG name and second is the KO
    this comes from the code mentioned in mag_phylogenies.smk - collects this info from concat dram output
'''

# first get the list of KOs and name the file to have the handmade species name and MAG name (do not add any file extensions)

os.makedirs(f'results/figures/visualize_temp/KO_collections/by_mag/', exist_ok=True)

mags = set()
with open('results/09_MAGs_collection/functions_list/all_kos.txt', 'r') as in_fh:
    for line in in_fh:
        mag = line.split('\t')[0]
        ko = line.split('\t')[1].strip()
        if mag in mags:
            continue
        else:
            mags.add(mag)

for mag in mags:
    print(mag)
    with open(f'results/figures/visualize_temp/KO_collections/by_mag/{mag}', 'w+') as out_fh:
        with open('results/09_MAGs_collection/functions_list/all_kos.txt', 'r') as in_fh:
            for line in in_fh:
                mag = line.split('\t')[0]
                if mag == mag:
                    ko = line.split('\t')[1].strip()
                    success = out_fh.write(f'{ko}\n')
                else:
                    continue

os.makedirs(f'results/figures/visualize_temp/microbeannotator_out/', exist_ok=True)
input_list = ' '.join([f'results/figures/visualize_temp/KO_collections/by_mag/{mag}' for mag in mags])
command = f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_mag/by_mag_all'
# print(f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_mag')
cluster_header = '''#!/bin/bash\n
######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name all_mags_by_mag
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 250G
#SBATCH --time 15:00:00 
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.err
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotater.out\n

source ~/.bashrc
conda activate 20230313_scripts_env\n

cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison\n
'''

with open('results/figures/visualize_temp/KO_collections/by_mag_command.sh', 'w+') as out_fh:
    success = out_fh.write(cluster_header)
    success = out_fh.write(command)


# only bombilacto
input_list = ' '.join([f'results/figures/visualize_temp/KO_collections/by_mag/{mag}' for mag in mags if "Bombilacto" in mag])
command = f'python3 scripts/MicrobeAnnotator/microbeannotator/pipeline/ko_mapper.py -i {input_list} -p results/figures/visualize_temp/microbeannotator_out/by_mag_bombi'
cluster_header = '''#!/bin/bash\n
######### SLURM OPTIONS
#SBATCH --partition cpu
#SBATCH --account pengel_spirit
#SBATCH --job-name cdhit-clustering
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 150G
#SBATCH --time 05:00:00 
#SBATCH --error /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotator_by_mag_bombi.err
#SBATCH --output /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/figures/visualize_temp/KO_collections/logs/microbe_annotator_by_mag_bombi.out\n

source ~/.bashrc
conda activate 20230313_scripts_env\n

cd /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison\n
'''

with open('results/figures/visualize_temp/KO_collections/by_mag_command_bombi.sh', 'w+') as out_fh:
    success = out_fh.write(cluster_header)
    success = out_fh.write(command)


'''
quick processing of the dram file to map gene id to function name when found
'''

df_functions = pd.read_csv('data/dram_ function_heatmap_form.tsv', sep ='\t')
dict_func = {}
for func_list, name in zip(df_functions['function_ids'], df_functions['function_name']):
    func_list = func_list.split(',')
    for func in func_list:
        if func not in dict_func.keys():
            dict_func[func] = name
with open('data/dram_function_map.tsv', 'w+') as out_fh:
    for func, name in dict_func.items():
        success = out_fh.write(f'{func}\t{name}\n')