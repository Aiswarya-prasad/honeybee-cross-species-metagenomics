import os
from Bio import SeqIO
import pandas as pd

mag_info = pd.read_csv('results/09_MAGs_collection/MAGs_metadata_summary.tsv', sep = '\t', index_col = 0)
mag_info.columns
# get handmade species names of each magOTU, mag,magotu,handmadespeciesname
handmade_info = pd.read_csv('config/Species_MAG_Cluster-names.txt', sep = '\t')
# magotu = '53_1'
def get_handmade_speciesname(magotu):
    if pd.isna(magotu):
        return 'NA'
    else:
        name = handmade_info[handmade_info['cluster'] == magotu]['MAG_species_name_final']
        if len(name) == 0:
            return 'NA'
        else:
            return name.values[0]

class genome:
    def __init__(self, header):
        self.id = ''
        self.org_id = '' # id including the GB_ or RS_ part not informative for MAGs
        self.type = ''
        self.tax_string = ''
        self.species = ''
        self.genus = ''
        self.family = ''
        self.order = ''
        self.tax_class = ''
        self.phylum = ''
        self.domain = ''

        if len(header.split(' ')) > 1 and ('GB_' in header or 'RS_' in header):
            if 'RS_' in header:
                self.id = header.split(' ')[0].split('RS_')[1]
                self.org_id = header.split(' ')[0]
            else:
                self.id = header.split(' ')[0].split('GB_')[1]
                self.org_id = header.split(' ')[0]
            self.type = 'GTDB_ref'
            self.tax_string = ' '.join(header.split(' ')[1:])
            self.species = self.tax_string.split('s__')[1]
            self.genus = self.tax_string.split('g__')[1].split(';')[0]
            self.family = self.tax_string.split('f__')[1].split(';')[0]
            self.order = self.tax_string.split('o__')[1].split(';')[0]
            self.tax_class = self.tax_string.split('c__')[1].split(';')[0]
            self.phylum = self.tax_string.split('p__')[1].split(';')[0]
            self.domain = self.tax_string.split('d__')[1].split(';')[0]
        else:
            self.id = header
            self.org_id = self.id
            self.type = 'MAG'
            species = mag_info.loc[self.id, 'Species']
            if species != 's__':
                species = str(species)
                self.species = species.split('s__')[1]
            genus = mag_info.loc[self.id, 'Genus']
            if genus != 'g__':
                genus = str(genus)
                self.genus = genus.split('g__')[1]
            family = mag_info.loc[self.id, 'Family']
            if family != 'f__':
                family = str(family)
                self.family = family.split('f__')[1]
            order = mag_info.loc[self.id, 'Order']
            if order != 'o__':
                order = str(order)
                self.order = order.split('o__')[1]
            tax_class = mag_info.loc[self.id, 'Class']
            if tax_class != 'c__':
                tax_class = str(tax_class)
                self.tax_class = tax_class.split('c__')[1]
            phylum = mag_info.loc[self.id, 'Phylum']
            if phylum != 'p__':
                phylum = str(phylum)
                self.phylum = phylum.split('p__')[1]
            domain = mag_info.loc[self.id, 'Domain']
            if domain != 'd__':
                domain = str(domain)
                self.domain = domain.split('d__')[1]

'''
This script subsets the GTDB MSA to only include one MAG per species
and explores and manually as specified in this script, subsets the
references that are relevant otherwise there are way too many lines in
the alignment
The input is the unizpped version of the MSA found in:
 gtdb_output/align/*.bac120.msa.fasta.gz

Also make metadata tables containing type and taxa information
'''

mag_phyla = set()
mag_classes = set()
mag_orders = set()
mag_families = set()
mag_genera = set()
mag_species = set()
genomes = []
# Read in the MSA
for record in SeqIO.parse('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs.bac120.msa.fasta', 'fasta'):
    genomes.append(genome(record.description))


for gen in genomes:
    if gen.type == 'MAG':
        mag_phyla.add(gen.phylum)
        mag_classes.add(gen.tax_class)
        mag_orders.add(gen.order)
        mag_families.add(gen.family)
        mag_genera.add(gen.genus)
        mag_species.add(gen.species)

ids_to_keep = set()
for gen in genomes:
    if gen.family in mag_families:
        ids_to_keep.add(gen.org_id)
len(ids_to_keep)
out_info_table = f'results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family-info.tsv'
with open(out_info_table, 'w+') as f:
    f.write(f'id\ttype\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmagotu\n')
    for gen in genomes:
        if gen.org_id in ids_to_keep:
            if gen.type == 'MAG':
                magOTU = mag_info.loc[gen.id]['magOTU']
                name = get_handmade_speciesname(magOTU)
            else:
                name = gen.species
            line = '\t'.join([gen.org_id, gen.type, gen.phylum, gen.tax_class, gen.order, gen.family, gen.genus, gen.species, str(name)])
            f.write(f'{line}\n')
out_msa = 'results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_family.bac120.msa.fasta'
with open(out_msa, 'w+') as out:
    for record in SeqIO.parse('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs.bac120.msa.fasta', 'fasta'):
        if genome(record.description).id in ids_to_keep:
            out_header = record.description
            out.write('>' + out_header + '\n' + str(record.seq) + '\n')

ids_to_keep = set()
for gen in genomes:
    if gen.order in mag_orders:
        ids_to_keep.add(gen.org_id)
len(ids_to_keep)
out_info_table = 'results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_order-info.tsv'
with open(out_info_table, 'w+') as f:
    f.write(f'id\ttype\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmagotu\n')
    for gen in genomes:
        if gen.org_id in ids_to_keep:
            if gen.type == 'MAG':
                magOTU = mag_info.loc[gen.id]['magOTU']
                name = get_handmade_speciesname(magOTU)
            else:
                name = gen.species
            line = '\t'.join([gen.org_id, gen.type, gen.phylum, gen.tax_class, gen.order, gen.family, gen.genus, gen.species, str(name)])
            f.write(f'{line}\n')
out_msa = 'results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_order.bac120.msa.fasta'
with open(out_msa, 'w') as out:
    for record in SeqIO.parse('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs.bac120.msa.fasta', 'fasta'):
        if genome(record.description).id in ids_to_keep:
            out_header = record.description
            out.write('>' + out_header + '\n' + str(record.seq) + '\n')

ids_to_keep = set()
for gen in genomes:
    if gen.tax_class in mag_classes:
        ids_to_keep.add(gen.org_id)
len(ids_to_keep)
out_info_table = 'results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_class-info.tsv'
with open(out_info_table, 'w+') as f:
    f.write(f'id\ttype\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmagotu\n')
    for gen in genomes:
        if gen.org_id in ids_to_keep:
            if gen.type == 'MAG':
                magOTU = mag_info.loc[gen.id]['magOTU']
                name = get_handmade_speciesname(magOTU)
            else:
                name = gen.species
            line = '\t'.join([gen.org_id, gen.type, gen.phylum, gen.tax_class, gen.order, gen.family, gen.genus, gen.species, str(name)])
            f.write(f'{line}\n')
out_msa = 'results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs_class.bac120.msa.fasta'
with open(out_msa, 'w') as out:
    for record in SeqIO.parse('results/11_phylogenies/04_MAGs_gtdb/20230313_MAGs.bac120.msa.fasta', 'fasta'):
        if genome(record.description).id in ids_to_keep:
            out_header = record.description
            out.write('>' + out_header + '\n' + str(record.seq) + '\n')