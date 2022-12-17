#!/usr/bin/env python3

import os

counter = 0

genomes = []

phylotypes = ["api", "bapis", "bifido", "bom", "com", "firm4", "firm5", "fper", "gilli", "lkun", "snod"]

for x in os.listdir("../database/faa_files"):
    if x.endswith(".faa"):
        x = x.split(".")[0]
        genomes.append(x)

phylo_dict = {key: None for key in genomes}

with open("../database/genome_db_210402_metafile.txt") as metafile:
    for line in metafile:
        line = line.split()
        locus_tag = line[0]
        phylotype = line[1]
        phylo_dict[locus_tag] = phylotype

for key in phylo_dict.keys():
    if phylo_dict[key] is None:
        print(key + " not listed in metafile")
        exit()

for phylo in phylotypes:
    directory = os.path.join("../database/faa_files", phylo)
    if not os.path.exists(directory):
        os.mkdir(directory)

for genome in genomes:
    path = "../database/faa_files/"
    os.rename(path+genome+".faa", path+phylo_dict[genome]+"/"+genome+".faa")
