#!/usr/bin/env python3

import os

def get_sdp_phylotype(genome_check):
    # return a dict of phylotype ans sdp info
    info = dict()
    with open(os.path.join(os.getcwd(), metafile), "r") as metafile_fh:
        for line in metafile_fh:
            genome = line.split()[0]
            phylotype = line.split()[1]
            sdp = line.split()[2]
            if genome == genome_check:
                info["sdp"] = sdp
                info["phylotype"] = phylotype
                return(info)

def get_gene_id(OG_id):
    with open(orthofile, "r") as orthofile_fh:
        for line in orthofile_fh:
            line = line.strip()
            split_line = line.split()
            OG_id = split_line.pop(0).split(":")[0]
            for gene in split_line:
                if genome == gene.split("_")[0]:
                    return(gene)

genome = snakemake.wildcards["genome"]
metafile = os.path.join(os.getcwd(), snakemake.input["metafile"])
bedfile = os.path.join(os.getcwd(), snakemake.input["bedfile"])
sdp = get_sdp_phylotype(genome)["sdp"]
phylotype = get_sdp_phylotype(genome)["phylotype"]

for file in snakemake.input["sdp_coords"]:
    core_cov_dir = "/".join(file.split("/")[:-1])
    continue

coord_file = os.path.join(os.getcwd(), core_cov_dir, phylotype+"_corecov_coord.txt")
corecov_file = os.path.join(os.getcwd(), core_cov_dir, phylotype+"_corecov.txt")
orthofile = os.path.join(os.getcwd(), snakemake.params["ortho_dir"], phylotype+"_single_ortho_filt.txt")

outfile = os.path.join(os.getcwd(), snakemake.output["core_red_bed"])

#Open the *coord-file, containing the terminus coverage per sample,
# save names of samples with at least 20x ter-cov
high_cov_samples = dict()
with open(coord_file, "r") as coord_fh:
    for line in coord_fh:
        if line.startswith("SDP"):
            continue
        else:
            sample = line.split()[1]
            ter_cov = float(line.split()[2])
            if ter_cov >= 20:
                high_cov_samples[sample] = ter_cov

#Open the *corecov-file, containing the gene coverage per core gene,
# for each sample. For the samples with more than 20x ter-cov, check if the
# gene coverage is at least 10x. If it is, save the gene-id.
low_cov_genes = dict()
with open(corecov_file, "r") as corecov_fh:
    for line in corecov_fh:
        if line.startswith("SDP"):
            continue
        else:
            sample = line.split()[1]
            if sample in high_cov_samples.keys():
                gene_cov = float(line.split()[4])
                OG_id = line.split()[2]
                gene_id = get_gene_id(OG_id)
                if gene_cov < 10:
                    low_cov_genes[gene_id] = 1

#Open the bed-file containing the core-genes, and create a reduced version
# containing genes with sufficient coverage (i.e. the real core-genes)
with open(bedfile, "r") as bedfile_fh:
    with open(outfile, "w+") as outfile_fh:
        for line in bedfile_fh:
            gene_id = line.split()[3]
            if gene_id in low_cov_genes:
                continue
            else:
                outfile_fh.write(line)
