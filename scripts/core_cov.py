#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess

def make_g_magOTU_dict(path):
    """
    !!! modified to work with ref info (minimal information) rather than checkpoint output
    return the a dict with values as magOTU of a genome which are the keys
    by reading the checkpoint file from the path provided
    """
    g_dict = {}
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            magOTU = line.split("\t")[1]
            g_dict[genome] = magOTU
        return(g_dict)

def make_magOTU_g_dict(path):
    """
    !!! modified to work with ref info (minimal information) rather than checkpoint output
    return a dict with keys as magOTU magOTUs with list of corresponding MAGs
    as values by reading from checkpoint output where First column is MAGs names
    and 12th column is affiliated magOTU magOTU
    """
    magOTU_list_dict = {}
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r", encoding='utf-8-sig') as info_fh:
        for line in info_fh:
            line = line.strip()
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            magOTU = line.split("\t")[1]
            # only include groups of interest!
            if magOTU not in magOTU_list_dict.keys():
                magOTU_list_dict[magOTU] = [genome]
            else:
                magOTU_list_dict[magOTU].append(genome)
    return(magOTU_list_dict)

def get_bedcov_genome(bedfile, bamfile):
    if not os.path.exists(bamfile+".bai"):
        subprocess.run(['samtools', 'index', bamfile])
    bedcov_run = subprocess.run(['samtools', 'bedcov', bedfile, bamfile], universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bedcov = bedcov_run.stdout
    bedcov = bedcov.strip()
    bedcov_list = bedcov.split("\n")
    genome_gene_cov = dict()
    for line in bedcov_list:
        line = line.strip()
        split_line = line.split("\t")
        print(split_line)
        gene_id = split_line[3]
        if gene_id not in gene_OG: continue
        OG_id = gene_OG[gene_id]
        gene_start = int(split_line[1])
        gene_end = int(split_line[2])
        gene_length = gene_end - gene_start
        mapped_reads = int(split_line[4])
        gene_cov = mapped_reads/gene_length
        genome_gene_cov[gene_id] = gene_cov
    return(genome_gene_cov)

def check_file_exists(file):
    if os.path.exists(file):
        pass
    else:
        print(f"Cant find this file: {file}")
        print("Exiting script!")
        exit()

def print_to_file(outfile, magOTU,magOTU_genefam_cov, sample):
    fh_outfile = open(outfile, 'a')
    magOTU_ref_pos = OG_ref_pos[magOTU]
    sorted_pos = sorted(magOTU_ref_pos.items(), key=lambda x: x[1]) #Tuple (OG-pos pairs)
    for ele in sorted_pos:
        OG_id = ele[0][0:-1] #Trim off colon from id
        start_pos = str(ele[1])
        line_out = [magOTU, sample, OG_id, start_pos, str(magOTU_genefam_cov[magOTU][ele[0]])]
        line_out_str = "\t".join(line_out)
        fh_outfile.write(line_out_str + "\n")
    fh_outfile.close()

#Check software requirements
python_version_major = sys.version_info.major
if python_version_major != 3:
    print("This script required python3! You are running:")
    print("You are using Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    print("Exiting script!")
    sys.exit(1)
python_version_minor = sys.version_info.minor
if python_version_minor < 6:
    print("You are running and old version of python (< 3.6). This may cause trouble with the subprocess module")
    print("Exiting script!")
    sys.exit(1)
try:
    samtools_check = subprocess.run(['samtools', 'bedcov'], universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except:
    print("ERROR: subprocess module calling 'samtools bedcov' returned a bad exit code")
    print("Check the samtools installation from the command-line:")
    print("samtools bedcov")
    print("Exiting script!")
    sys.exit(1)

#Parse input options
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--info',metavar="db_metafile",required=True, help="File from checkpoint output with all mag details", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
requiredNamed.add_argument('--ortho', metavar='orthofinder_file', required=True, help="Filtered single-copy ortholog gene file (orthofinder format)", action="store")
requiredNamed.add_argument('--beddir', metavar="bedfiles_dir", required=True, help="Directory containing bed-files for genomes in db", action="store")
requiredNamed.add_argument('--bamfile', metavar="bamfile_notlistanymore", help="List of bam-files, one line per file", action="store")
requiredNamed.add_argument('--group', metavar="group", required=True, help="Name of sample being evaluated", action="store")
requiredNamed.add_argument('--sample', metavar="sample", required=True, help="Name of genus / group being evaluated", action="store")
requiredNamed.add_argument('--outfile', metavar="out_file", required=True, help="Output file containing coverage of genefamilies in each magOTU", action="store")
args = parser.parse_args()

database_metafile = args.info
ortho_file = args.ortho
bedfiles_dir = args.beddir
outfile = args.outfile
bamfile = args.bamfile
group = args.group
sample = args.sample

#Read the database metafile, store magOTU-affiliation and magOTU-ref ids in dictionaries
magOTU_ref = make_magOTU_g_dict(database_metafile)
genome_magOTU_ref = make_g_magOTU_dict(database_metafile)

#Read the orthofinder-file, get genome-ids per magOTU and gene-family members per magOTU. Store OG-family affiliation for all gene-ids.
with open(ortho_file, "r") as fh_orthofile:
    magOTU_genomes = dict() #Genome-ids contained within each magOTU in orthofinder file (magOTU - genome_id - 1)
    gene_OG = dict() #Family affiliation of all genes in orthofinder file
    magOTU_OG_genes = dict() #Gene-members for each OG, for each magOTU (magOTU - OG_id -gene_list)
    for line in fh_orthofile:
        line = line.strip()
        split_line = line.split()
        OG_id = split_line.pop(0)
        for gene in split_line:
            gene_OG[gene] = OG_id
            split_gene = gene.split('_')
            genome_id = "_".join(split_gene[:-1])
            magOTU = genome_magOTU_ref[genome_id]
            if magOTU not in magOTU_genomes:
                magOTU_genomes[magOTU] = dict()
            magOTU_genomes[magOTU][genome_id] = 1
            if magOTU not in magOTU_OG_genes:
                magOTU_OG_genes[magOTU] = dict()
            if OG_id not in magOTU_OG_genes[magOTU]:
                magOTU_OG_genes[magOTU][OG_id] = list()
            magOTU_OG_genes[magOTU][OG_id].append(gene)

magOTU_ref = dict()
with open(database_metafile, "r") as database_metafile_fh:
    for line in database_metafile_fh:
        if line.startswith("genome"):
            continue
        line = line.strip()
        split_line = line.split("\t")
        genome_id = split_line[0]
        ref_status = split_line[2]
        magOTU = split_line[1]
        if ref_status == "1":
            magOTU_ref[magOTU] = genome_id

#Read the bed-files for the magOTU reference genomes, get the start-position for each gene-family

# For each OG,
    # check if it exists in each ref genoeme of a magOTU / SDP
    # get the start position of that OG in that genome (start position and chromosome?)
    # if it does not exist in that genome (probably most complete genome of that SDP),
    # OG_ref_pos contains for each magOTU a dict of OG_id and start_pos
    # OGs missing in the ref genome will just not be in the OG_ref_pos dictionary

OG_ref_pos = dict()
for magOTU in magOTU_genomes.keys():
    ref_genome = magOTU_ref[magOTU]
    bedfile = os.path.join(bedfiles_dir, ref_genome + ".bed")
    check_file_exists(bedfile)
    with open(bedfile, "r") as fh_bedfile:
        if magOTU not in OG_ref_pos:
            OG_ref_pos[magOTU] = dict()
        for line in fh_bedfile:
            line = line.strip()
            print(line)
            split_line = line.split("\t")
            print(split_line)
            gene_id = split_line[3]
            start_pos = split_line[1]
            if gene_id not in gene_OG: continue
            OG_id = gene_OG[gene_id]
            OG_ref_pos[magOTU][OG_id] = int(start_pos)

#Prepare for outfile, print the header
fh_outfile = open(outfile, 'w')
with open(outfile, 'w') as fh_outfile:
    header = ["magOTU","Sample","OG", "Ref_pos", "Coverage"]
    header_str = "\t".join(header)
    fh_outfile.write(header_str + "\n")

#Get the magOTU-coverage of all gene-families in orthofile, for all listed bam-files
magOTU_genefam_cov = dict()  #magOTU - bamfile - OG - summed_cov
gene_cov = dict()
magOTU_list = list(magOTU_genomes.keys())
magOTU_list.sort()

for magOTU in magOTU_list:
    print("Working on magOTU:", magOTU)
    magOTU_genefam_cov[magOTU] = dict()
    genomes = magOTU_genomes[magOTU]
    print("\tProcessing bamfile:", bamfile)
    magOTU_genefam_cov[magOTU] = dict()
    #For each genome, get coverage on OG genes, store in gene-cov dict
    for genome in genomes.keys():
        bedfile = os.path.join(bedfiles_dir, genome + ".bed")
        check_file_exists(bedfile)
        fh_bedfile = open(bedfile)
        print("\t\tGenome:", genome)
        print(f"\t\tEstimating bedcov for:{bedfile}, {bamfile}")
        gene_cov_genome = get_bedcov_genome(bedfile, bamfile)
        gene_cov.update(gene_cov_genome)
    print("\tSumming up magOTU-coverage for each gene-family..")
    OG_fams = magOTU_OG_genes[magOTU]
    for OG_id in OG_fams.keys():
        OG_genes = OG_fams[OG_id]
        genefam_cov = 0
        for gene in OG_genes:
            genefam_cov += gene_cov[gene]
            magOTU_genefam_cov[magOTU][OG_id]= round(genefam_cov,2)
    print_to_file(outfile, magOTU, magOTU_genefam_cov, sample)
