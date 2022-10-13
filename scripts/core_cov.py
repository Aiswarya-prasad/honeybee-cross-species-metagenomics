#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess

def make_g_cluster_dict(path):
    """
    return the a dict with values as cluster of a genome which are the keys
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
            cluster = line.split("\t")[11]
            g_dict[genome] = cluster
        return(g_dict)

def make_cluster_g_dict(path):
    """
    return a dict with keys as magOTU clusters with list of corresponding MAGs
    as values by reading from checkpoint output where First column is MAGs names
    and 12th column is affiliated magOTU cluster
    """
    cluster_list_dict = {}
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
            cluster = line.split("\t")[11]
            # only include groups of interest!
            if cluster not in cluster_list_dict.keys():
                cluster_list_dict[cluster] = [genome]
            else:
                cluster_list_dict[cluster].append(genome)
    return(cluster_list_dict)

def get_bedcov_genome(bedfile, bamfile):
    bedcov_run = subprocess.run(['samtools', 'bedcov', bedfile, bamfile], universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bedcov = bedcov_run.stdout
    bedcov = bedcov.strip()
    bedcov_list = bedcov.split("\n")
    genome_gene_cov = dict()
    for line in bedcov_list:
        line = line.strip()
        split_line = line.split("\t")
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

def print_to_file(outfile, cluster,cluster_bam_genefam_cov):
    fh_outfile = open(outfile, 'a')
    cluster_ref_pos = OG_ref_pos[cluster]
    sorted_pos = sorted(cluster_ref_pos.items(), key=lambda x: x[1]) #Tuple (OG-pos pairs)
    bamfile = "what"
    for ele in sorted_pos:
        OG_id = ele[0][0:-1] #Trim off colon from id
        start_pos = str(ele[1])
        if bamfile.find('/') != -1: #Check if a directory was provided for the bam-file. If so, remove the path (so the sample-id in the outfile wont contain the path)
            split_filename = bamfile.split('/')
            sample = split_filename[-1][0:-4]
        else:
            sample = bamfile[0:-4]
        line_out = [cluster, sample, OG_id, start_pos, str(cluster_bam_genefam_cov[cluster][bamfile][ele[0]])]
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
requiredNamed.add_argument('--info',metavar="db_metafile",required=True, help="File detailing genome-id and cluster affiliation", action="store") # same as 06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv
requiredNamed.add_argument('--ortho', metavar='orthofinder_file', required=True, help="Filtered single-copy ortholog gene file (orthofinder format)", action="store")
requiredNamed.add_argument('--beddir', metavar="bedfile_dir", required=True, help="Directory containing bed-files for genomes in db", action="store")
requiredNamed.add_argument('--bamfile', metavar="bamfile_notlistanymore", help="List of bam-files, one line per file", action="store")
requiredNamed.add_argument('--outfile', metavar="out_file", required=False, help="Output file containing coverage of genefamilies in each cluster", action="store")

database_metafile = args.info
ortho_file = args.ortho
bedfile_dir = args.beddir
outfile = args.outfile

#Read the database metafile, store cluster-affiliation and cluster-ref ids in dictionaries
cluster_ref = make_cluster_g_dict(database_metafile)
genome_cluster_ref = make_g_cluster_dict(database_metafile)

#Read the orthofinder-file, get genome-ids per cluster and gene-family members per cluster. Store OG-family affiliation for all gene-ids.
with open(ortho_file, "r") as fh_orthofile:
    cluster_genomes = dict() #Genome-ids contained within each cluster in orthofinder file (cluster - genome_id - 1)
    gene_OG = dict() #Family affiliation of all genes in orthofinder file
    cluster_OG_genes = dict() #Gene-members for each OG, for each cluster (cluster - OG_id -gene_list)
    for line in fh_orthofile:
        line = line.strip()
        split_line = line.split()
        OG_id = split_line.pop(0)
        for gene in split_line:
            gene_OG[gene] = OG_id
            split_gene = gene.split('_')
            genome_id = "_".join(split_gene[:-1])
            cluster = genome_cluster_ref[genome_id]
            if cluster not in cluster_genomes:
                cluster_genomes[cluster] = dict()
            cluster_genomes[cluster][genome_id] = 1
            if cluster not in cluster_OG_genes:
                cluster_OG_genes[cluster] = dict()
            if OG_id not in cluster_OG_genes[cluster]:
                cluster_OG_genes[cluster][OG_id] = list()
            cluster_OG_genes[cluster][OG_id].append(gene)

# For each OG,
    # check if it exists in each genoeme
    # get the start position of that OG in that genome (start position and chromosome?)

cluster_ref = dict()
with open(ref_info, "r") as ref_info_fh:
    for line in ref_info_fh:
        if line.startswith("genome"):
            continue
        line = line.strip()
        split_line = line.split("\t")
        genome_id = split_line[0]
        ref_status = split_line[2]
        cluster = split_line[1]
        if ref_status == "1":
            print(genome_id)
            cluster_ref[genome_id] = cluster

#Read the bed-files for the cluster reference genomes, get the start-position for each gene-family
OG_ref_pos = dict()
for cluster in cluster_genomes.keys():
    ref_genome = cluster_ref[cluster]
    bedfile = ref_genome + ".bed"
    check_file_exists(bedfile)
    fh_bedfile = open(bedfile)
    if cluster not in OG_ref_pos:
        OG_ref_pos[cluster] = dict()
    for line in fh_bedfile:
        line = line.strip()
        split_line = line.split("\t")
        gene_id = split_line[3]
        start_pos = split_line[1]
        if gene_id not in gene_OG: continue
        OG_id = gene_OG[gene_id]
        OG_ref_pos[cluster][OG_id] = int(start_pos)
    fh_bedfile.close()
os.chdir(cwd)

#Prepare for outfile, print the header
outfile = out_name + "_corecov.txt"
outfile = os.path.join(out_dir, outfile)
fh_outfile = open(outfile, 'w')
header = ["cluster","Sample","OG", "Ref_pos", "Coverage"]
header_str = "\t".join(header)
fh_outfile.write(header_str + "\n")
fh_outfile.close()

#Get the cluster-coverage of all gene-families in orthofile, for all listed bam-files
cluster_bam_genefam_cov = dict()  #cluster - bamfile - OG - summed_cov
gene_cov = dict()
cluster_list = list(cluster_genomes.keys())
cluster_list.sort()
for cluster in cluster_list:
    print("Working on cluster:", cluster)
    cluster_bam_genefam_cov[cluster] = dict()
    genomes = cluster_genomes[cluster]
    print("\tProcessing bamfile:",bamfile)
    cluster_bam_genefam_cov[cluster][bamfile] = dict()
    #For each genome, get coverage on OG genes, store in gene-cov dict
    for genome in genomes.keys():
        bedfile = bedfile_dir + '/' + genome + '.bed'
        check_file_exists(bedfile)
        fh_bedfile = open(bedfile)
        print("\t\tGenome:", genome)
        gene_cov_genome = get_bedcov_genome(bedfile, bamfile)
        gene_cov.update(gene_cov_genome)
    print("\tSumming up cluster-coverage for each gene-family..")
    OG_fams = cluster_OG_genes[cluster]
    for OG_id in OG_fams.keys():
        OG_genes = OG_fams[OG_id]
        genefam_cov = 0
        for gene in OG_genes:
            genefam_cov += gene_cov[gene]
            cluster_bam_genefam_cov[cluster][bamfile][OG_id]= round(genefam_cov,2)
    print_to_file(outfile, cluster, cluster_bam_genefam_cov)







"""
OLD CODE
"""
# def get_bedcov_genome(bedfile, bamfile):
#     bedcov_run = subprocess.run(['samtools', 'bedcov', bedfile, bamfile], universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     bedcov = bedcov_run.stdout
#     bedcov = bedcov.strip()
#     bedcov_list = bedcov.split("\n")
#     genome_gene_cov = dict()
#     for line in bedcov_list:
#         line = line.strip()
#         split_line = line.split("\t")
#         gene_id = split_line[3]
#         if gene_id not in gene_OG: continue
#         OG_id = gene_OG[gene_id]
#         gene_start = int(split_line[1])
#         gene_end = int(split_line[2])
#         gene_length = gene_end - gene_start
#         mapped_reads = int(split_line[4])
#         gene_cov = mapped_reads/gene_length
#         genome_gene_cov[gene_id] = gene_cov
#     return(genome_gene_cov)
#
# def check_file_exists(file):
#     try:
#         fh = open(file)
#     except:
#         print("Cant find/open this file: ", file)
#         print("Exiting script!")
#         exit()
#
# def print_to_file(outfile, sdp,sdp_bam_genefam_cov):
#     fh_outfile = open(outfile, 'a')
#     sdp_ref_pos = OG_ref_pos[sdp]
#     sorted_pos = sorted(sdp_ref_pos.items(), key=lambda x: x[1]) #Tuple (OG-pos pairs)
#     for bamfile in bamfiles:
#         for ele in sorted_pos:
#             OG_id = ele[0][0:-1] #Trim off colon from id
#             start_pos = str(ele[1])
#             if bamfile.find('/') != -1: #Check if a directory was provided for the bam-file. If so, remove the path (so the sample-id in the outfile wont contain the path)
#                 split_filename = bamfile.split('/')
#                 sample = split_filename[-1][0:-4]
#             else:
#                 sample = bamfile[0:-4]
#             line_out = [sdp, sample, OG_id, start_pos, str(sdp_bam_genefam_cov[sdp][bamfile][ele[0]])]
#             line_out_str = "\t".join(line_out)
#
#             fh_outfile.write(line_out_str + "\n")
#     fh_outfile.close()
#
# #Check software requirements
# python_version_major = sys.version_info.major
# if python_version_major != 3:
#     print("This script required python3! You are running:")
#     print("You are using Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
#     print("Exiting script!")
#     sys.exit(1)
# python_version_minor = sys.version_info.minor
# if python_version_minor < 6:
#     print("You are running and old version of python (< 3.6). This may cause trouble with the subprocess module")
#     print("Exiting script!")
#     sys.exit(1)
# try:
#     samtools_check = subprocess.run(['samtools', 'bedcov'], universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# except:
#     print("ERROR: subprocess module calling 'samtools bedcov' returned a bad exit code")
#     print("Check the samtools installation from the command-line:")
#     print("samtools bedcov")
#     print("Exiting script!")
#     sys.exit(1)
#
# #Parse input options
# parser = argparse.ArgumentParser()
# requiredNamed = parser.add_argument_group('required arguments')
# requiredNamed.add_argument('-d',metavar="db_metafile",required=True, help="File detailing genome-id and SDP affiliation")
# requiredNamed.add_argument('-l', metavar="bamfile_list",required=True, help="List of bam-files, one line per file")
# requiredNamed.add_argument('-g', metavar='orthofinder_file', required=True, help="Filtered single-copy ortholog gene file (orthofinder format)")
# requiredNamed.add_argument('-b', metavar="bedfile_dir", required=True, help="Directory containing bed-files for genomes in db")
# parser.add_argument('-o', metavar="out_dir", required=False, help="Directory to write outputs to")
# parser.add_argument('-n', metavar="out_name", required=False, help="name to append to output file")
# args = vars(parser.parse_args())
# if (args['d']):
#     database_metafile = args['d']
# if (args['b']):
#     bamfile_list = args['l']
# if (args['g']):
#     ortho_file = args['g']
# if (args['b']):
#     bedfile_dir = args['b']
# if not os.path.isdir(bedfile_dir):
#    print("The specified bed-file directory does not exist:", bedfile_dir)
#    print("Exiting script")
#    exit()
# if (args['o']):
#     out_dir = args['o']
#     if not os.path.isdir(out_dir):
#         if os.path.isabs(out_dir):
#             os.makedirs(out_dir)
#         else:
#             out_dir = os.path.join(os.getcwd(), out_dir)
#     else:
#         print(out_dir + " already exists. Existing files will be over-written")
# else:
#     out_dir = os.getcwd()
# print("outputs will be written to " + out_dir)
# if (args['n']):
#     out_name = args['n']
# else:
#     if ortho_file.find('/') != -1: #Check if a directory was provided for the orthofinder-file. If so, remove the path (so outfile will be printed in the run-dir)
#         split_filename = ortho_file.split('/')
#         name_split = split_filename[-1]
#     split_filename = name_split.split('_')
#     out_name = split_filename[0]
#
# #Read the database metafile, store sdp-affiliation and sdp-ref ids in dictionaries
# db_dict = dict() #genome - sdp
# sdp_ref = dict()
# check_file_exists(database_metafile)
# fh_metafile = open(database_metafile)
# for line in fh_metafile:
#     line = line.strip()
#     split_line = line.split("\t")
#     genome_id = split_line[0]
#     sdp = split_line[2]
#     ref_stat = split_line[3]
#     db_dict[genome_id] = sdp
#     if ref_stat == "Ref":
#         sdp_ref[sdp] = genome_id
# fh_metafile.close()
#
# #Read the bam-file list, save filenames in list
# bamfiles = list()
# check_file_exists(bamfile_list)
# fh_bamlist = open(bamfile_list)
# for line in fh_bamlist:
#     line = line.strip()
#     if len(line) == 0: continue #Check for empty lines
#     check_file_exists(line)
#     bamfiles.append(line)
# fh_bamlist.close()
#
# #Read the orthofinder-file, get genome-ids per SDP and gene-family members per SDP. Store OG-family affiliation for all gene-ids.
# check_file_exists(ortho_file)
# fh_orthofile = open(ortho_file)
# sdp_genomes = dict() #Genome-ids contained within each SDP in orthofinder file (SDP - genome_id - 1)
# gene_OG = dict() #Family affiliation of all genes in orthofinder file
# sdp_OG_genes = dict() #Gene-members for each OG, for each SDP (SDP - OG_id -gene_list)
# for line in fh_orthofile:
#     line = line.strip()
#     split_line = line.split()
#     OG_id = split_line.pop(0)
#     for gene in split_line:
#         gene_OG[gene] = OG_id
#         split_gene = gene.split('_')
#         genome_id = split_gene[0]
#         sdp = db_dict[genome_id]
#         if sdp not in sdp_genomes:
#             sdp_genomes[sdp] = dict()
#         sdp_genomes[sdp][genome_id] = 1
#         if sdp not in sdp_OG_genes:
#             sdp_OG_genes[sdp] = dict()
#         if OG_id not in sdp_OG_genes[sdp]:
#             sdp_OG_genes[sdp][OG_id] = list()
#         sdp_OG_genes[sdp][OG_id].append(gene)
# fh_orthofile.close()
#
# #Read the bed-files for the SDP reference genomes, get the start-position for each gene-family
# cwd = os.getcwd()
# os.chdir(bedfile_dir)
# OG_ref_pos = dict()
# for sdp in sdp_genomes.keys():
#     ref_genome = sdp_ref[sdp]
#     bedfile = ref_genome + ".bed"
#     check_file_exists(bedfile)
#     fh_bedfile = open(bedfile)
#     if sdp not in OG_ref_pos:
#         OG_ref_pos[sdp] = dict()
#     for line in fh_bedfile:
#         line = line.strip()
#         split_line = line.split("\t")
#         gene_id = split_line[3]
#         start_pos = split_line[1]
#         if gene_id not in gene_OG: continue
#         OG_id = gene_OG[gene_id]
#         OG_ref_pos[sdp][OG_id] = int(start_pos)
#     fh_bedfile.close()
# os.chdir(cwd)
#
# #Prepare for outfile, print the header
# outfile = out_name + "_corecov.txt"
# outfile = os.path.join(out_dir, outfile)
# fh_outfile = open(outfile, 'w')
# header = ["SDP","Sample","OG", "Ref_pos", "Coverage"]
# header_str = "\t".join(header)
# fh_outfile.write(header_str + "\n")
# fh_outfile.close()
#
# #Get the SDP-coverage of all gene-families in orthofile, for all listed bam-files
# sdp_bam_genefam_cov = dict()  #sdp - bamfile - OG - summed_cov
# gene_cov = dict()
# sdp_list = list(sdp_genomes.keys())
# sdp_list.sort()
# for sdp in sdp_list:
#     print("Working on SDP:", sdp)
#     sdp_bam_genefam_cov[sdp] = dict()
#     genomes = sdp_genomes[sdp]
#     for bamfile in bamfiles:
#         print("\tProcessing bamfile:",bamfile)
#         sdp_bam_genefam_cov[sdp][bamfile] = dict()
#         #For each genome, get coverage on OG genes, store in gene-cov dict
#         for genome in genomes.keys():
#             bedfile = bedfile_dir + '/' + genome + '.bed'
#             check_file_exists(bedfile)
#             fh_bedfile = open(bedfile)
#             print("\t\tGenome:", genome)
#             gene_cov_genome = get_bedcov_genome(bedfile, bamfile)
#             gene_cov.update(gene_cov_genome)
#         print("\tSumming up SDP-coverage for each gene-family..")
#         OG_fams = sdp_OG_genes[sdp]
#         for OG_id in OG_fams.keys():
#             OG_genes = OG_fams[OG_id]
#             genefam_cov = 0
#             for gene in OG_genes:
#                 genefam_cov += gene_cov[gene]
#                 sdp_bam_genefam_cov[sdp][bamfile][OG_id]= round(genefam_cov,2)
#     print_to_file(outfile, sdp, sdp_bam_genefam_cov)
