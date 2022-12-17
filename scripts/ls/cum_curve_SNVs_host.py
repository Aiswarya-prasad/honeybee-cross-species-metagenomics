#!/usr/bin/env python3

import os
import random

table_length = os.path.join(os.getcwd(), snakemake.input["lengths"])
freq = os.path.join(os.getcwd(), snakemake.input["freq"])
metadata = os.path.join(os.getcwd(), snakemake.input["metadata"])

cum_curve = os.path.join(os.getcwd(), snakemake.output["cum_curve"])
cum_curve_by_group = os.path.join(os.getcwd(), snakemake.output["cum_curve_by_group"])
cum_curve_by_colony = os.path.join(os.getcwd(), snakemake.output["cum_curve_by_colony"])

nb_curves = snakemake.params["nb_curves"]

sdp = snakemake.wildcards["sdp"]

sample_host = {}
sample_location = {}
sample_colony = {}
host_ids = {}
host_ids_group = {}
host_ids_colony = {}

def get_colony(sampleID):
    """
    This function may need to be changed if sample names are very different
    from the names in 211200_Meliifera_cerana_comparison
    """
    if "Am" == sampleID[:2] or "Ac" == sampleID[:2]:
        colony = sampleID[2:4]
    else:
        colony = sampleID[:2]
    return(colony)

with open(metadata, "r") as metadata_fh:
    header = metadata_fh.readline()
    # print(f"ignoring header line: {header}")
    for line in metadata_fh:
        # if line.startswith("ID"):
        # FOR SOME REASON THIS DOESNT WORK!
        #     continue
        sample = line.split(",")[0]
        host = line.split(",")[4]
        sample_host[sample] = host
        location = line.split(",")[6]
        colony = get_colony(sample)
        # sample_location[sample] = location
        sample_colony[sample] = colony
        host_ids[sample] = host
        host_ids_group[sample] = host+"_"+location
        host_ids_colony[sample] = host+"_"+colony

core_lengths = {}

with open(table_length, "r") as lengths_fh:
    for line in lengths_fh:
        line = line.strip("\n")
        sdp_key = line.split("\t")[0]
        length = line.split("\t")[2]
        core_lengths[sdp_key] = length

var_pos_sample = {}

with open(freq, "r") as freq_fh:
    for line in freq_fh:
        line = line.strip("\n")
        if line.startswith("SNP_info"):
            samples = line.split("\t")[1:]
            for sample in samples:
                var_pos_sample[sample] = set()
        else:
            if not line:
                # empty lines are considered false
                continue
            msnp_desc = line.split("\t")[0]
            msnp_pos = msnp_desc.split(":")[1]
            if len(line.split("\t")[1:]) != len(samples):
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("Number of samples is not equal to number of entries for m/snp: "+msnp_desc+"!!!")
                print("Results may be INVALID!")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            for num, sample in zip(line.split("\t")[1:], samples):
                num = float(num)
                # snps and mnps are counted as the same
                # ie. 1 unique variant so length not considered
                # this is why sets can be used
                # using dict or list here is too time consuming to run!
                # not using:
                # msnp_len = len(msnp_desc.split(":")[2].split(">")[0])
                if num > 0 and num < 1:
                    if msnp_pos not in var_pos_sample[sample]:
                        var_pos_sample[sample].add(msnp_pos)
print(f"#####################")
print(f"Working on cumulative curves per sample")
print(f"#####################")
with open(cum_curve, "w+") as cum_curve_fh:
    for host_id in set(host_ids.values()):
        host_samples = []
        nb_samples = 0
        for sample in host_ids.keys():
            if host_ids[sample] == host_id:
                if sample in var_pos_sample:
                    nb_samples += 1
                    host_samples.append(sample)
        print(f"#####################")
        print(f"\nhost id: {host_id}")
        print(f"#samples: {nb_samples}")
        print(f"sample names: {host_samples}")
        print(f"#####################")
        for curve_nb in range(1, nb_curves+1): # want numbers from 1 to n
            curve_id = "Curve_"+str(curve_nb)
            print(f"Cureve_id: {curve_id}")
            for i in range(1, nb_samples+1): # want numbers from 1 to n
                var_sites = 0
                pos_seen = []
                random_sample_list = random.sample(host_samples, k = i)
                print(f"***iteration with {i} samples: {random_sample_list}")
                for sample in random_sample_list:
                    new_sites = 0
                    discovered = set()
                    discovered = var_pos_sample[sample].difference(pos_seen)
                    new_sites = len(discovered)
                    pos_seen.extend(discovered)
                    var_sites += new_sites
                    print(f"Number of new sites in {sample}: {new_sites}")
                print(f"Total sites: {var_sites}")
                fraction = var_sites/float(core_lengths[sdp])*100
                cum_curve_fh.write(f"{host_id}\t{curve_id}\t{i}\t{fraction}\t{sdp}\n")
print(f"#####################")
print(f"Working on cumulative curves per group")
print(f"#####################")
with open(cum_curve_by_group, "w+") as cum_curve_by_group_fh:
    for host_id in set(host_ids_group.values()):
        host_samples = []
        nb_samples = 0
        for sample in host_ids_group.keys():
            if host_ids_group[sample] == host_id:
                if sample in var_pos_sample:
                    nb_samples += 1
                    host_samples.append(sample)
        # print(f"#####################")
        # print(f"\nhost id: {host_id}")
        # print(f"#samples: {nb_samples}")
        # print(f"sample names: {host_samples}")
        # print(f"#####################")
        for curve_nb in range(1, nb_curves+1): # want numbers from 1 to n
            curve_id = "Curve_"+str(curve_nb)
            # print(f"Cureve_id: {curve_id}")
            for i in range(1, nb_samples+1): # want numbers from 1 to n
                var_sites = 0
                pos_seen = []
                random_sample_list = random.sample(host_samples, k = i)
                # print(f"***iteration with {i} samples: {random_sample_list}")
                for sample in random_sample_list:
                    new_sites = 0
                    discovered = set()
                    discovered = var_pos_sample[sample].difference(pos_seen)
                    new_sites = len(discovered)
                    pos_seen.extend(discovered)
                    var_sites += new_sites
                #     print(f"Number of new sites in {sample}: {new_sites}")
                # print(f"Total sites: {var_sites}")
                fraction = var_sites/float(core_lengths[sdp])*100
                cum_curve_by_group_fh.write(f"{host_id}\t{curve_id}\t{i}\t{fraction}\t{sdp}\n")
print(f"#####################")
print(f"Working on cumulative curves per colony")
print(f"#####################")
with open(cum_curve_by_colony, "w+") as cum_curve_by_colony_fh:
    for host_id in set(host_ids_colony.values()):
        host_samples = []
        nb_samples = 0
        for sample in host_ids_colony.keys():
            if host_ids_colony[sample] == host_id:
                if sample in var_pos_sample:
                    nb_samples += 1
                    host_samples.append(sample)
        # print(f"#####################")
        # print(f"\nhost id: {host_id}")
        # print(f"#samples: {nb_samples}")
        # print(f"sample names: {host_samples}")
        # print(f"#####################")
        for curve_nb in range(1, nb_curves+1): # want numbers from 1 to n
            curve_id = "Curve_"+str(curve_nb)
            # print(f"Cureve_id: {curve_id}")
            for i in range(1, nb_samples+1): # want numbers from 1 to n
                var_sites = 0
                pos_seen = []
                random_sample_list = random.sample(host_samples, k = i)
                # print(f"***iteration with {i} samples: {random_sample_list}")
                for sample in random_sample_list:
                    new_sites = 0
                    discovered = set()
                    discovered = var_pos_sample[sample].difference(pos_seen)
                    new_sites = len(discovered)
                    pos_seen.extend(discovered)
                    var_sites += new_sites
                #     print(f"Number of new sites in {sample}: {new_sites}")
                # print(f"Total sites: {var_sites}")
                fraction = var_sites/float(core_lengths[sdp])*100
                cum_curve_by_colony_fh.write(f"{host_id}\t{curve_id}\t{i}\t{fraction}\t{sdp}\n")
