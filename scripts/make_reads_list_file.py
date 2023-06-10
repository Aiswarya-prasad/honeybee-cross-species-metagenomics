import os
import sys
import glob
import yaml

"""
use this script to make the list of raw files to be accessed via config.yaml
"""
# get this list from config file
file_path = "config/raw_file_paths.yaml"
with open("config/config.yaml", "r") as config_fh:
    config_dict = yaml.safe_load(config_fh)
    samples = config_dict["SAMPLES_KE"] + config_dict["SAMPLES_INDIA"] + config_dict["SAMPLES_MY"]
    rawpaths = config_dict["RawdataPaths"]

paths_dict = {}
for sample in samples:
    paths_dict[sample] = {}
    paths_dict[sample]["R1"] = []
    paths_dict[sample]["R2"] = []

for dir in rawpaths:
    files = glob.glob(dir+"/*")
    for file in files:
        if file.endswith("fastq.gz"):
            for sample in samples:
                if sample in file:
                    if "_R1" in file:
                        paths_dict[sample]["R1"].append(file)
                    if "_R2" in file:
                        paths_dict[sample]["R2"].append(file)
    else:
        continue
    
with open(file_path, "w") as out_fh:
    out_fh.write(yaml.dump(paths_dict))
# raw_paths_dict_all = yaml.safe_load(open("config/raw_file_paths.yaml", "r"))
# len(get_list_of_values(get_renamed_input_files(raw_paths_dict_all)))