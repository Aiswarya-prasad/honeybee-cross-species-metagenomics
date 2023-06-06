def convertToMb(string):
    """
    This function can convert text in the form
    xxG to mb
    If it does not end with G, it returns the string
    It does not handle other cases of invalid input
    """
    if string.endswith("G"):
        number = int(string.split("G")[0])
        return(number*1000)
    else:
        return(string)

def convertToSec(string):
    """
    This function can convert text in the form
    D-hh:mm:ss to seconds
    D - # days
    hh # hours
    mm # mins
    ss # secs
    """
    days = string.split("-")[0]
    hrs = string.split("-")[1].split(":")[0]
    min = string.split("-")[1].split(":")[1]
    sec = string.split("-")[1].split(":")[2]
    total = int(sec)
    total = total + 60*int(min)
    total = total + 60*60*int(hrs)
    total = total + 24*60*60*int(days)
    return(total)

def get_all_input_files(paths_dict):
    """
    get list of all paths
    """
    all_paths = []
    for sample in paths_dict.keys():
        R1_reads = paths_dict[sample]["R1"]
        R2_reads = paths_dict[sample]["R2"]
        all_reads = [x for x in chain(R1_reads, R2_reads)]
        all_paths = [x for x in chain(all_paths, all_reads)]
    return(all_paths)


def get_input_file(paths_dict, sample, read, lane, run):
    """
    runs are identified by their date (first 8 digits) - if this is different edit function accordingly
    """
    if lane == "L0":
        the_path = [x for x in paths_dict[sample][read] if run in x]
    else:
        the_path = [x for x in paths_dict[sample][read] if lane in x and run in x]
    return(the_path)

def get_list_of_values(dict_in):
    """
    for a dictonary of lists, this function returns
    a list of all the elements of the lists in the
    dict as one list
    """
    values_list = []
    for key in dict_in.keys():
        values_list = [x for x in chain(dict_in[key], values_list)]
    return(values_list)

def get_renamed_input_files(paths_dict):
    """
    get list of all paths
    this function may need to be adapted if raw file paths change
    """
    renamed_paths_dict = {}
    for sample in paths_dict.keys():
        sample_paths = []
        renamed_sample_paths = []
        renamed_paths_dict[sample] = []
        for a_path in paths_dict[sample]["R1"]:
            sample_paths.append(a_path)
        for a_path in paths_dict[sample]["R2"]:
            sample_paths.append(a_path)
        for a_path in sample_paths:
            split_path = [x for x in a_path.split("/") if "_" in x and x not in ["NGS_data", "general_data"]]
            run_name = split_path[0].split("_")[0]
            file_name = split_path[-1].strip(".fastq.gz")
            lane = "L0"
            # this will be an issue if anything else starts with "_L" but so far not the case
            if "_L" in file_name:
                lane = [x for x in file_name.split("_") if x.startswith("L") and len(x) == 2][0]
                read = [x for x in file_name.split("_") if x.startswith("R") and len(x) == 2][0]
                sample = file_name.split("_L")[0]
            else:
                read = [x for x in file_name.split("_") if x.startswith("R") and len(x) == 2][0]
                sample = file_name.split("_R")[0]
            renamed_path = f"{sample}_{lane}_{read}_{run_name}.fastq.gz"
            renamed_sample_paths.append(renamed_path)
        renamed_paths_dict[sample] = renamed_sample_paths
    return(renamed_paths_dict)

def get_all_mags():
    with open(metadata, "r") as f:
        header = f.readline()
        id_ind = header.split("\t").index("ID")
        all_mags = [line.split("\t")[id_ind] for line in f.readlines() if "unbinned" not in line.split("\t")[id_ind]]
    return all_mags

def get_rep_mags(metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        id_ind = header.split("\t").index("ID")
        reference_ind = header.split("\t").index("reference")
        rep_mags = [line.split("\t")[id_ind] for line in f.readlines() if str(line.split("\t")[reference_ind]) == "1"]
    return rep_mags

def get_medium_mags(metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        id_ind = header.split("\t").index("ID")
        quality_ind = header.split("\t").index("Quality")
        all_mags = [line.split("\t")[id_ind] for line in f.readlines() if str(line.split("\t")[quality_ind]) in ["medium", "high"]]
    return all_mags

def get_mags_of_genus(genus, metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        id_ind = header.split("\t").index("ID")
        genus_ind = header.split("\t").index("Genus")
        mags = [line.split("\t")[id_ind] for line in f.readlines() if str(line.split("\t")[genus_ind]) == genus_name]
    return mags

def get_genus_of_mag(mag, metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        id_ind = header.split("\t").index("ID")
        genus_ind = header.split("\t").index("Genus")
        for line in f:
            if line.split("\t")[id_ind] == mag:
                return line.split("\t")[genus_ind]

def get_significant_genera_list(metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        genus_ind = header.split("\t").index("Genus")
        genera = [line.split("\t")[genus_ind] for line in f.readlines() if str(line.split("\t")[quality_ind]) in ["medium", "high"]]
        genera_significant = [x for x in set(genera) if genera.count(x) > 3]
    return genera_significant

# check if the genera names all make sense for the phylogenies - rename those