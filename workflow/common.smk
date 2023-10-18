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
        header = header.strip()
        id_ind = header.split("\t").index("ID")
        all_mags = [line.split("\t")[id_ind] for line in f.readlines() if "unbinned" not in line.strip().split("\t")[id_ind]]
    return all_mags

def get_rep_mags(metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        header = header.strip()
        id_ind = header.split("\t").index("ID")
        reference_ind = header.split("\t").index("Representative")
        rep_mags = [line.split("\t")[id_ind] for line in f.readlines() if str(line.strip().split("\t")[reference_ind]) == "1"]
    return rep_mags

def get_medium_mags(metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        header = header.strip()
        id_ind = header.split("\t").index("ID")
        quality_ind = header.split("\t").index("Quality")
        all_mags = [line.split("\t")[id_ind] for line in f.readlines() if str(line.strip().split("\t")[quality_ind]) in ["medium", "high"]]
    return all_mags

def get_mags_of_genus(genus, metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        header = header.strip()
        id_ind = header.split("\t").index("ID")
        genus_ind = header.split("\t").index("Genus")
        mags = [line.strip().split("\t")[id_ind] for line in f.readlines() if str(line.strip().split("\t")[genus_ind]) == genus]
    return mags

def get_genus_of_mag(mag, metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        header = header.strip()
        id_ind = header.split("\t").index("ID")
        genus_ind = header.split("\t").index("Genus")
        for line in f:
            if line.split("\t")[id_ind] == mag:
                return line.split("\t")[genus_ind]

def get_significant_genera(metadata):
    """
    reads the metadata file and returns a list of genera that are present in at least 3 medium or high quality MAGs
    """
    with open(metadata, "r") as f:
        header = f.readline()
        header = header.strip()
        genus_ind = header.split("\t").index("Genus")
        quality_ind = header.split("\t").index("Quality")
        genera = [line.split("\t")[genus_ind] for line in f.readlines() if str(line.split("\t")[quality_ind]) in ["medium", "high"]]
        genera_significant = [x for x in set(genera) if genera.count(x) > 3]
        if "g__" in genera_significant:
            genera_significant.remove("g__")
    return genera_significant

# check if the genera names all make sense for the phylogenies - rename those

def four_digit(n):
    n = int(n)
    if n < 10:
        return f"000{n}"
    elif n < 100:
        return f"00{n}"
    elif n < 1000:
        return f"0{n}"
    else:
        if n > 9999:
            print("Gene number seems too high. Check your genome!")
        return f"{n}"

def get_mags_for_genus_phylogeny(genus, metadata):
    with open(metadata, "r") as f:
        header = f.readline()
        header = header.strip()
        id_ind = header.split("\t").index("ID")
        genus_ind = header.split("\t").index("Phylogeny_group")
        mags = [line.strip().split("\t")[id_ind] for line in f.readlines() if str(line.strip().split("\t")[genus_ind]) == genus]
    return mags

def get_genome_path(genome, metadata, type):
    """
    gets a genome name and the mag metadata file as input
    if the genome is in the mag metadata file, it returns the path to the fna and faa files
    if the genome is not in the mag metadata file, it downloads the genome from NCBI
    and returns the fna path so that it can be annotated by prodigal inside the rule
    returns a dictionary with keys path and type
        path is the path to the fna or faa file 
        but if asked for faa but it was a downloaded so the path is to the fna, it will
        print a message and return the fna path 
        and set type to fna else the type will be what was asked for
    """
    with open(metadata, "r") as f:
        header = f.readline()
        ID_ind = header.split("\t").index("ID")
        IDs = [line.split("\t")[ID_ind] for line in f.readlines()]
        if genome not in IDs:
            print(f"{genome} not found in {metadata}, will be downloaded from NCBI")
            fna_path = f"results/11_phylogenies/downloads/{genome}.fna"
            if os.path.exists(fna_path):
                print(f"{genome}.fna already exists, skipping download")
                returned_path = fna_path
                type = "fna"
            else:
                ftp_link = ''
                with open("results/11_phylogenies/assembly_summary.txt", 'r') as f:
                    for line in f:
                        if line.startswith('#') and 'ftp_path' in line:
                            ind_link = line.strip().split('\t').index('ftp_path')
                        else:
                            line_split = line.strip().split('\t')
                            if line_split[0] == genome:
                                ftp_link = line_split[ind_link]
                if ftp_link == '':
                    print(f"ftp link for {genome} not found in assembly_summary.txt, please check {ftp_link}")
                    sys.exit(1)
                ftp_id=''.join(ftp_link.split("/")[-1])
                if os.path.exists("results/11_phylogenies/downloads"):
                    pass
                else:
                    shell("mkdir -p results/11_phylogenies/downloads")
                # faa_link = f'{ftp_link}/{ftp_id}_protein.faa.gz'
                # faa_path = f"results/11_phylogenies/downloads/{genome}.faa"
                fna_link = f'{ftp_link}/{ftp_id}_genomic.fna.gz'
                try:
                    print(f"for {genome}, fna path will be returned")
                    print(f"downloading {fna_link} to {fna_path}")
                    shell(f"wget -c {fna_link} -O {fna_path}.gz")
                    print(f"gunzipping {fna_path}.gz")
                    shell(f"gunzip {fna_path}.gz")
                    print(f"returning {fna_path} for genome. Annotate for faa file")
                    returned_path = fna_path
                    type = "fna"
                    print(f"type set to {type}")
                except:
                    print(f"could not downolad from NCBI, {fna_path} ")
                    print(f"please check the {genome} name and the assembly_summary.txt file")
                    sys.exit(1)
        else:
            fna_path = f"results/09_MAGs_collection/MAGs/{genome}.fa"
            # faa_path = f"results/09_MAGs_collection/prodigal_output/from_checkm/{genome}.faa"
            faa_path = f"results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{genome}/{genome}.faa",
            ffn_path = f"results/09_MAGs_collection/prodigal_output/renamed_for_pangenome/{genome}/{genome}.ffn",
            if type == "fna":
                returned_path = fna_path
            if type == "faa":
                returned_path = faa_path
            if type == "ffn":
                returned_path = ffn_path
        genome_path_dict = {"path": returned_path, "type": type}
        return genome_path_dict

def get_species_from_rep_mag(mag):
    if os.path.isfile('config/Species_MAG_Cluster-names.txt'):
        pass
    else:
        print("Species_MAG_Cluster-names.txt not found in config folder")
        sys.exit(1)
    with open('config/Species_MAG_Cluster-names.txt', 'r') as f:
    # MAG_species_name_final	MAG_species_Name_in_analysis	cluster	RepMAG  MAG_species_name_final_nospace
        header = f.readline()
        header = header.strip()
        mag_ind = header.split("\t").index("RepMAG")
        species_ind = header.split("\t").index("MAG_species_name_final_nospace")
        for line in f:
            if line.split("\t")[mag_ind] == mag:
                return line.split("\t")[species_ind]