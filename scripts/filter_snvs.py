#!/usr/bin/env python3
import os

input_vcf = os.path.join(os.getcwd(), snakemake.input["vcf"])
output = os.path.join(os.getcwd(), snakemake.output["freq"])

summary = os.path.join(os.getcwd(), snakemake.params["summary"])

# Step 1: Using the sub-routine "extract_dp_data", the relative coverage of a candidate SNV is calculated for each sample. The coverage is set to zero if it less than 0.1 (10%). Additionally, the fraction of samples with missing data (due to low coverage) is calculated. If the missing fraction is larger than 0.1 (10%), the SNV is skipped

# Step 2: Using the sub-routine "is_polymorphic", the script will check whether the candidate SNV is still polymorphic after rare SNVs are removed. If all samples have the values 0 or 1 (i.e. no polymorphism, or only polymorphic relative to the reference genome), the SNV is skipped

# Step 3: If an SNV passes the first two steps, it is printed to the output file

# The output file is named [SDP]"_filt.freq". In the first column, the reference genome, SNV positiion and base change relative to the reference is indicated. The remaining columns display the relative proportion of the SNV within each sample. In cases where there is more than one two alleles for the same position, multiple lines are generated.


def extract_info(data):
    info = {"missing": 0, "data_AO": []}
    if "\n" in data:
        data.remove("\n")
    num_tabs = len(data)
    missing_data_count = 0
    data_AO = []
    for tab in data:
        try:
            AO = tab.split(":")[5]
        except:
            AO = "."
        if "." in tab and AO == ".":
            missing_data_count += 1
            data_AO.append(-1)
        else:
            AO = tab.split(":")[5]
            DP = tab.split(":")[1]
            allele_relcov = int(AO)/int(DP)
            if allele_relcov <= 0.1:
                data_AO.append(0)
            else:
                data_AO.append(allele_relcov)
    missing_fraction = missing_data_count/num_tabs
    info["missing"] = missing_fraction
    info["data_AO"] = data_AO
    return(info)

def is_polymorphic(all_data):
    data = []
    for tab in all_data:
        if int(tab) == -1:
            pass
        else:
            data.append(tab)
    nb_data = len(data)
    sum_data = sum(data)
    result = True
    if sum_data == nb_data or sum_data == 0:
        result = False
    return(result)


with open(output, "w+") as out_fh:
    with open(input_vcf, "r") as vcf_fh:
        num_polymorphic_sites = 0
        num_skipped = 0
        for line in vcf_fh:
            if "##" in line:
                continue
            if line.split("\t")[0] == "#CHROM":
                samples = line.split("\t")[9:]
                out_fh.write("SNP_info")
                for sample in samples:
                    out_fh.write(f"\t{sample}")
                out_fh.write("\n")
            else:
                data = line.split("\t")[9:]
                data_info = extract_info(data)
                missing_fraction = data_info["missing"]
                data_AO = data_info["data_AO"]
                if missing_fraction >= 0.1:
                    num_skipped += 1
                    continue
                if is_polymorphic(data_AO):
                    num_polymorphic_sites += 1
                else:
                    continue
                genome_id = line.split("\t")[0]
                pos = line.split("\t")[1]
                ref_base = line.split("\t")[3]
                alt_base = line.split("\t")[4]
                snv_info_line = genome_id+":"+pos+":"+ref_base+">"+alt_base
                out_string = snv_info_line
                for num in data_AO:
                    out_string = out_string + "\t" + str(num)
                out_fh.write(f"{out_string}\n")

summary_has_header = False

if not os.path.isfile(summary):
    summary_fh = open(summary, "w")
    summary_fh.write("Input\tPolymorphic\tMissing\n")
    summary_has_header = True
    summary_fh.close()
else:
    summary_fh = open(summary, "r")
    for line in summary_fh:
        if line.startswith("Sample"):
            summary_has_header = True
    summary_fh.close()

with open(summary, "a") as summary_fh:
    if not summary_has_header:
        summary_fh.write("Input\tPolymorphic\tMissing\n")
    input_name = input_vcf.split("/")[-1].split(".")[0]
    summary_fh.write(f"{input_name}\t{num_polymorphic_sites}\t{num_skipped}\n")
