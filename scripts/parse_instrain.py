
import os
import glob
import pandas as pd

header = False
for file in glob.glob("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/*"):
    with open(file, "r") as file_fh:
        if file.endswith("SNVs.tsv"):
            print(f"{file}")
            for line in file_fh:
                line = line.strip()
                if not header:
                    header = line
                else:
                    print(header)
                    print(line)
                    break
snvs_df = pd.read_csv("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/C1.4_profile.IS_SNVs.tsv", sep = "\t")
num_SNVs = {}
for i, gene in enumerate(snvs_df["gene"]):
    alleles = snvs_df["allele_count"].iloc[i]
    class_name = snvs_df["class"].iloc[i]
    if "SNV" in class_name:
        if gene in num_SNVs.keys():
            num_SNVs[gene] += 1
        else:
            num_SNVs[gene] = 1

genes_df = pd.read_csv("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/C1.4_profile.IS_SNVs.tsv", sep = "\t")

glob.glob("/scratch/aprasad/211018_Medgenome_india_samples/10_instrain/C1.4_profile.IS/output/*")