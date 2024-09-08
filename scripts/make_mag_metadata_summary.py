import os
import sys
import argparse
import pandas as pd
import numpy as np

"""
reads the outputs of gtdb, drep and checkm and makes a summary table of the MAGs with the following columns:
ID - mag id
Domain
Phylum
Class
Order
Family
Genus
Species
magOTU
Sample
Completeness
Contamination
Size - Genome size (bp)
scaffolds - number of scaffolds
contigs - number of contigs
N50_scaffolds - (scaffolds)
N50_contigs - (contigs)
Mean_scaffold_length - Mean scaffold length (bp)
Mean_contig_length - Mean contig length (bp)
Longest_scaffold - Longest scaffold (bp)
Longest_contig - Longest contig (bp)
genes - # predicted genes
Quality - high, medium, low 
          (high - completeness > 90% and contamination < 5%
           medium - completeness > 50% and contamination < 10%
           low - completeness < 50% and contamination > 10%)

--gtdb {input.gtdb} \
--checkm {input.checkm} \
--drep_gI {input.drep_gI} \
--drep_Wi {input.drep_Wi} \
--drep_S {input.drep_S} \
--outfile {output.metadata}

"""


# add parser arguments
parser = argparse.ArgumentParser(description='make a summary table of the MAGs')
parser.add_argument('--gtdb', help='gtdb output file')
parser.add_argument('--checkm', help='checkm output file')
parser.add_argument('--drep_gI', help='drep output file with gANI')
parser.add_argument('--drep_Wi', help='drep output file with wANI')
parser.add_argument('--drep_S', help='drep output file with secondary clustering')
parser.add_argument('--outfile', help='output file name')
args = parser.parse_args()

# read in the gtdb output file

gtdb = args.gtdb
checkm = args.checkm
drep_gI = args.drep_gI
drep_Wi = args.drep_Wi
drep_S = args.drep_S
outfile = args.outfile
gtdb_ar = "".join(gtdb.split(".bac120.summary.tsv")[0]) + ".ar53.summary.tsv"

# gtdb = "/...<project_dir_path>.../20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/gtdb_output/classify/All_mags_sub.bac120.summary.tsv"
# gtdb_ar = "".join(gtdb.split(".bac120.summary.tsv")[0]) + ".ar53.summary.tsv"
# checkm = "results/09_MAGs_collection/All_mags_sub/checkm_merged.tsv"
# drep_gI = "/...<project_dir_path>.../20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/drep_results/data_tables/genomeInfo.csv"
# drep_Wi = "/...<project_dir_path>.../20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/drep_results/data_tables/Widb.csv"
# drep_S = "/...<project_dir_path>.../20230313_apis_species_comparison/results/09_MAGs_collection/All_mags_sub/drep_results/data_tables/Sdb.csv"

info_df = pd.DataFrame()
checkm_info = pd.read_csv(checkm, sep="\t", header=0)
gtdb_info = pd.read_csv(gtdb, sep="\t", header=0)
gtdb_info_ar = pd.read_csv(gtdb_ar, sep="\t", header=0)
drep_gI_info = pd.read_csv(drep_gI, sep=",", header=0)
drep_S_info = pd.read_csv(drep_S, sep=",", header=0)
drep_C_info = pd.read_csv(drep_S.split("Sdb.csv")[0] + "Cdb.csv", sep=",", header=0)
drep_W_info = pd.read_csv(drep_S.split("Sdb.csv")[0] + "Wdb.csv", sep=",", header=0)
drep_Wi_info = pd.read_csv(drep_Wi, sep=",", header=0)
# drep_M_info = pd.read_csv(drep_S.split("Sdb.csv")[0] + "Mdb.csv", sep=",", header=0)
# drep_N_info = pd.read_csv(drep_S.split("Sdb.csv")[0] + "Ndb.csv", sep=",", header=0)
# # drep_M_info["ID"] = drep_M_info["genome"].apply(lambda x: x.split(".fa")[0])
# # drep_N_info["reference"] = drep_N_info["reference"].apply(lambda x: x.split(".fa")[0])
# # drep_N_info["querry"] = drep_N_info["querry"].apply(lambda x: x.split(".fa")[0])

info_df_cols = ['ID', 
          'Completeness', 'Contamination',
          'Size', 'scaffolds', 'contigs',
          'N50_scaffolds', 'N50_contigs', 
          'Mean_scaffold_length',
          'Mean_contig_length',
          'Longest_scaffold',
          'Longest_contig',
          'genes']
checkm_df_cols = ['Bin Id', 
             'Completeness', 'Contamination',
             'Genome size (bp)', '# scaffolds', '# contigs',
             'N50 (scaffolds)', 'N50 (contigs)', 
             'Mean scaffold length (bp)',
             'Mean contig length (bp)',
             'Longest scaffold (bp)',
             'Longest contig (bp)',
             '# predicted genes']
info_df[info_df_cols] = checkm_info[checkm_df_cols]
info_df["Size_readable"] = info_df["Size"].apply(lambda x: str(round(x/1000000, 2))+"M" )
info_df["Sample"] = info_df["ID"].apply(lambda x: x.split("_")[0])
info_df["Type"] = info_df["ID"].apply(lambda x: "unbinned" if "unbinned" in x else "MAG")
info_df["Quality"] = info_df.apply(lambda x: "high" if x["Completeness"] > 90 and x["Contamination"] < 5 else "medium" if x["Completeness"] > 50 and x["Contamination"] < 5 else "low", axis=1)
classification_info = gtdb_info["classification"].str.split(";", expand=True)
classification_info["ID"] = gtdb_info["user_genome"]
classification_info.rename(columns={0:'Domain', 1:'Phylum', 2:'Class', 3:'Order', 4:'Family', 5:'Genus', 6:'Species'},inplace=True)
info_df = pd.merge(info_df, classification_info, left_on = "ID", right_on = "ID", how="left")
classification_info_ar = gtdb_info_ar["classification"].str.split(";", expand=True)
classification_info_ar["ID"] = gtdb_info_ar["user_genome"]
classification_info_ar.rename(columns={0:'Domain', 1:'Phylum', 2:'Class', 3:'Order', 4:'Family', 5:'Genus', 6:'Species'},inplace=True)

for column in classification_info_ar.columns:
    info_df[column].fillna(classification_info_ar[column], inplace=True)


drep_gI_info["ID"] = drep_gI_info["genome"].apply(lambda x: x.split(".fa")[0])
drep_S_info["ID"] = drep_S_info["genome"].apply(lambda x: x.split(".fa")[0])
drep_C_info["ID"] = drep_C_info["genome"].apply(lambda x: x.split(".fa")[0])
drep_W_info["ID"] = drep_W_info["genome"].apply(lambda x: x.split(".fa")[0])
drep_Wi_info["ID"] = drep_Wi_info["genome"].apply(lambda x: x.split(".fa")[0])

info_df["Score"] = info_df["ID"].apply(lambda x: drep_S_info[drep_S_info["ID"] == x]["score"].values[0] if x in drep_S_info["ID"].values else np.nan)
info_df["magOTU"] = info_df["ID"].apply(lambda x: drep_C_info[drep_C_info["ID"] == x]["secondary_cluster"].values[0] if x in drep_C_info["ID"].values else np.nan)
info_df["closest_cluster_member"] = info_df["ID"].apply(lambda x: drep_Wi_info[drep_Wi_info["ID"] == x]["closest_cluster_member"].values[0] if x in drep_Wi_info["ID"].values else np.nan)
info_df["furthest_cluster_member"] = info_df["ID"].apply(lambda x: drep_Wi_info[drep_Wi_info["ID"] == x]["furthest_cluster_member"].values[0] if x in drep_Wi_info["ID"].values else np.nan)

max_score_indices = info_df[info_df["Quality"].isin(["medium", "high"])].groupby("magOTU")["Score"].idxmax()
print(f"Identified {len(max_score_indices)} representative MAGs")
info_df["Representative"] = 0
info_df.loc[max_score_indices, "Representative"] = 1
info_df.to_csv(outfile, sep="\t", index=False)