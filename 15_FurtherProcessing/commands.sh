# This is a command from GTDB tk to make a tree of all the MAGs
# to use for UniFrac later
# Later incorporate this in the Snakemake script
conda activate 211018_Medgenome_india_samples-gtdb-env
cd /scratch/aprasad/211018_Medgenome_india_samples/
cd 15_FurtherProcessing
export GTDBTK_DATA_PATH="/work/FAC/FBM/DMF/pengel/spirit/aprasad/gtdb/release207_v2/"
gtdbtk infer --msa_file "../06_MAG_binning/gtdbtk_out_dir/align/gtdbtk.bac120.user_msa.fasta.gz" --out_dir "gtdb_inferred_tree"


# test out motupan and compare with orthofinder results
conda activate 211018_Medgenome_india_samples-mags-env

mkdir -p motupan_test

checkm_out="motupan_test/generic_completeness.tsv"
cat ../06_MAG_binning/evaluate_bins/checkm_drep_summary.txt | tr "," "\t" | sed -e 's/.fa//' | sed -e 's/genome/Bin Id/' | sed -e 's/comp/Comp/' | sed -e 's/cont/Cont/' > ${checkm_out}
input_faa_dir="../database/MAGs_database_Orthofinder/g__Lactobacillus/"
out_dir="motupan_test/g__Lactobacillus"
mkdir -p ${out_dir}
OG_file="motupan_test/g__Lactobacillus/motupan_cog.tsv" # using scripts/preapare_motupan_OG_file.py
mOTUpan.py --gene_clusters_file ${OG_file} --boots 100 -o ${out_dir}/mOTUpan.tsv --checkm ${checkm_out} | tee ${out_dir}/mOTUpan.log

