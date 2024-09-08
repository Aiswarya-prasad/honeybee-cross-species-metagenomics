# make sure not to run on front end node!
# ran on 2023-12-20 on Sinteractive
source ~/.bashrc
conda activate macse-env

# align nucleotides
root="...<project_dir_path>.../20230313_apis_species_comparison"
cd $root
input_dir="results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences"
output_dir="results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned"
output_dir_n="results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/nuc"
mkdir -p $output_dir_n
output_dir_a="results/11_phylogenies/05_MAG_bac120_nucleotide_trees/bac120_sequences_aligned/aa"
mkdir -p $output_dir_a

# default genetic code (for bacteria) in prodigal
for file in ${input_dir}/*.fa
do
    marker=$(basename ${file} .fa)
    if [ -f ${output_dir_n}/${marker}_aligned.aln ]; then
        echo "Skipping ${marker}"
        continue
    fi
    echo ${marker}
    macse -prog alignSequences \
          -seq ${file} \
          -gc_def 11 \
          -out_NT ${output_dir_n}/${marker}_aligned.aln \
          -out_AA ${output_dir_a}/${marker}_aligned.aln | tee ${output_dir}/${marker}_aligned.log
done