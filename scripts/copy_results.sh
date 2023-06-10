nas_path=/users/aprasad/nas_recherche/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/20230313_apis_species_comparison
project_path=/scratch/aprasad/20230313_apis_species_comparison
cd ${project_path}
mkdir -p ${nas_path}/results/
mkdir -p ${nas_path}/results/00_rawreads
mkdir -p ${nas_path}/results/00_rawreads/fastqc
rsync -aviP results/00_rawreads/fastqc ${nas_path}/results/00_rawreads/
mkdir -p ${nas_path}/results/00_trimmedreads
mkdir -p ${nas_path}/results/01_cleanreads
mkdir -p ${nas_path}/results/01_trimmedconcatreads
mkdir -p ${nas_path}/results/02_motus_profile
mkdir -p ${nas_path}/results/03_host_mapping
rsync -aviP results/03_host_mapping/*.t* ${nas_path}/results/03_host_mapping
mkdir -p ${nas_path}/results/05_assembly
mkdir -p ${nas_path}/results/06_metagenomicORFs
mkdir -p ${nas_path}/results/07_MAG_binng_QC
mkdir -p ${nas_path}/results/08_gene_content
mkdir -p ${nas_path}/results/08_gene_content/00_cdhit_clustering
mkdir -p ${nas_path}/results/08_gene_content/01_profiling
mkdir -p ${nas_path}/results/08_gene_content/02_DRAM_annotations
mkdir -p ${nas_path}/results/08_gene_content/03_gene_counts/
mkdir -p ${nas_path}/results/09_MAGs_collection
mkdir -p ${nas_path}/results/11_phylogenies

rsync -aviP ${nas_path}/results/