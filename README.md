[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13732977.svg)](https://doi.org/10.5281/zenodo.13732977)

This pipeline requires prior set-up and the installation of appropriate tools and packages to be installed (typically using conda or mamba which is much faster). If you would like to use it or its parts and have questions or require clarification, contact Aiswarya Prasad. To explore the data or replot the figures, intermediate files and a .RData file containing R objects including functions and dataframes can be found in the [Zenodo](https://doi.org/10.5281/zenodo.13732978) repository. The code is provided for transparency and reproducibility of the analysis. The raw data files are publicly available at NCBI SRA ([PRJNA1157353](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1157353/)).

# Honeybee Cross-Species Metagenomics

The code referred to below can be found in this [GitHub repository](https://github.com/Aiswarya-prasad/honeybee-MAGs) The aim of this repository is to document the steps used to process raw shotgun metagenomic data (R1 and R2 fastq reads), assemble and bin scaffolds into MAGs, cluster them into magOTUs, estimate the coverage of each of the magOTUs and then finally also SNP profiling across samples based on the high quality MAGs chosen. Downstream analysis is performed using independent scripts and documented in their respective directories.

The pipeline outlined here was run using snakemake v7.28.1 and parallelized on the UNIL computing cluster, Curnagl (kernel: Linux 4.18.0-477.64.1.el8_8.x86_64). Information about the packages and versions used for each part of the pipeline are found in the environment specification files and the detailed list of packages in the conda environments used for downstream processing and plotting are found in the `config/envs` directory described in _Config_ section below.

## Overview

The code used can be found in the `scripts` and `workflow` directory. Most of the bioinformatic processing was done using the Snakemake rules found in the files of the `workflow` directory. Outputs from the workflow were processed using python and R scripts found in the `scripts` directory. Figures were made using the R and Rmd scripts found within `scripts/visualization`.

The `config` directory contains files including metadata and parameter specifications used in various scripts. It also contains the `config/envs` directory with YAML files specifying the packages used as conda environment specifications.

The `data` directory contains host genomes, database information etc. used by various scripts. All outputs are written among the various subdirectories of `results/`.

All paths are specified with respect to the working directory unless absolute paths are necessary (sometimes within some scripts).

Workflow was parallelized at the computation cluster at UNIL with the setup described here: [how to setup Snakemake profile](https://github.com/RomainFeron/snakemake-slurm)

Refer to the `targets` rule in the Snakefile to get the list of the outputs generated using the workflow. Further processing was done by the scripts found in `scripts/`. Refer to the respective scripts of interest for further information about the parameters used. The final summaries of outputs used for Figures are loaded into R by `scripts/visualization/Load_data.Rmd` where the paths to all the files of interest can be found.

## Raw data

Raw data for this project was sequenced in two different facilities. For 50 samples from India, all collection, analysis and sequencing was done at the collaborating lab in India and sequenced at Medgenome, India. The other 150 samples collected from Malaysia were sequenced at the GTF facility at UNIL. The files can be accessed through NCBI SRA (PRJNA1157353). Sample names reflect the species (first letter) followed by colony number and then individual number. For example, A1-1 is the first individual from the first colony of A. andreniformis and C2-3 is the third individual from the second colony of A. cerana. The files name then includes R1 or R2 indicating the forward or reverse read set and for some samples, there are multiple files for each read set from different runs and lanes as indicated in the file name. All further files (including MAGs) and figures can be created using the raw reads and the snakemake pipeline. However, for most purposes of exploration this may not be feasible. Hence MAGs, trees and some other intermediate files are included in the [zenodo](https://zenodo.org/doi/10.5281/zenodo.13732977) repository.

## Code Overview

### Snakemake workflow

* `workflow/Snakefile` - This Snakefile orchestrates the honeybee-MAGs pipeline, handling QC, assembly, binning, and downstream steps by defining the target files.
* `workflow/common.smk` - This file contains utility functions for the Snakemake pipeline, handling conversions, file path manipulations, and metadata processing.
* `workflow/trim-qc.smk` - This Snakemake file handles trimming and quality control (QC) of sequencing reads before and after trimming. It includes rules for concatenating trimmed reads and mapping reads to host and MAGs.
* `workflow/motus-profiling.smk` - This Snakemake file runs mOTUs profiling on cleaned reads to get an initial picture of the community profile.
* `workflow/assemble-qc.smk` - This Snakemake file orchestrates the assembly of metagenomes, mapping of reads, and various quality control (QC) steps. It includes rules for normalizing read coverage, mapping reads to the host genome, extracting non-host reads, and summarizing assembly statistics. Additionally, it performs taxonomic classification of contigs using Kaiju and Kraken2.
* `workflow/annotate_profile_orfs.smk` - This Snakemake file annotates and profiles open reading frames (ORFs) predicted from the assembled scaffolds. It includes rules for predicting ORFs, annotating them using various databases, and profiling their abundance across samples. The workflow ensures comprehensive functional annotation and quantification of ORFs, facilitating downstream analysis of metagenomic data.
* `workflow/backmapping-binning.smk` - This Snakemake file handles backmapping of reads to assemblies and binning of metagenomic contigs into metagenome-assembled genomes (MAGs). It includes rules for backmapping reads to assemblies, binning contigs using MetaBAT2, and summarizing binning statistics.
* `workflow/binning_summary_annotation.smk` - This Snakemake file summarizes binning results, annotates MAGs, and evaluates their quality. It includes rules for summarizing binning statistics, evaluating MAG quality using CheckM, and annotating MAGs using GTDB-Tk and DRAM.
* `workflow/mag_db_instrain.smk` - This Snakemake file profiles MAGs using InStrain to identify strain-level variation and track MAGs across samples. It includes rules for preparing MAG databases, mapping reads to MAGs, and profiling strain-level variation using InStrain.
* `workflow/mag_phylogenies.smk` - This Snakemake file infers phylogenetic trees of MAGs and external genomes to understand their evolutionary relationships. It includes rules for preparing genomes, running OrthoFinder, inferring orthologous groups, aligning orthologs, and building phylogenetic trees.

The rulegraph below summarized some of the rules including the ones involved in the steps after assembly up to taxonomic profiling.

![rulegraph](rulegraph.pdf)

### Scripts

* `scripts/` - This directory contains scripts for use in the Snakemake workflow and for processing and analyzing data generated by the Snakemake pipeline.
* `scripts/visualization/` - This subdirectory contains R scripts others for visualizing data, generating figures, and summarizing results. It includes scripts for loading data, creating plots, and generating tables for the final report and some python scripts for parsing outputs and summarizing results.

### Config

* `config/` - This directory contains configuration files used by the Snakemake pipeline and scripts.
* `config/envs` - This directory contains conda environment specifications used by the Snakemake pipeline build using mamba 1.4.2 and conda 23.3.1
* `config/config.yaml` - This file contains configuration parameters used by the Snakemake pipeline, including file paths, database locations, and metadata information.
* `config/NexteraPE-PE.fa` and `config/Adapters-PE.fa` - These files contains adapter sequences used for trimming reads in the Snakemake pipeline.
* `config/index_table.csv` - This file contains information about sample indices used in library preparation, which is used to generate adapter sequences for trimming reads.
* `config/Metadata_211018_Medgenome_india_samples.csv` - This file contains metadata information for samples sequenced in India, including sample names, species, and other details.
* `config/IsolateGenomeInfo.csv` - This file contains information about isolate genomes used in the phylogenetic analysis, including genome names, accession numbers, and other details.
* `config/Species_MAG_Cluster-names.txt` - This file contains information about species, MAGs, and clusters as defined for downstream analysis after the MAG annotation step outputs were used to determine final species names to be used in the Snakemake pipeline.

### Data

* `data/` - This directory contains data files used by the Snakemake pipeline and scripts.
* `data/host_database/` - This subdirectory contains host genomes and other reference databases used in the Snakemake pipeline.


#### Host Database

The host database is named **host_database/apis_bees_db.fasta**.

A [paper](https://academic.oup.com/gbe/article/12/1/3677/5682415) published in Dec. 2019 a high quality [_Apis dorsata_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCA_009792835.1/) as an improvement over a previous submission in 2013. The paper also mentioned studies that had previously sequenced the [_Apis florea_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCA_000184785.2) in 2012, [_Apis cerana_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_001442555.1) in 2015 (other assemblies submitted found [here](https://www.ncbi.nlm.nih.gov/assembly/organism/7460/latest/)) and [_Apis mellifera_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_003254395.2) in 2018 (other assemblies submitted listed here). So far I have not found any whole genome assemblies of _Apis andreniformis_.

These assemblies were downloaded and concatenated to make the **4_host_db**. It contains,

+ `>apis_mellifera_2018 PRJNA471592 version Amel_Hac3.1`
+ `>Apis_cerana  PRJNA235974`
+ `>Apis_cerana_mitochondrion PRJNA235974`
+ `>Apis_florea PRJNA45871`
+ `>Apis_dorsata PRJNA174631`


## Some useful intermediate files

Several large files and directories containing results are referenced in the code but not included in the repository. The repository includes the code and scripts used to process the data and generate the results. The code is provided for transparency and reproducibility of the analysis. The raw data files are publicly available and some useful intermediate files for exploration and can be found in the accompanying [Zenodo](https://doi.org/10.5281/zenodo.13732978) repository.

+ Assembled scaffold are in `results/05_assembly/all_reads_assemblies/{sample}_scaffolds.fasta` and `results/05_assembly/contig_fates/` contains the output of whokaryote, Kaiju etc. on the scaffolds
+ The directory `results/06_metagenomicORFs` contains the output of prodigal gene prediction on the scaffolds and filtered ORFs using the script `scripts/filt_orfs.py` which in turn uses files in `results/05_assembly/contig_fates/`
+ All MAGs are in `results/09_MAGs_collection/MAGs/{mag}.fa` with other directories within `results/09_MAGs_collection/` containing various checkM, GTDB annotation and dRep depreplication results
  - `results/09_MAGs_collection/drep_output` contains the output of drep dereplication
  - `results/09_MAGs_collection/gtdb_output` contains the output of GTDB annotation
  - `results/09_MAGs_collection/MAGs`contains all the MAGs
  - `results/09_MAGs_collection/dram_distill` contains the output of DRAM distill
  - `results/09_MAGs_collection/dram_output` contains the output of DRAM including a genes file (.faa)
  - `results/09_MAGs_collection/functions_list` contains a parsed list of all KOs identified by DRAM listed on column 2 with the species name appended to MAG ID as column 1
  - `results/09_MAGs_collection/MAGs_metadata_summary.tsv` contains the metadata of all MAGs collected in one file and used in downstream analysis
+ The `results/figures` directory contains all the figures many of which were further edited and arranged to make final figures for the report and tables for summaries and further processing
  - The results of filtering and annotation of ORFs combined with gene coverage information obtained by mapping clean reads back to the assemblies can be found in `results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_info.csv` and `results/figures/08-summarize_functions/gene_info_tables/{sample}_df_detected_genes_cov.csv` by the scripts `scripts/visualization/make_gene_info_summaries.py` and `scripts/visualization/summarize_gene_functions.py` respctively
  - Tip: If exploring this directory, search for the file name in the GitHub repository to find the scripts that processed and generated the files.
+ Phylogenies are in `results/11_phylogenies/03_iqtree_trees` and its neighbouring and subdirectories
  - The directory 11_phylogenies contains `results/11_phylogenies/00_genomes` with nucleotide (.ffn) and aminoacid (.faa) sequences of genes from all the MAGs and external genomes previously published (mostly used as outgroups) as listed in `results/11_phylogenies/phylo_genomes_metadata.tsv`
  - `results/11_phylogenies/01_orthofinder_input` contains the files for genomes corresponding to each genus collected together and their OrthoFinder results are in `results/11_phylogenies/02_orthofinder_results`
  - The directory `results/11_phylogenies/03_iqtree_trees` contains concate aminoacid sequence trees of core genes for each genus
  - `results/11_phylogenies/04_MAGs_gtdb` contains a tree of all MAGs made with the aminoacid sequences of the 120 bacterial genes identified and aligned by GTDB-tk
  - `results/11_phylogenies/05_MAG_bac120_nucleotide_trees`contain the trees made with the nucleotide sequences of the 120 bacterial genes identified by GTDB-tk collected and aligned using a codon-aware aligner (macse)


Below is a collection of manuscript figures and the corresponding plots that were used in making them (searching for this path / filename in the GitHub repository will provide to the script that generated the plot):

* Fig. 1 (Relevant scripts are `scripts/visualization/09-plot_community_composition.Rmd` and `scripts/visualization/06-plot_and_summarize_MAGs.Rmd` and )
  - B: `results/figures/06-figures/06-mag_by_host_ref_vs_nonref_monochrome.pdf`
  - C: `results/figures/09-figures/09-genus_prevalence_heatmap.pdf`
  - D: `results/figures/09-figures/09-Species_per_sample_boxplot.pdf`
  - E: `results/figures/00-figures/00b-qPCR_copy_numbers.pdf`
  - F: `results/figures/10-figures/SNV_mean_perc_per_host.pdf`
  - G: `results/figures/09-figures/09-PCoA_sorensen_microbiome_turnover.pdf`
  - H: `results/figures/09-figures/09-Microbe_host_dissimilarities_sorenson_turnover.pdf`
* Fig. 2 (Relevant scripts are `scripts/visualization/09-plot_community_composition.Rmd` and` scripts/visualization/10-plot_instrain_results.Rmd`)
  - A: `results/figures/09-figures/3-magotu_rohdes_scatter_pie_radius.pdf` and `results/figures/09-figures/3-magotu_rohdes_by_magotu_across.pdf`
  - B: `results/figures/10-figures/barplot_strain_level_p_value_by_genus.pdf`
  - C: `"results/figures/10-strain_nmds_IN_MY/", spec,"_nmds.pdf"` # where spec is the species name
* Fig. 3 (Relevant script is `scripts/visualization/12-phylogenies.py` and the phylogeny rules of the Snakemake workflow)
  - trees may be found in `results/11_phylogenies/03_iqtree_trees/\*/\*.treefile`
  - codiversity tests were using the code in `scripts/visualization/12-phylogenies.py`
* Fig. 4 (Relevant script `scripts/visualization/11-gene_content.Rmd`)
  - A: `results/figures/08-gene_content_plots/Number_of_kos_per_sample.pdf`
  - B: `results/figures/08-gene_content_plots/KO_pcoa_host_species_norm_aitchison.pdf`
  - C: `results/figures/11-figures/tax_func_plots/all_kos_by_host_spec_circle_marked.pdf`
  - D: `results/figures/08-gene_content_plots/cazyme_figure/CAZyme_dots.pdf` and `results/figures/08-gene_content_plots/cazyme_figure/CAZyme_tiles.pdf`
  - E: `results/figures/08-gene_content_plots/cazyme_figure/CAZyme_bars.pdf`
  - F: `"results/figures/08-gene_content_plots/gene_maps/Dysgonomonas_cazyme_genes_pathway_", sample, ".pdf"` # where sample is the sample name

For other plots not used in figures, since the figures in results/figures were moved into subfolders within after they were made by the code, the script might not match the path exactly. However, the name of the subdirectory that houses it is named with the number corresponsing to the script that generated it. For example, the script `scripts/visualization/09-plot_community_composition.Rmd` generated the plots in `results/figures` and the plots were moved to the subdirectory `results/figures/09-figures/` and so on.

* Fig. S1:
  - Map created using the scripts in `scripts/R_Shiny_map/`
* Fig. S2:
  - `06-mags_filt_contam_vs_completeness_w_shape.pdf`
  - `06-Genome_sizes_by_genus.pdf`
* Fig. S3:
  - `"results/figures/06-ANI_heatmaps/06-", genus_iter, "_ANI_heatmap_values.pdf"`
* Fig. S4:
  - `09-Cumulative_curve_by_location.pdf`
  - `09-Cumulative_curve.pdf`
* Fig. S5:
  - `results/figures/3-magotu_shared.pdf`
* Fig. S6:
  - `results/figures/10-figures/SNV_perc_all_species.pdf`
* Fig. S7:
  - `results/figures/09-relative_abundance_sep.pdf` and `results/figures/09-relative_abundance_sep_legend.pdf`
* Fig. S8:
  - `results/figures/Species_Relative_abundance/*.pdf`
* Fig. S9:
  - `09-PCoA_sorensen_microbiome.pdf`
  - `09-PCoA_sorensen_microbiome_w_colony_and_country.pdf`
  - `09-PCoA_sorensen_microbiome_nestedness.pdf`
  - `09-PCoA_sorensen_microbiome_turnover.pdf`
  - `09-PCoA_sorensen_microbiome_location_noAndreMelli.pdf`
  - `09-PCoA_sorensen_microbiome_location.pdf`
  - `09-PCoA_jaccard_microbiome_turnover.pdf`
  - `09-PCoA_jaccard_microbiome.pdf`
  - `09-PCoA_jaccard_microbiome_nestedness.pdf`
  - `09-PCoA_atchinson_relative_ab_microbiome.pdf`
  - `09-PCoA_atchinson_microbiome.pdf`
* Fig. S10:
  - `08-gene_content_plots/KO_*.pdf`
* Fig. S11:
  - `results/figures/11-figures/tax_func_plots/*.pdf`
* Fig. S12:
  - `08-gene_content_plots/cazyme_*.pdf`
* Fig. S13:
  - `results/figures/01-figures/01-Sequencing_depth_summary_IN_MY`

For any dataframes or files not loaded in the script of interest, the data is likely loaded in the script `scripts/visualization/Load_data.Rmd` which loads all the data used in the figures and tables of the manuscript. The saved `results/figures/workspace_generaldata_chunks_20230611.RData` should contain all the important dataframes used in the manuscript and allow the repetition of most of the code in the other scripts. The environment as described in `config/envs/rmd-env.yaml` was used for R and Rmd code and the one in `config/envs/scripts-env.yaml` was used for python code.