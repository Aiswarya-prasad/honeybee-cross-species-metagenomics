Note that this pipeline is a work in progress so some parts of the documentation may not be up-to-date yet. If you would like to have questions or require clarification, contact Aiswarya Prasad. (aiswarya.prasad@unil.ch)

# honeybee-MAGs

The aim of this pipeline is to document the steps used to process raw shotgun metagenomic data (R1 and R2 fastq reads), assemble and bin scaffolds into MAGs, cluster them into magOTUs, estimate the "core" coverage of each of the magOTUs and then finally also SNP profiling across samples based on the high quality MAGs chosen. The pipeline also anotates ORFs and MAGs and includes a section for SDP validation (from Kirsten_Ellegaard's work) and making phylogenies of MAGs. Downstream analysis is performed using independent scripts and documented in their respective directories. All scripts are present in the scripts directory in general.

In this file, each of the rules in the Snakefile and the associated are explained in more detail. Every step is taken care of by Snakemake as instructed in the Snakefile. None of the scripts are expected to be run independently. The overall approach is to have to the ultimate rule named `backup` which will ask for all the desired outputs and copy important files (currently it copies the entire working directory but the backup script will be modified to only copy important checkpoints) to a specified local or remote backup location.

If you wish to not run any rules, remove the entries correponding to their outputs from the list of inputs to these rules.

## Getting started

All paths are specified with respect to the working directory. Snakemake also happily understands this. In the process of writing this pipeline "211018_Medgenome_india_samples" was the name of my project and the working directory, `/scratch/aprasad/211018_Medgenome_india_samples` on curnagl. So if you see this appear anywhere, you know what to do!

To get started with this pipeline you will need the following files.

+ config/config.yaml
  - This is used by Snakemake to get some important information. It is best to create it manually or using a script that will read a metadatasheet compiled manually.
    + The list of samples is used by snakemake to decide which samples are to be processed
    + The Adapters section specifies the adapters.fa file to be given to trimmmomatic for each group of samples
    + The section READS and TYPES are used to infer the names/paths of certain files
    + GENOME_DBs specifies the path to various databases if they are used in any rules
    + ProjectPath is the absolute path to the project working directory and may be used sometimes as a prefix to provide the absolute path in some places.
    + BackupPath is the string describing the path to where the backups should be written and if it is remote, the <username>@xx.x.xx.xx: part as well
    + LocalBackup is the variable that will used to decide whether to submit the backup job to the cluster or run it locally. Eg. if you want backup on the /nas partition on curnagl, it can only be accessed from login nodes and if you do not run it locally, the job will fail as /nas will not be found.
+ config/NexteraPE-PE.fa
  - This is the file to provide to trimmomatic for trimming
+ config/index_table.csv (optional)
  - In some cases the library prep may hae involved complicated combinations of indices and the adapters (eg. Kapa HyperPrep kit) to be trimmed will have to be inferred from this information. There used to be a script for this but it was lost. So it will have to be re-written at some point. config/Adapters-PE.fa is the output made by the script that parses the index table. For now, since the script does not exist, Adapters.fa is just kept as initial input an the rule to make it is commented out.
+ config/Adapters-PE.fa (not needed if to be created from index_table.csv)
+ config/Metadata_211018_Medgenome_india_samples.csv
+ config/initialize_project-211018_Medgenome_india_samples.sh
  - This might eventually be expanded into a complete scripts to enable quick and pain-free project initialisation. Currently it is more or less a bunch of `rsync` commands to copy raw files and important starting point files from a remote location and into the working directory to get started.
+ config/IsolateGenomeInfo.csv
  - This file is required for the part of the pipeline performing the inference of orthologs and phylogenetic trees of the MAGs. The trees includes some reference genomes from the honeybee microbiome database and outgroups for the genuses expected to be found. This however, depends on the dataset being used and may have to be modified or skipped depending on the analysis.
  - The file contains the following fields comma-separated without spaces padding the commas
    + ID,Accession,Locus_tag_KE_DB,Strain_name,Phylotype,SDP,Species,Host,Study,Origin,Source_database,Cluster,Group,Genus,Strain_name_alt
    + ID can be the same as Strain_name. In cases where the strain on ncbi contains spaces, they are replaced by '_'. Accession is the genbank ID which can be used to find the genome online (eg. GCA_019469245.1). The locus tag is only present to relate it to the genomes (if needed) in Kirsten's database which mostly are named by IMG locus tags. The SDP is 'NA' if it is a MAG or if the SDP is unknown. Some of these columns are not used but if they are to be removed or reordered, the scripts parsing this file must be modified accordingly.
+ config/reinitialize_project-211018_Medgenome_india_samples.sh
  - A list of rsync commands to copy files back from a remote backup to resume the analyses.
+ 00_RawData/<sample>_R*.fastq.gz
  - Raw data in fastq.gz format kept in the directory called 00_RawData
+ 211018_Medgenome_india_samples_Report.Rmd
+ database/4_host_db
  - This is required for host mapping. It is a fasta file with host genomes concatenated. Currently the pipeline only assembles host-filtered reads. But if the host genomes is unavailable or you would like to skip this step, ignore these steps and change the inputs of dependent rules to take the trimmed files directly.
+ envs/<env>.yaml
  - The conda environments needed for each rule will be created based on the yaml files in this directory. So at least the files corresponding to the conda environments that your rules will use must be present.
+ scripts/<script>.xxx
  - The scripts used by various rules are all to be in this directory. The list of all scripts and theit use can be found later in this file.

## Checkpoints

Snakemake uses the file/path provided to the final rule (targets) and then looks for rules with output paths that match the pattern of wildcards in the respective target. Checkpoints are a provision of snakemake as described [here](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) to allow for rules to be defined with outputs that will depend on another rule. For example, the list of MAGs on which to perform downstream steps on will only be known after the binning step has been completed. So far, there are two checkpoint as described below.

+ `checkpoint make_phylo_table`
  - The inputs that this rule takes for are "config/IsolateGenomeInfo.csv", "06_MAG_binning/all_genomes.csv" (provided by rule extract_mag_lists), "06_MAG_binning/gtdbtk_out_dir/classify/gtdbtk.bac120.summary.tsv" (provided by rule gtdb_annotate) and "06_MAG_binning/drep_results/data_tables/Cdb.csv" (provided by rule drep).
  - It uses the script `scripts/make_phylo_table.R` and `scripts/csv_to_tsv.py` to make tsv files which will be parsed later to obtain information about which MAGs downstream processes should run on.
  - The various outputs made by this rule are as follows:
    - They all have at least the headers:
      - ID - name of the MAG or stain name (modified) of the genome.
      - Accession - Genbank ID (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Locus_tag (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Strain_name (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Phylotype (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - SDP (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Species (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Host - for MAGs, the host from which the MAG was assembled and for Isolates, the source from which it was isolated.
      - Study - (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Origin - (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Source_database (for isolate genomes only - same as in config/IsolateGenomeInfo.csv)
      - Cluster - The name of the 95% dereplicated cluster that it was assigned to by drep (refer to <rule drep>)
      - Group - This field was used for assigning isolates to appropriate genus-level groups with MAGs for the sake of the tree. This tree is filled manually for isolated genomes when making "config/IsolateGenomeInfo.csv"
      - Genus, Family, Order, Class - Taxonomic information assigned for MAGs by gtdb (see <rule gtdb_annotate>) (only for MAGs)
      - Sample - The sample from which the MAG was assembled (only for MAGs)
      - Group_auto - The Group assigned to each row automatically by the script. For isolates, it is inferred either from the information from the Group column of "config/IsolateGenomeInfo.csv" file or the phylotype. Make sure that the way this is done (based on a pre-defined dictionary written in the script) makes sense for your inputs. For the MAGs, it is the Genus assigned by gtdb.
      - Ref_status - 1 if the MAG is a representative genome for its respective drep cluster, 0 otherwise. The MAG with the highest score within a given drep cluster is designated as the representative genome of the cluster (no matter how low the score). The method in which this is currently done assumes that, upon filtering by completeness and contamination, the MAG with the highest score will only be filtered out if all of the other members of its cluster are also filtered out. This makes sense because the score is weighted in such a way that completeness and contamination are given heavy weightage.
    - "06_MAG_binning/all_GenomeInfo_auto.tsv" : This file contains the information for all Isolates and all MAGs
    - "06_MAG_binning/ForTree_GenomeInfo_auto.tsv" : This file contains the information for all Isolates and filtered MAGs (Completeness > 50, Contamination < 5)
    - "06_MAG_binning/MAGs_GenomeInfo_auto.tsv" : This file contains the information for and all MAGs
    - "06_MAG_binning/MAGs_filt_GenomeInfo_auto.tsv" : This file contains the information filtered MAGs (Completeness > 50, Contamination < 5)


## Utility functions

Some functions in the script are used to provide inputs to various rules and to parse the checkpoint outputs so the list of wildcards can be obtained.

## Rules

  + `rule targets`
    + This is the master rule that specifies all the output files / directories desired and accoringly runs the rules that provide the specified output. Most of the outputs from the different rules are used for making the report. Just the ones that are not used by the report are mentioned here.
  + `rule raw_qc`
    + Uses raw reads and produces fastqc outputs
    + This rule runs [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on raw files and saves the output in "/fastqc/raw".
  + `rule make_adapters`
    + Uses the "config/index_table.csv" and the script `scripts/write_adapters.py` (absent for now) and produces the "Adapters-PE.fa" file containing indexed adapters.
  + `rule trim`
    + This rules trims reads using [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).
    + The _Adapters-PE.fa_ files is used.
  + `rule trim_qc`
    + This rule runs [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on trimmed files and saves the output in `fastqc/trim`.
  + `rule index_bwa`
    + Indexes genomes in `database` for use by [bwa](http://bio-bwa.sourceforge.net/) using [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml#3).
  + `rule index_samtools`
    + Indexes genomes in `database` for use by [samtools](http://www.htslib.org/doc/#manual-pages).
  + `rule make_genome_list`
    + Creates a text file corresponding to each set of genomes in `database` to be used when we need to know which genomes are present in given genome database.
  + `rule run_motus`
    + Uses the trimmed reads "01_Trimmed/{sample}_R*_trim.fastq.gz" to produce an [motus](https://github.com/motu-tool/mOTUs) output.
  + `rule merge_motus`
    + merges all the motu outputs from the previous rule and creates the merged output "08_motus_profile/samples_merged.motus" to be parsed later for visualization.
  + `rule host_mapping`
    + Uses [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml#:~:text=BWA%20is%20a%20software%20package,such%20as%20the%20human%20genome.&text=BWA%2DMEM%20also%20has%20better,genome%20(the%20index%20command)) to map reads for each sample to a database containing host genomes, `database/4_host_db`.
    + Unmapped alignments identified by samtools with the option `-f4` are stored in a seperate bam file to be used later. All others will be deleted.
    + The bam file with all alignments is used to create a flagstat output which will be parsed later. It also creates a coverage histogram and a tsv file summarising coverage.
  + `rule host_mapping_extract_host_filtered_reads`
    + Reads that did not map to the host database are extracted and then mapped to the microbiome database.
    + They are extracted using [picard](https://broadinstitute.github.io/picard/).
    + The option `-Xmx8g` ensures that java is given 8 GB memory. If suffecient memory is not allocated, the job will fail. The extracted reads will be used later "02_HostMapping/{sample}_R*_host_unmapped.fastq" and also referred to as host-filtered reads.
    + This is done in order to make the assembly a little easier.
  + `rule microbiomedb_mapping`
    + The host unmapped reads extracted earlier are mapped to the microbiome database.
    + Mapped reads are extracted using a perls script as follows. First, unmapped reads are excluded using `-F4` and then supplementary reads are excluded `-F 0x800`. Finally, the remaining reads are sent through `scripts/filter_sam_aln_length.pl`. The script filters away reads that have less than 50bps matching in the alignment.
    + It keeps the [samtools](http://www.htslib.org/doc/samtools.html) [flagstat](https://www.biostars.org/p/84396/) result but removes all the others.
  + `rule microbiomedb_direct_mapping`
    + trimmed reads are directly mapped to the microbiome database to see how much of the reads map in this way as compared to after host-filtering. Some reads that are of microbial origin may map to the host database due to contamination in the host database or even non-specifically by chance. To ensure that we do not lose too many reads in this way, upon host filtering, the flagstat output from this mapping is useful.
  + `rule assemble_host_unmapped`
    + Uses the host-filtered reads to give to spades using the meta option to assemble.
    + The bash script is eritten to make sure assemblies are resumed if an older output directory exists.
    + The script `scripts/parse_spades_metagenome.py` is then used to extract only scaffolds that pass the coverage and length threshold that can be specified in the parameters section of the rule.
    + The scffolds, log and assembly graphs are copied out of the output directory and kept but the rest of the directory is deleted.
  + `rule map_to_assembly`
    + This rule maps the host-filtered reads to the assembly and creates the flagstat output used to visualise how well the assembly is represented by the reads for each sample.
  + `rule summarize_mapping_assembly`
    + It used the scaffolds and the flagstat output from mapping create a summary across all samples that can be parsed and visualized by the script making the report.
  + `rule backmapping`
    + Carries out NxN mapping of all N samples and only keeps the depth file that needs to be supplied to metabat2 for binning after the mapping for each sample (to be merged before being given to metabat2). It works on 1 sample at a time and maps the assembly of that sample against the reads from all the samples. Each iteration of the rule will result in N files. For example, when there are 50 samples in total and the rule runs on sample X, there will 50 files called "<sample>_mapped_to_X.depth". The rule only creates one bam file at a time. However, if the rule is run on 10 samples simultaneously, there will be about 10 bam files at a time.
  + `rule merge_depths`
    + All depth files of all the samples mapped to a given sample will be merged into 1 file that will be given to metabat2.
  + `rule binning`
    + For each sample, this rule provides the set of scaffolds assembled and the merged depth file with the coverages in each sample of these scaffolds to metabat2.
    + The output is a directory containing as many MAG.*.fa files as there are bins with the scaffolds distributed across them.
  + `rule process_metabat2`
    + Since the names given to bins by metabat2 is duplicated across samples, this rule renames them by including the name of the sample from which the bin was made.
  + `rule summarize_metabat2_contig_fates`
    + It creates sample-wise summaries with sample, contig_name, length, coverage, passed_filter (whether the scaffold was >1000 length and >1 coverage), binned (whether it was binned) and bin_name.
  + `rule summarize_metabat2_contig_coverages`
    + It creates sample-wise summaries with the columns bin_name , contig (or scaffold name), avg_coverage (from depth file made by jgi_summarize_bam_contig_depths for metabat2. it is actually  median not mean), variance avg_coverage (from depth file made by jgi_summarize_bam_contig_depths for metabat2), sample
    + More information about how these depths are calculated can be found [here](https://gitlab.com/robegan21/MetaBAT). Briefly (as summarised [here](https://bitbucket.org/berkeleylab/metabat/issues/48/jgi_summarize_bam_contig_depths-coverage)),
      + The edges of a contig are generally excluded from the coverage counts up to a default of 75 bases or the average read length (--includeEdgeBases, --maxEdgeBases). This is because, generally mappers have a difficult time aligning a partial read to a contig when it would extend off edge and the coverage ramps up from 0 to the true coverage in this region
      + reads that map imperfectly are excluded when the %ID of the mapping drops below a threshold (--percentIdentity=97). MetaBAT is designed to resolve strain variation and mapping reads with low %ID indicate that the read actually came from a different strain/species.
      + clips/insertions/deletions/mismatches are excluded from the coverage count -- only the read bases that exactly match the reference are counted as coverage. This generally has a small effect, except in the case of long reads from PacBio and Nanopore.
    + `rule checkm_evaluation`
      + This rule creates a checkm summary of each bin (MAG) from each sample, one sample at a time. The output contains an entry for each MAG from the sample with the columns, Bin Id, Marker lineage , # genomes, # markers, # marker sets, 0 , 1, 2, 5+, Completeness, Contamination, Strain heterogeneity
      + A detailed description of each column can be found [here](https://github.com/Ecogenomics/CheckM/wiki/Reported-Statistics#qa)
    + `rule prepare_info_for_drep`
      + drep requires a file with the columns genome, completeness and contamination. This rule compiles this information from across samples to create the file "06_MAG_binning/evaluate_bins/checkm_drep_summary.txt".
    + `rule drep`
      + Runs drep dereplicate on **all** the MAGs. It creates various outputs that are described [here](https://drep.readthedocs.io/en/latest/advanced_use.html) and will used by various rules downstream.
    + `rule extract_mag_lists`
    + `rule prep_gtdb_annotate`
    + `rule gtdb_annotate`
      + Need to make sure that the database is set up as explained [here](https://ecogenomics.github.io/GTDBTk/installing/index.html) and provide the path to the database in the params section to the rule.
    + `rule prepare_genomes`
    + `rule annotate`
    + `rule prepare_faa`
      + For phylogenies (includes specified isolate genomes)
    + `rule run_orthofinder_phylo`
      + For phylogenies (includes specified isolate genomes)
    + `rule summarise_orthogroups`
      + For phylogenies (includes specified isolate genomes)
    + `rule get_single_ortho_phylo`
      + For phylogenies (includes specified isolate genomes)
    + `rule extract_orthologs_phylo`
      + For phylogenies (includes specified isolate genomes)
    + `rule align_orthologs`
      + For phylogenies (includes specified isolate genomes)
    + `rule prune_and_concat`
      + For phylogenies (includes specified isolate genomes)
    + `rule make_tree`
      + For phylogenies (includes specified isolate genomes)
    + `rule concat_all_mags`
      + Create and index "database/MAGs_database" and index to map reads to see how much the MAGs recruit and how well they are covered
    + `rule copy_mag_database_annotation`
      + To create a compiled directory of annotation outputs for all MAGs which Kirsten's scripts expect (also it is neat if the other previous directories are then removed - To be implemented later)
    + `rule prepare_faa_mag_database`
      + Prepare faa files in groups (genus-level) for orthofinder to run.
    + `rule run_orthofinder_mag_database`
    + `rule summarise_orthogroups_mag_database`
    + `rule get_single_ortho_mag_database`
    + `rule extract_orthologs_mag_database`
    + `rule calc_perc_id_mag_database`
    + `rule filter_orthologs_mag_database`
    + `rule summarise_orthogroups_filtered_mag_database`
    + `rule make_bed_files_mag_database`
    + `rule map_to_MAGs`
    + `rule core_cov`
    + `rule merge_core_cov`
    + `rule core_cov_plots`
    + `rule make_MAG_reduced_db`
      + Read the checkpoint output and use biopython to only include rep genomes in the final database called "database/MAGs_rep_database" which only contains the representative genomes.
    + `rule prep_for_instrain`
    + `rule map_to_rep_MAGs`
    + `rule instrain_profile`
    + `rule instrain_compare`
<!-- + "Continue readme update here - need to finish function to get rep genomes and then rerun everything after backmapping!" -->
  + `rule run_orthofinder`
    + Runs [Orthofinder](https://github.com/davidemms/OrthoFinder) for each phylotype.
    + **Before** running this, group genomes by phylotype in directories for Orthofinder to be able to get which groups to consider together. When the genomes for the database are downloaded at `database/faa_files/{genome}`, they are all in one directory. Grouping was done using `scripts/rearange_faa.py`. As written, it is to be run from the scripts directory in which it resides (!! it uses relative paths !!).
    + faa files for each genome comes from the respective databese (NCBI for example)
    + When orthofinder finishes, the following file will be generated and used for the following steps, `database/faa_files/{phylotype}/OrthoFinder/Results_dir/Orthogroups/Orthogroups.txt`.
    + The file _Orthogroups.txt_ contains a list of orthogroups. Eg, each line would look like
      - **OG0000003**: C4S76_01365 C4S76_01370 C4S76_01375 C4S77_06100 C4S77_06130 C4S77_06135 C4S77_06775 C4S77_06780 C4S77_06785 C4S77_06790 C4S77_06795 C4S77_06800 C4S77_06805 C4S77_06810 C4S77_09595 C4S77_09600 C4S77_09605 C4S77_09610 **C4S77_09615 C4S77_09620 C4S77_10540 Ga0307799_111506**
      - where, **OG0000003** is an orthogroup for this group of genomes (phylotype) and **C4S77**, **Ga0307799** etc. are genomes that belong to that group. **09615, 09620, 10540** are genes from **C4S77** and **111506** from **Ga0307799** that belong to orthogroup OG0000003.
  + `rule get_single_ortho`
    + The files _Orthogroups.txt_ is parsed by `scripts/get_single_ortho.py` and single-copy orthologs are written to `database/Orthofinder/{phylotype}_single_ortho.txt`
    + The script reads each orthogroup and counts the number of genomes present in genes of that orthogroup. If the number of genes in the orthogroup and the number of genomes in the orthogroup are the same as the total number of genomes in the database for said phylotype, the genes in the group are considered single-copy core genes and included for core coverage estimation.
  + `rule extract_orthologs`
    + This rule prepares files with sequences of orthologs in order to calculate percentage identity (perc_id).
    + First, it reads the file `database/Orthofinder/{phylotype}_single_ortho.txt` and gets all the genome-ids present in the ortholog-file, and all the gene-ids associated with each gene-family. Using this list it extracts and stores the sequences of each of the genes of an orthogroup in an faa file and ffn file corresponding to each group in the directory`database/Orthofinder/{phylotype}_ortho_sequences/`.
    EXAMPLE:
    + `cat ./database/Orthofinder/firm5_ortho_sequences/OG0001034.faa` \
      **\>Ga0061073_1479** \
      MTKYQTLIFVPEGSLLNEKTAEQVALRQTLKELGHDFGPAERLKYSSLQGQVKMMGFSER \
      IALTLQNFCTDDLAEAEKIFKTKLGGQRQLVKDAIPFLDQITNQVKLILLAKEERELISA \
      RLSDSELLNYFSASYFKEDFADPLPNKNVLFQIIKEQELDPDNCLVIGTDLVEEIQGAEN \
      AGLQSLWIAPKKVKMPISPRPTLHLTKLNDLLFYLELN \
      **\>Ga0070887_12184** \
      MKGKVHLAKYETLIFILEGSLLNEKVAEQNALRQTLKLTGREYGPAERIQYNSLQEKIKL \
      LGFDERIKLTLQEFFKNDWISAKGTFYNQLQKQDQLNKDVVPFLDEVKNKVNLVLLTKEK \
      KDVASYRMQNTELINYFSAVYFKDDFACKFPNKKVLITILQQQNLTPATCLVIGTNLVDE \
      IQGAENANLDSLWLAPKKVKMPISPRPTLHLNKLTDLLFYLELS \
      **\>Ga0072400_11133** \
      LAKFQTLIFILEGSLLDEKIAEQSALKQTLKSTGRDFGPSERLKYNSVRENNKLLGFEDR \
      IQLILQTFFHENWQDAGQIFIKELQKQNRLNKEVLPFLNKVNCKVKLILLAKENKKVALQ \
      RMKNTELVNYFPFAYFKDDFTEKLPHKKVLTTILQKQNLAFATSLVIGTDLADEIQAAEN \
      AKIQSLWLAPKKVKMPISPHPTLHLNKLNDLLFYLELS
    + `cat ./database/Orthofinder/firm5_ortho_sequences/OG0001034.ffa` \
      **\>Ga0061073_1479** \
      GTGACTAAATATCAAACGTTAATTTTTGTTCCTGAAGGTAGTTTATTAAATGAAAAAACG \
      GCTGAACAAGTCGCACTCAGGCAAACTTTAAAAGAACTCGGACATGATTTTGGACCAGCT \
      GAACGCCTAAAATATTCTAGCTTACAAGGACAAGTTAAAATGATGGGTTTCAGCGAGCGC \
      ATTGCACTAACCCTGCAAAATTTTTGTACCGACGATTTGGCTGAGGCCGAAAAAATTTTC \
      AAAACAAAATTAGGAGGTCAGCGACAACTAGTCAAAGATGCTATTCCATTTCTTGACCAA \
      ATAACAAACCAAGTTAAGCTAATTCTCCTTGCCAAAGAAGAACGTGAACTAATCTCAGCT \
      CGCCTATCTGATAGCGAACTACTTAACTATTTTTCTGCTTCCTATTTTAAAGAAGATTTT \
      GCTGATCCTTTGCCAAATAAAAATGTCCTGTTTCAAATTATAAAAGAGCAAGAATTAGAT \
      CCAGATAATTGCCTAGTTATCGGCACAGATTTAGTTGAAGAAATTCAAGGAGCAGAAAAC \
      GCTGGCTTGCAATCATTATGGATTGCACCAAAAAAGGTTAAAATGCCAATTAGTCCTCGA \
      CCTACTCTGCATTTAACTAAACTCAATGACTTGCTTTTTTATCTTGAATTAAACTAG \
      **\>Ga0070887_12184** \
      ATGAAAGGAAAAGTACACTTGGCAAAATATGAAACTTTAATTTTTATTCTTGAAGGAAGC \
      TTATTAAACGAAAAAGTTGCAGAACAAAATGCACTTAGGCAAACTTTGAAATTAACTGGC \
      AGAGAATATGGTCCAGCTGAGCGCATACAATATAATTCATTACAAGAAAAGATTAAATTA \
      CTAGGATTTGATGAGCGCATTAAATTAACTTTGCAGGAATTCTTTAAAAATGACTGGATT \
      TCTGCGAAAGGCACTTTTTATAACCAGTTGCAAAAACAAGATCAGTTAAATAAAGATGTA \
      GTGCCCTTTTTAGATGAGGTGAAAAACAAAGTTAACTTGGTTTTGCTGACGAAAGAGAAA \
      AAAGATGTGGCTTCATACCGCATGCAAAATACAGAGCTAATAAATTATTTTTCCGCAGTT \
      TATTTTAAAGACGATTTTGCATGTAAGTTTCCAAATAAGAAGGTTTTGATAACAATATTG \
      CAGCAGCAGAATCTGACGCCAGCCACTTGTCTTGTAATTGGGACAAACTTAGTCGATGAA \
      ATTCAGGGTGCCGAAAATGCTAACCTGGATTCTCTTTGGCTAGCGCCCAAGAAAGTAAAA \
      ATGCCAATTAGTCCACGTCCAACTTTACATTTAAATAAATTAACTGATTTATTATTTTAC \
      CTAGAATTAAGCTAG \
      **\>Ga0072400_11133** \
      TTGGCAAAATTTCAAACATTAATTTTTATTCTTGAGGGCAGTTTATTAGATGAAAAGATT \
      GCTGAACAAAGTGCATTAAAGCAAACTTTAAAGTCAACTGGCAGAGATTTTGGTCCCAGT \
      GAACGTTTAAAATATAATTCTGTACGAGAAAATAATAAGTTGCTTGGCTTTGAAGACCGC \
      ATACAATTAATTTTACAAACATTTTTTCATGAAAATTGGCAAGATGCAGGGCAGATTTTT \
      ATCAAAGAATTACAAAAGCAAAATCGCTTGAATAAAGAAGTATTGCCATTTTTAAACAAA \
      GTTAACTGCAAGGTTAAACTAATTCTGCTGGCAAAAGAGAACAAAAAAGTAGCATTACAG \
      CGCATGAAGAACACAGAGTTGGTAAATTATTTTCCGTTTGCTTATTTTAAAGATGACTTT \
      ACGGAAAAATTGCCACATAAAAAAGTTTTGACCACCATTTTGCAGAAACAAAACTTGGCG \
      TTCGCAACTAGTTTAGTAATCGGAACTGACTTAGCAGATGAAATTCAGGCTGCAGAGAAT \
      GCCAAAATACAGTCACTCTGGCTAGCGCCTAAGAAAGTAAAAATGCCGATTAGCCCGCAC \
      CCAACTTTACATTTAAATAAATTAAACGATTTATTATTTTACCTAGAATTAAGCTAG
  + `rule calc_perc_id`
    + this rule relies on various tools and scripts tied together by `scripts/aln_calc.sh`. The scripts are:
      + [mafft](https://mafft.cbrc.jp/alignment/software/manual/manual.html#index) for alignment
        + The aligned result is in a multi-fasta file called {OrthogroupID}_aln.fasta.
        EXAMPLE:
        + `cat ./database/Orthofinder/firm5_ortho_sequences/OG0001034_aln.fasta` \
          **\>Ga0061073_1479** \
          -------MTKYQTLIFVPEGSLLNEKTAEQVALRQTLKELGHDFGPAERLKYSSLQGQVK \
          MMGFSERIALTLQNFCTDDLAEAEKIFKTKLGGQRQLVKDAIPFLDQIT---NQVKLILL \
          AKEERELISARLSDSELLNYFSASYFKEDFADPLPNKNVLFQIIKEQELDPDNCLVIGTD \
          LVEEIQGAENAGLQSLWIAPKKVKMPISPRPTLHLTKLNDLLFYLELN \
          **\>Ga0070887_12184** \
          M-KGKVHLAKYETLIFILEGSLLNEKVAEQNALRQTLKLTGREYGPAERIQYNSLQEKIK \
          LLGFDERIKLTLQEFFKNDWISAKGTFYNQLQKQDQLNKDVVPFLDEVK---NKVNLVLL \
          TKEKKDVASYRMQNTELINYFSAVYFKDDFACKFPNKKVLITILQQQNLTPATCLVIGTN \
          LVDEIQGAENANLDSLWLAPKKVKMPISPRPTLHLNKLTDLLFYLELS \
          **\>Ga0072400_11133** \
          -------LAKFQTLIFILEGSLLDEKIAEQSALKQTLKSTGRDFGPSERLKYNSVRENNK \
          LLGFEDRIQLILQTFFHENWQDAGQIFIKELQKQNRLNKEVLPFLNKVN---CKVKLILL \
          AKENKKVALQRMKNTELVNYFPFAYFKDDFTEKLPHKKVLTTILQKQNLAFATSLVIGTD \
          LADEIQAAENAKIQSLWLAPKKVKMPISPHPTLHLNKLNDLLFYLELS
      + **aln_aa_to_dna.py**
        + This scripts converts the alignments into nucleotide sequences rather than amino acid sequences
        EXAMPLE:
        + `cat ./database/Orthofinder/firm5_ortho_sequences/OG0001034_aln.fasta`
        **\>Ga0061073_1479** \
        ---------------------GTGACTAAATATCAAACGTTAATTTTTGTTCCTGAAGGT \
        AGTTTATTAAATGAAAAAACGGCTGAACAAGTCGCACTCAGGCAAACTTTAAAAGAACTC \
        GGACATGATTTTGGACCAGCTGAACGCCTAAAATATTCTAGCTTACAAGGACAAGTTAAA \
        ATGATGGGTTTCAGCGAGCGCATTGCACTAACCCTGCAAAATTTTTGTACCGACGATTTG \
        GCTGAGGCCGAAAAAATTTTCAAAACAAAATTAGGAGGTCAGCGACAACTAGTCAAAGAT \
        GCTATTCCATTTCTTGACCAAATAACA---------AACCAAGTTAAGCTAATTCTCCTT \
        GCCAAAGAAGAACGTGAACTAATCTCAGCTCGCCTATCTGATAGCGAACTACTTAACTAT \
        TTTTCTGCTTCCTATTTTAAAGAAGATTTTGCTGATCCTTTGCCAAATAAAAATGTCCTG \
        TTTCAAATTATAAAAGAGCAAGAATTAGATCCAGATAATTGCCTAGTTATCGGCACAGAT \
        TTAGTTGAAGAAATTCAAGGAGCAGAAAACGCTGGCTTGCAATCATTATGGATTGCACCA \
        AAAAAGGTTAAAATGCCAATTAGTCCTCGACCTACTCTGCATTTAACTAAACTCAATGAC \
        TTGCTTTTTTATCTTGAATTAAAC \
        **\>Ga0070887_12184** \
        ATG---AAAGGAAAAGTACACTTGGCAAAATATGAAACTTTAATTTTTATTCTTGAAGGA \
        AGCTTATTAAACGAAAAAGTTGCAGAACAAAATGCACTTAGGCAAACTTTGAAATTAACT \
        GGCAGAGAATATGGTCCAGCTGAGCGCATACAATATAATTCATTACAAGAAAAGATTAAA \
        TTACTAGGATTTGATGAGCGCATTAAATTAACTTTGCAGGAATTCTTTAAAAATGACTGG \
        ATTTCTGCGAAAGGCACTTTTTATAACCAGTTGCAAAAACAAGATCAGTTAAATAAAGAT \
        GTAGTGCCCTTTTTAGATGAGGTGAAA---------AACAAAGTTAACTTGGTTTTGCTG \
        ACGAAAGAGAAAAAAGATGTGGCTTCATACCGCATGCAAAATACAGAGCTAATAAATTAT \
        TTTTCCGCAGTTTATTTTAAAGACGATTTTGCATGTAAGTTTCCAAATAAGAAGGTTTTG \
        ATAACAATATTGCAGCAGCAGAATCTGACGCCAGCCACTTGTCTTGTAATTGGGACAAAC \
        TTAGTCGATGAAATTCAGGGTGCCGAAAATGCTAACCTGGATTCTCTTTGGCTAGCGCCC \
        AAGAAAGTAAAAATGCCAATTAGTCCACGTCCAACTTTACATTTAAATAAATTAACTGAT \
        TTATTATTTTACCTAGAATTAAGC \
        **\>Ga0072400_11133** \
        ---------------------TTGGCAAAATTTCAAACATTAATTTTTATTCTTGAGGGC \
        AGTTTATTAGATGAAAAGATTGCTGAACAAAGTGCATTAAAGCAAACTTTAAAGTCAACT \
        GGCAGAGATTTTGGTCCCAGTGAACGTTTAAAATATAATTCTGTACGAGAAAATAATAAG \
        TTGCTTGGCTTTGAAGACCGCATACAATTAATTTTACAAACATTTTTTCATGAAAATTGG \
        CAAGATGCAGGGCAGATTTTTATCAAAGAATTACAAAAGCAAAATCGCTTGAATAAAGAA \
        GTATTGCCATTTTTAAACAAAGTTAAC---------TGCAAGGTTAAACTAATTCTGCTG \
        GCAAAAGAGAACAAAAAAGTAGCATTACAGCGCATGAAGAACACAGAGTTGGTAAATTAT \
        TTTCCGTTTGCTTATTTTAAAGATGACTTTACGGAAAAATTGCCACATAAAAAAGTTTTG \
        ACCACCATTTTGCAGAAACAAAACTTGGCGTTCGCAACTAGTTTAGTAATCGGAACTGAC \
        TTAGCAGATGAAATTCAGGCTGCAGAGAATGCCAAAATACAGTCACTCTGGCTAGCGCCT \
        AAGAAAGTAAAAATGCCGATTAGCCCGCACCCAACTTTACATTTAAATAAATTAAACGAT \
        TTATTATTTTACCTAGAATTAAGC
      + **trim_aln.py** and `sed` to simplify headers to contain just genome ID and leave out gene identifier (as they are all single copy core genes).
        + This script trims out all the sections that do not align by counting which positions have "-" and removing those from all the members of the orthogroup.
        + `cat ./database/Orthofinder/firm5_ortho_sequences/OG0001034_aln_trim.fasta` \
          **\>Ga0061073** \
          GTGACTAAATATCAAACGTTAATTTTTGTTCCTGAAGGTAGTTTATTAAATGAAAAAACG \
          GCTGAACAAGTCGCACTCAGGCAAACTTTAAAAGAACTCGGACATGATTTTGGACCAGCT \
          GAACGCCTAAAATATTCTAGCTTACAAGGACAAGTTAAAATGATGGGTTTCAGCGAGCGC \
          ATTGCACTAACCCTGCAAAATTTTTGTACCGACGATTTGGCTGAGGCCGAAAAAATTTTC \
          AAAACAAAATTAGGAGGTCAGCGACAACTAGTCAAAGATGCTATTCCATTTCTTGACCAA \
          ATAACAAACCAAGTTAAGCTAATTCTCCTTGCCAAAGAAGAACGTGAACTAATCTCAGCT \
          CGCCTATCTGATAGCGAACTACTTAACTATTTTTCTGCTTCCTATTTTAAAGAAGATTTT \
          GCTGATCCTTTGCCAAATAAAAATGTCCTGTTTCAAATTATAAAAGAGCAAGAATTAGAT \
          CCAGATAATTGCCTAGTTATCGGCACAGATTTAGTTGAAGAAATTCAAGGAGCAGAAAAC \
          GCTGGCTTGCAATCATTATGGATTGCACCAAAAAAGGTTAAAATGCCAATTAGTCCTCGA \
          CCTACTCTGCATTTAACTAAACTCAATGACTTGCTTTTTTATCTTGAATTAAAC \
          **\>Ga0070887** \
          TTGGCAAAATATGAAACTTTAATTTTTATTCTTGAAGGAAGCTTATTAAACGAAAAAGTT \
          GCAGAACAAAATGCACTTAGGCAAACTTTGAAATTAACTGGCAGAGAATATGGTCCAGCT \
          GAGCGCATACAATATAATTCATTACAAGAAAAGATTAAATTACTAGGATTTGATGAGCGC \
          ATTAAATTAACTTTGCAGGAATTCTTTAAAAATGACTGGATTTCTGCGAAAGGCACTTTT \
          TATAACCAGTTGCAAAAACAAGATCAGTTAAATAAAGATGTAGTGCCCTTTTTAGATGAG \
          GTGAAAAACAAAGTTAACTTGGTTTTGCTGACGAAAGAGAAAAAAGATGTGGCTTCATAC \
          CGCATGCAAAATACAGAGCTAATAAATTATTTTTCCGCAGTTTATTTTAAAGACGATTTT \
          GCATGTAAGTTTCCAAATAAGAAGGTTTTGATAACAATATTGCAGCAGCAGAATCTGACG \
          CCAGCCACTTGTCTTGTAATTGGGACAAACTTAGTCGATGAAATTCAGGGTGCCGAAAAT \
          GCTAACCTGGATTCTCTTTGGCTAGCGCCCAAGAAAGTAAAAATGCCAATTAGTCCACGT \
          CCAACTTTACATTTAAATAAATTAACTGATTTATTATTTTACCTAGAATTAAGC \
          **\>Ga0072400** \
          TTGGCAAAATTTCAAACATTAATTTTTATTCTTGAGGGCAGTTTATTAGATGAAAAGATT \
          GCTGAACAAAGTGCATTAAAGCAAACTTTAAAGTCAACTGGCAGAGATTTTGGTCCCAGT \
          GAACGTTTAAAATATAATTCTGTACGAGAAAATAATAAGTTGCTTGGCTTTGAAGACCGC \
          ATACAATTAATTTTACAAACATTTTTTCATGAAAATTGGCAAGATGCAGGGCAGATTTTT \
          ATCAAAGAATTACAAAAGCAAAATCGCTTGAATAAAGAAGTATTGCCATTTTTAAACAAA \
          GTTAACTGCAAGGTTAAACTAATTCTGCTGGCAAAAGAGAACAAAAAAGTAGCATTACAG \
          CGCATGAAGAACACAGAGTTGGTAAATTATTTTCCGTTTGCTTATTTTAAAGATGACTTT \
          ACGGAAAAATTGCCACATAAAAAAGTTTTGACCACCATTTTGCAGAAACAAAACTTGGCG \
          TTCGCAACTAGTTTAGTAATCGGAACTGACTTAGCAGATGAAATTCAGGCTGCAGAGAAT \
          GCCAAAATACAGTCACTCTGGCTAGCGCCTAAGAAAGTAAAAATGCCGATTAGCCCGCAC \
          CCAACTTTACATTTAAATAAATTAAACGATTTATTATTTTACCTAGAATTAAGC
      + **calc_perc_id_orthologs.py**
        + Uses as input, trimmed aligned sequences and a metafile (`database/genome_db_210402_metafile.txt`) which is a tab-delim file with genome-id in tab1 and SDP-affiliation in tab 3
        + First, it checks the number of SDPs contained within the alignment. If more than one, it continues by calculating alignment percentage identity stats across SDPs. If only one SDP, exits script.
        + Next, it Compares the genomes in each SDP to all other genomes in the alignment: calculates percentage identity for all pairwise combinations. Calculates the max, min, and mean values, prints to file `database/Orthofinder/{phylotype}_perc_id.txt` showing one orthogroup per line.
        EXAMPLE:
        + `cat ./database/Orthofinder/firm5_perc_id.txt` \
          ... \
          OG0001034	0.674	0.586	0.972 \
          ... \

  + `rule filter_orthogroups`
    + orthogroups are filtered based on:
      - Minimum gene-length 300bp (applied to all members of each gene-family)
      - Inter-SDP max alignment identity 95% (only if the phylotype contain multiple SDPs)
    + Short genes are filtered off, because they are likley to be less reliable for accurate coverage estimates. Similarly, the inter-SDP similarity threshold is used to ensure that there is enough divergence between the SDPs for reliable mapping (at least as estimated from the currently availabe genomes). It is worthwhile to check the number of gene-families before/after filtering. If a lot of gene-families were filtered off, this could be an indication that the SDPs are not properly discrete.
    + This finally results in the single-copy core genes that have been filtered to be used for core coverage estimation. `database/Orthofinder/{phylotype}_single_ortho_filt.txt`.
  + `rule core_cov`
    + takes as input bam files with the alignments for each sample to be considered (as a text file containing a list of these files) and the _\_single_ortho_filt_ file. Outputs are written to `04_CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov.txt`.
    + The script reads the filtered orthofile `database/Orthofinder/{phylotype}_single_ortho_filt.txt` and gets the gene-famililes and genome-ids for each SDP.
    + Then, from bed files, it finds the start and end positions of each of the genes in an orthogroup for each of the genomes of the orthogroup. It writes these to the file, `04_CoreCov_*/{phylotype}_corecov.txt` each SDP in the phylotype, start position for each gene family in the genome marked reference for that SDP.
    + The coverage is also written to this file.
    Example:
      + SDP	Sample	OG	Ref_pos	Coverage \
        firm5_1	DrY1_N1_microbiome_mapped	OG0000932	448	18.81 \
        firm5_1	DrY1_N1_microbiome_mapped	OG0000931	1991	23.34 \
        ... \
        firm5_2	M1.5_microbiome_mapped	OG0000935	1852405	12.29 \
        firm5_2	M1.5_microbiome_mapped	OG0000934	1853270	9.95 \
        firm5_3	DrY1_N1_microbiome_mapped	OG0000932	1	501.61 \
        firm5_3	DrY1_N1_microbiome_mapped	OG0000931	1542	534.77 \
        ... \
        firm5_bombus	M1.5_microbiome_mapped	OG0000936	1674767	0.0 \
        firm5_bombus	M1.5_microbiome_mapped	OG0000935	1676256	0.0 \
        firm5_bombus	M1.5_microbiome_mapped	OG0000934	1677124	0.0
    + SDP abundance is estimated based on mapped read coverage of core genes. It sums up gene coverages of all the genes og OG families associated with said SDP across genomes belogining to the SDP.
    + It also reports PTR (Peak-Trough Ratio).
    + Most species in the database are represented by multiple genomes (< 98.5% gANI between genomes). Core genes are inferred at the phylotype. More accurate estimates can be obtained by using a large number (+700) of core genes.
  + `rule core_cov_plots`
    + This R-script will estimate the coverage at the terminus, using the summed core gene family coverages. If the cov-ter cannot be properly estimated (fx. due to draft genome status or lack of replication), an estimate will be generated using the median coverage across core gene families, and the PTR is set to NA. If more than 20\% of the core gene families have no coverage, the abundance will be set to zero. As output, a tabular file is generated (including the cov-ter/median cov, and PTR), and a pdf-file with plots for visual validation.
    + First, filter for samples with coverage of at least 1 on > 80% of the core genes. Next, values that are deviating no more than 2 times the median are kept others are discarded as outliers.
    + Next, gets fitted coordinates and append values to coord-table. It does this by using the segmented package. As explained below,

  ```{r eval=FALSE}
      x <- data_filt$Ref_pos
   	  y <- data_filt$Coverage
   	  psi_est <- max(x)/2
      lin.mod = y ~ x
      segmented(lin.mod, seg.Z=~x, psi=psi_est)
  ```
  where, `lin.mod` is a simple linear model that was made by base R. The R package `Segmented` supports breakpoint analysis. The methods used by this package are applicable when segments are (nearly continuous) so this means that for the regression to make sense the core gene families selected should cover the reference genome well and without too many huge gaps. `psi`, is a starting value of the breakpoint. Example of a model fit using segemented,
  ```{r eval=FALSE}
      Call: segmented.lm(obj = lin.mod, seg.Z = ~x, psi = psi_est)

      Meaningful coefficients of the linear terms:
      (Intercept)            x         U1.x
        1.332e+02   -5.984e-05    1.243e-04

      Estimated Break-Point(s):
      psi1.x
      860065
  ```
  `x` is the slope of the first segment and `U1.x` is the difference in slopes between the first and second segment. `psi_est` is the newly estimated breakpoint. This along with the slopes
      - The summary function shows:
  ```{r eval=FALSE}
  ***Regression Model with Segmented Relationship(s)***

  Call:
  segmented.lm(obj = lin.mod, seg.Z = ~x, psi = psi_est)

  Estimated Break-Point(s):
      Est.   St.Err
  psi1.x 860065 9421.565

  Meaningful coefficients of the linear terms:
      Estimate Std. Error t value Pr(>|t|)
  (Intercept)  1.332e+02  9.672e-01  137.75   <2e-16 ***
  x           -5.984e-05  1.782e-06  -33.57   <2e-16 ***
  U1.x         1.243e-04  2.688e-06   46.24       NA
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

  Residual standard error: 8.804 on 750 degrees of freedom
  Multiple R-Squared: 0.7417,  Adjusted R-squared: 0.7407

  Convergence attained in 3 iter. (rel. change 2.8661e-06)
  ```
  Finally, from the segmented model the ptr is calculated as follows:
  ```{r eval=FALSE}
  cov_ter <- round(slope1*psi + intercept1, digits=1)
  cov_ori2  <- slope2*(tail(x, n=1)) + intercept2
  max_ori_cov <- max(intercept1, cov_ori2)
  min_ori_cov <- min(intercept1, cov_ori2)
  if ((psi<psi_min) || (psi>psi_max) || (min_ori_cov<cov_ter)){
      ptr <- NA
  }  else {
      ptr <- round(max_ori_cov/cov_ter, digits=2)
  }

  # where,
  psi <- (summary.segmented(seg.mod)[[12]])[2]
  psi_est <- max(x)/2
  psi_est <- max(x)/2
  psi_min <- psi_est-(0.5*psi_est)
  psi_max <- psi_est+(0.5*psi_est)
  ```
  For `cor_ter` is the coverage $y = ax + b$ where x is psi (the breakpoint on the x axis) and the `cor_ori2` is the coverage at the ori which is the section with maximum coverage and at the `tail`. The condition `(psi<psi_min) || (psi>psi_max) || (min_ori_cov<cov_ter)` checks: **First**; If the break-point is too far from the expected place (+/-50% of break-point estimate), ptr is set to `NA`. **Second**; If the coverage at ori (either beginning or end of dataframe) is lower than the estimated coverage at ter, ptr is also set to `NA`. Finally, if the coverage of the origin is not greater than the terminus, ptr is set to `NA`.

    + the PTR was set to `NA`, the median will be plotted and used for quantification. Else, the segmented regression line is plotted, and the terminus coverage is used for quantification.

  ## Scripts
<!-- This list is not exhaustive. Make a complete list #TODO -->
  + `aln_aa_to_dna.py`
  + `aln_calc.sh`
  + `calc_perc_id_orthologs.py`
  + `core_cov.py`
  + `core_cov.R`
  + `download_data.py`
  + `extract_orthologs.py`
  + `fasta_generate_regions.py`
  + `filter_bam.py`
  + `filter_orthologs.py`
  + `filter_sam_aln_length.pl`
  + `filter_sam_aln_length_unmapped.pl`
  + `filter_snvs.pl`
  + `filt_vcf_samples.pl`
  + `get_single_ortho.py`
  + `parse_core_cov.py`
  + `parse_spades_metagenome.pl`
  + `rearange_faa.py`
  + `subset_orthofile.py`
  + `trim_aln.py`
  + `scripts/write_adapters.py`

  ## Envs

`core-cov-env.yaml`
`mapping-env.yaml`
`rmd-env.yaml`
`snv-env.yaml`
`trim-qc-env.yaml`
...

`snakemake -p --use-conda --conda-prefix /scratch/aprasad/built-envs/ --conda-frontend mamba --profile slurm --restart-times 0 --cluster-cancel scancel --rerun-incomplete --keep-going --rerun-triggers mtime -r --scheduler greedy`

# Introduction

The analysis is split into multiple parts.

1. A database-dependent mapping and community profiling at the SDP and Phylotype level.
2. An assembly-based approach which involes MAG binning per sample and dereplication.
Requires manual annotation of input for the next section
3. Core genome phylogeny construction of MAGs and isolate genomes.
  (a) All high-quality MAGs
  (b) One MAG per magOTU cluster
4. Mapping and SNV calling of reads to
  (a) Just MAGs
  (b) A combination of isolates and MAGs

The samples used in this analysis are mentioned below.

* "M1.1", "M1.2", "M1.3", "M1.4", "M1.5",
  * _Apis mellifera_ from India
* "DrY2_F1","DrY2_F2",
  * _Apis mellifera_ from Switzerland
* "AmAi02","AmIu02",
  * _Apis mellifera_ from Japan
* "C1.1", "C1.2", "C1.3", "C1.4", "C1.5",
  * _Apis cerana_ from India, colony number 1
* "C2.1", "C2.2", "C2.3", "C2.4", "C2.5",
  * _Apis cerana_ from India, colony number 2
* "C3.1", "C3.2", "C3.3", "C3.4", "C3.5",
  * _Apis cerana_ from India, colony number 1
* "AcCh05","AcKn01",
  * _Apis cerana_ from Japan
* "D1.1","D1.2","D1.3","D1.4","D1.5",
  * _Apis dorsata_ from India, colony number 1
* "D2.1","D2.2","D2.3","D2.4","D2.5",
  * _Apis dorsata_ from India, colony number 2
* "D3.1","D3.2","D3.3","D3.4","D3.5",
  * _Apis dorsata_ from India, colony number 3
* "F1.1","F1.2","F1.3","F1.4","F1.5",
  * _Apis florea_ from India, colony number 1
* "F2.1","F2.2","F2.3","F2.4","F2.5",
  * _Apis florea_ from India, colony number 2
* "F3.1","F3.2","F3.3","F3.4","F3.5"
  * _Apis florea_ from India, colony number 3

## Data summary

50 samples were first mapped to a host database and then the unmapped reads to a microbiome database.

The raw reads and trimmed reads were checked for quality using the tool [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). The report (html format) also summarises basic statistics including number of reads. The QC results can be found in their respective folders at `fastqc/raw/{SAMPLE}_R*_fastqc.html` and `fastqc/trim/{SAMPLE}_R*_trim_fastqc.html`.

# Data description

The dataset comprises samples from 4 species of honey bees sampled in India.

+ _Apis mellifera_ (5 individuals from 1 colony)
+ _Apis cerana_ (5 individuals from 3 colonies)
+ _Apis dorsata_ (5 individuals from 3 colonies)
+ _Apis florea_ (5 individuals from 3 colonies)

Some samples from older publications of the lab are also included for _Apis mellifera_ and _Apis cerana_

The data is stored in various locations as described below and backed up on the NAS.

+ **Raw data**:

	- NAS recerche:
	`/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/211102_Medgenome_india_samples_resequenced`
	`/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/211018_Medgenome_india_samples`
	(cluster - aprasad@curnagl.dcsr.unil.ch)

  -	NAS:
	`/home/aiswarya/mnt/aprasad/SPIRIT_Project/Data/RawData/211018_Medgenome_india_samples.tar.gz`
	`/home/aiswarya/mnt/aprasad/SPIRIT_Project/Data/RawData/211102_Medgenome_india_samples_resequenced.tar.gz`
	`/home/aiswarya/mnt/lab_resources/NGS_data/20211018_A01223-105-HC32VDSX2/`
	`/home/aiswarya/mnt/lab_resources/NGS_data/20211102_A01223-105-HC32VDSX2/`
	(workstation - aiswarya@130.223.110.124)

+ **Trimmed data**:

	`/work/FAC/FBM/DMF/pengel/spirit/aprasad/211018_Medgenome_india_samples/01_Trimmed/`
	`/work/FAC/FBM/DMF/pengel/spirit/aprasad/211102_Medgenome_india_samples_resequenced/01_Trimmed/`
	(cluster - aprasad@curnagl.dcsr.unil.ch)

+ **Working directory backup**:
  (need to keep up-to-date using script on the cluster `bash ~/backup_workdir.sh` and logs are writted to ~/yymmdd_backup_log on the cluster)

	`/home/aiswarya/mnt/aprasad/Backups/working_dir_backup/Cluster/211102_Medgenome_india_samples_resequenced`
	`/home/aiswarya/mnt/aprasad/Backups/working_dir_backup/Cluster/211018_Medgenome_india_samples`
  (workstation - aiswarya@130.223.110.124)

+ **Results and important intermediate files**:

	`/home/aiswarya/mnt/aprasad/SPIRIT_Project/Data/211018_Medgenome_india_analysis`
  (workstation - aiswarya@130.223.110.124)


+ **Conda installation (cluster)**:

	`/work/FAC/FBM/DMF/pengel/spirit/aprasad/Miniconda3`
  (cluster - aprasad@curnagl.dcsr.unil.ch)

### Nomenclature

There are 56 samples at the moment.

+ M1.1 - M1.5 are 5 individuals of _Apis mellifera_ from colony 1
+ Cx.1 - Cx.5 are 5 individuals of _Apis cerana_ from colony x for 3 colonies (1 - 3)
+ Dx.1 - Dx.5 are 5 individuals of _Apis dorsata_ from colony x for 3 colonies (1 - 3)
+ Fx.1 - Fx.5 are 5 individuals of _Apis florea_ from colony x for 3 colonies (1 - 3)
+ DrY2_F1 and DrY2_F2 are samples from KE's 2015 paper. _Apis mellifera_ from switzerland (Les Droites)
+ AcCh05, AcKn01 and AmAi02, AmIu02 are two samples of _Apis cerana_ and _Apis mellifera_ from different apiaries in Japan

These samples were from earlier runs:

  + **20151119_WINDU89**	20151119	Kirsten_Ellegaard	6	GTF	Illumina	100	PE	HiSeq 2500	Genomic diversity landscape of the honey bee gut microbiota (2019, NatCom)	Nurses, Year 1, Les Droites
  + **20160415_OBIWAN225**	20160415	Kirsten_Ellegaard	12	GTF	Illumina	100	PE	HiSeq 2500	Genomic diversity landscape of the honey bee gut microbiota (2019, NatCom)	Foragers/Winterbees, Year 1, Les Droites \
  + **20161216_OBIWAN275**	20161216	Kirsten_Ellegaard	6	GTF	Illumina	100	PE	HiSeq 2500	Genomic diversity landscape of the honey bee gut microbiota (2019, NatCom)	Nurses, Year 2, Les Droites \
  + **20170310_WINDU179**	20170310	Kirsten_Ellegaard	12	GTF	Illumina	100	PE	HiSeq 2500	Genomic diversity landscape of the honey bee gut microbiota (2019, NatCom)	Foragers/Winterbees, Year 2, Les Droites (**INCLUDED FOR NOW**) \
  + **20170426_OBIWAN300**	20170426	Kirsten_Ellegaard	6	GTF	Illumina	100	PE	HiSeq 2500	Genomic diversity landscape of the honey bee gut microbiota (2019, NatCom)	Nurses, Year 2, Grammont \
  + **20170428_WINDU191**	20170428	Kirsten_Ellegaard	12	GTF	Illumina	100	PE	HiSeq 2500	Genomic diversity landscape of the honey bee gut microbiota (2019, NatCom)	Foragers/Winterbees, Year 2, Grammont \
  + **20180118_OBIWAN338-339**	20180118	Kirsten_Ellegaard	30	GTF	Illumina	100	PE	HiSeq 2500	Metagenomes of individual honey bees, subjected to dietary manipulation and kept in the lab \
  + **20180612_KE_japan_metagenomes**	20180612	Ryo_Miyasaki	40	Japan	Illumina	100	PE	HiSeq 2500	Vast differences in strain-level diversity in two closely related species of honey bees (2020, CurBiol)	Sampling and sequencing done in Japan (**INCLUDED FOR NOW**)

## Databases

### Host database

The database is named **4_host_db**.

A [paper](https://academic.oup.com/gbe/article/12/1/3677/5682415) published in Dec. 2019 a high quality [_Apis dorsata_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCA_009792835.1/) as an improvement over a previous submission in 2013. The paper also mentioned studies that had previously sequenced the [_Apis florea_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCA_000184785.2) in 2012, [_Apis cerana_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_001442555.1) in 2015 (other assemblies submitted found [here](https://www.ncbi.nlm.nih.gov/assembly/organism/7460/latest/)) and [_Apis mellifera_ genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_003254395.2) in 2018 (other assemblies submitted listed here). So far I have not found any whole genome assemblies of _Apis adreniformis_.

These assemblies were downloaded and concatenated to make the **4_host_db**. It contains,

+ `>apis_mellifera_2018 PRJNA471592 version Amel_Hac3.1`
+ `>Apis_cerana  PRJNA235974`
+ `>Apis_cerana_mitochondrion PRJNA235974`
+ `>Apis_florea PRJNA45871`
+ `>Apis_dorsata PRJNA174631`


### Microbiome database

The database is named **genome_db_210402**.

This database was set up by Dr Kirsten Ellegard (KE) and is described on zenodo. It uses NCBI and IMG genome assemblies. It is non-redundant and contains concatenated genomes. Located at in lab NAS directory at lab_resources/Genome_databse. In this pipeline so far, the version of the pipeline set up by KE’s community profiling pipeline.

It was downloaded by the script `download.py --genome_db` from [zenodo](https://zenodo.org/record/4661061#.YcGlSxPMK3I). This dowloads multiple directories. The Orthofider directory can be deleted as this is generated for the pipeline as needed. The bed files can be generated from gff files if desired but this was already done for the genomes of that database so was not repeated. The other files (ffn, gff) are found in the public repository from where the genome was downloaded. The faa files were reorganised in directories corresponding to their respective SDPs in order to allow the Orthofinder scripts to assign orthogroups per SDP.

These repositories follow their own annotation pipeline to generate these files. The database can also be found at `<NAS>/lab_resources/Genome_database/database_construction`. It contains 198 genomes identified by their locus tags and described in `<NAS>/lab_resources/Genome_database/database_construction/database_construction` in metadata sheets.

<!-- 2.0.5.2 honeybee_gut_microbiota_db_red
Generated by subsetting genomes and Orthofiles from genome_db_210402 for the sake of SNV analysis. So there is only one genome per SDP. The choice of which genome is made by Kirsten’s pipeline. later, review to see if it is the best choice for this pipleine as well. -->

<!-- 2.0.5.4 Orthofinder files
In preparation,
In the last step, all the database files will be generated, based on the database metadata-file (i.e. "Locustag_Phylotype_SDP_final.txt"). The following three steps must be performed:

1. Use the bash-script "get_selected_genome_files.sh" will copy the relevant genome data-files into four database dirs "faa_files", "ffn_files", "fna_files", "gff_files").
2. Generate concatenate versions of the each genome assembly fasta-file will be generated, and bed-files detailing the gene positions on these concatenate genomes are generated ("generate_bed_concat.py")
3. Infer single-copy core gene families

Make Orthofinder files from genomes in database grouped by phylotypes. There is a metafile in the database directory which was made by Gilles and was copied from his NAS directory. This is first used to rename genomes. The genomes in the database are named using locus tags. Make a script to rename them (and modify related files accordingly later).

First, get a list of all phylotypes using the output of cat all_genomes_metafile.tsv | cut -f7 | uniq. The output contains some SDP names. Replace this by the corresponding phylotype (choose more appropriate way to handle this later). The resulting list is declared as CorePhylos in the Snakefile (Not core phylotypes but rather list for core coverage calculation for all phylos!). This information is also mostly available in the metafile that comes with the database but the one copied from Gilles includes more information and also some other genomes (added by German and mentioned in all_genomes_metafile.tsv).

CorePhylos = ["Apibacter", "Bartonella", "Bifido", "Bombella", "Commensalibacter", "Firm4", "Firm5", "Frischella", "Gilliamella", "Lkunkeei", "Snodgrassella"]



Run orthofinder on the genomes grouped by phylotype. The -og flag says to stop after inferring orthogroups and avoids further unecessary computation. -f specifies to start analysis from directory containing FASTA files. Then get the single-copy core genes using the script get_single_copy.py

Next, filter and continue to re run core cov and then proceed to snv calling and filtering! -->

# Description of pipeline methods

Run entire snakemake pipeline using:

`snakemake -p --use-conda --conda-prefix /work/FAC/FBM/DMF/pengel/spirit/aprasad/Miniconda3/spirit_envs --conda-frontend conda --profile slurm --restart-times 0 --keep-going`

and if resuming a failed or stopped run, use:

`snakemake -p --use-conda --conda-prefix /work/FAC/FBM/DMF/pengel/spirit/aprasad/Miniconda3/spirit_envs --conda-frontend conda --profile slurm --restart-times 0 --keep-going --rerun-incomplete`

`snakemake -p --use-conda --conda-prefix /scratch/aprasad/built-envs/ --conda-frontend mamba --profile slurm --restart-times 0 -r --cluster-cancel scancel --keep-going --rerun-incomplete --rerun-triggers mtime` (as of 22/01/23)

conda environments are all specified in `envs/` and built by snakemake under various names in `/work/FAC/FBM/DMF/pengel/spirit/aprasad/Miniconda3/spirit_envs`

Run the pipeline in the conda environment called `snakmake_with_samtools` in the cluster. It is a clone of the snakemake environment made as recommended by Snakemake [docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba) followed by `conda install biopython` and later `conda install samtools` in it. This is so that Kirsten's core_cov script works (specific conda environments can only be specified for rules using bash).

## Description of directory structure

Directory names are largely self-explanatory.

>`00_rawdata`, `01_Trimmed`, `02_HostMapping`, `03_MicrobiomeMapping`
>`database` contains databases to be used for mapping. It also contains `Orthofinder` files. These are described later in the sections describing associated rules.
>`envs` contains all yaml files required for this pipeline. They contain a list of packages needed to specify the conda environment for various rules to work within.
>`logs` contains log files
>`scripts` contains all scripts needed for the snakemake pipeline. Many of these scripts are adapted from Kirsten's scripts from the zenodo directories, github or from the lab_resources directories.
The **results** of the core coverage estimation are stored in,
> `04_CoreCov_211018_Medgenome_india_samples`
> `07_SNVProfiling` is not fully implemented (yet) for these samples as it is not relevant at this time.
>`fastqc` contains fastqc **results** for trimmed and raw files
+ bamfile_list_red.txt - required by KE's core coverage pipeline
+ bamfile_list.txt - required by KE's core coverage pipeline
+ Adapters-PE.fa - is generated based on index sequences by the script `scripts/write_adapters.py` (was deleted earlier as it was on scratch. Needs to be re-written.)
+ config.yaml - comprises information including list of samples
+ index_table.csv - used by the script `scripts/write_adapters.py` to make indexed adapters
+ Mapping_summary.csv - result from the rule summarize_mapping
+ rulegraph.pdf - summary DAG of rules in the pipeline (made using `snakemake --forceall --rulegraph | dot -Tpdf > Figuers/rulegraph.pdf`)
+ Report.Rmd - this report !
+ Report.html - this report compiled !
+ Snakefile - the pipipeline !!!

# Next steps

It is clear that the database is not best suited for some SDPs found especially in host species other than _Apis mellifera_. So, the next step would be to implement a MAG based analysis to compare these samples. However, as the database was already shown to be well-suited for _Apis mellifera_ and _Apis cerana_, another set of analysis would compare these samples from the [publication](https://www.sciencedirect.com/science/article/pii/S0960982220305868) ([zenodo](https://zenodo.org/record/3747314#.YcGkvRPMK3I)) with the samples from India.

## Choosing representative genomes

Information in [drep](https://drep.readthedocs.io/en/latest/choosing_parameters.html) documentation.

$$A * Completeness - B * Contamination + C * (Contamination * frac{strainheterogeneity}{100}) + D * log(N50) + E * log(size) + F * (centrality - S_{a}ni) $$s

"It makes no sense to perform bootstrap with less than 4 sequences" from IQTree

## extra info to be curated later

meaning of benchmark output

colname | type (unit) | description
s | float (seconds) | Running time in seconds
h:m:s | string (-) | Running time in hour, minutes, seconds format
max_rss | float (MB) | Maximum "Resident Set Size”, this is the non-swapped physical memory a process has used.
max_vms | float (MB) | Maximum “Virtual Memory Size”, this is the total amount of virtual memory used by the process
max_uss | float (MB) | “Unique Set Size”, this is the memory which is unique to a process and which would be freed if the process was terminated right now.
max_pss | float (MB) | “Proportional Set Size”, is the amount of memory shared with other processes, accounted in a way that the amount is divided evenly between the processes that share it (Linux only)
io_in | float (MB) | the number of MB read (cumulative).
io_out | float (MB) | the number of MB written (cumulative).
mean_load | float (-) | CPU usage over time, divided by the total running time (first row)
cpu_time | float(-) | CPU time summed for user and system

[how to setup snakemake profile](https://github.com/RomainFeron/snakemake-slurm)

### Information about MAGs to report (MIMAG)

* type_of_genome: MAG or Isolate
* assembly_software: metaspades #get version
* annotation_tool: prokka #get version
* assembly_quality:
  * finished
  * high_quality_draft
  * medium_quality_draft
  * low_quality_draft
* completeness_score:
* contamination_score:
* completeness_software:
  * checkm
* number_of_contigs:
* 16_recovered:
  * yes/no
* 16_recovery_software:
  * checkm
<!-- * number_of_standard_tRNAs_extracted -->
<!-- * tRNA_extraction_software -->
* completeness_approach:
  * marker gene based
  * reference genome based
  * other
*  bin_parameters:
  *  kmer+coverage (?)
*  binning_software
  *  metabat2 #get version
*  reassembly_post_binning
  *  no
*  mag_coverage_software
  *  bwa #get version

### Logging and redirecting stdout and stderr

Some scripts or tools write the output to stdout and others write outputs and errors while writing the output to a specified file. Logging will have to be done accordingly.

If the tool is writing the desired output to stdout, only the stderr shout be redirected to the log file. This can be done using `>2` for example, `command 2> log.txt 1> output.txt`

If the tool is writing the desired output to a file and stdout is only going to include messages for logging, stderr and stdout can both be redirected to the same log file using `&>` for example, `command &> log.txt`

Look [here](https://www.gnu.org/software/bash/manual/html_node/Redirections.html#:~:text=Redirection%20allows%20commands'%20file%20handles,the%20current%20shell%20execution%20environment.) for more information.
