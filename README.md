# Summary of contents

# Directory organization

The chunk below describes

```{bash}
├── [ 237G]  00_RawData
│   ├── [ 1.8G]  ${SAMPLE}_R1.fastq.gz
│   └── [ 1.8G]  ${SAMPLE}_R2.fastq.gz
├── [ 210G]  01_Trimmed
│   ├── [ 1.6G]  ${SAMPLE}_R1_trim.fastq.gz
│   └── [ 1.6G]  ${SAMPLE}_R2_trim.fastq.gz
├── [ 1.7T]  02_HostMapping
│   ├── [ 6.6K]  ${SAMPLE}_coverage_histogram.txt
│   ├── [  425]  ${SAMPLE}_coverage.tsv
│   ├── [  541]  ${SAMPLE}_flagstat.tsv
│   ├── [ 3.3G]  ${SAMPLE}_R1_host_unmapped.fastq
│   └── [ 3.2G]  ${SAMPLE}_R2_host_unmapped.fastq
├── [ 150G]  03_MicrobiomeMapping
│   ├── [  541]  ${SAMPLE}_flagstat.tsv
│   ├── [ 3.1G]  ${SAMPLE}_R1_microbiome_mapped.fastq
│   └── [ 3.1G]  ${SAMPLE}_R2_microbiome_mapped.fastq
├── [  30K]  04_MicrobiomeMappingDirect
│   └── [  542]  ${SAMPLE}_flagstat.tsv
├── [  62G]  05_Assembly
│   └── [  62G]  host_unmapped
│       ├── [  88M]  ${SAMPLE}_assembly_graph.fastg
│       ├── [  90M]  ${SAMPLE}_graph.fastg
│       ├── [  34M]  ${SAMPLE}_scaffolds.fasta
│       ├── [ 8.6K]  ${SAMPLE}_scaffolds.fasta.amb
│       ├── [ 223K]  ${SAMPLE}_scaffolds.fasta.ann
│       ├── [  34M]  ${SAMPLE}_scaffolds.fasta.bwt
│       ├── [ 8.4M]  ${SAMPLE}_scaffolds.fasta.pac
│       ├── [  17M]  ${SAMPLE}_scaffolds.fasta.sa
│       ├── [  42M]  ${SAMPLE}_scaffolds_unparsed.fasta
│       ├── [ 168K]  ${SAMPLE}_spades.log
│       ├── [  35G]  check_assembly
│       │   ├── [  670]  Assembly_mapped_counts.txt
│       │   ├── [ 3.5K]  Assembly_mapping_summary.csv
│       │   ├── [ 5.0G]  ${SAMPLE}_assembly.bam
│       │   ├── [   14]  ${SAMPLE}_assembly_mapped_counts.txt
│       │   ├── [  30G]  ${SAMPLE}_assembly.sam
│       ├── [  97M]  ${SAMPLE}_assembly_graph.fastg
│       ├── [  34M]  ${SAMPLE}_scaffolds.fasta
│       ├── [ 8.5K]  ${SAMPLE}_scaffolds.fasta.amb
│       ├── [ 329K]  ${SAMPLE}_scaffolds.fasta.ann
│       ├── [  34M]  ${SAMPLE}_scaffolds.fasta.bwt
│       ├── [ 8.4M]  ${SAMPLE}_scaffolds.fasta.pac
│       ├── [  17M]  ${SAMPLE}_scaffolds.fasta.sa
│       ├── [  45M]  ${SAMPLE}_scaffolds_unparsed.fasta
│       └── [ 174K]  ${SAMPLE}_spades.log
├── [  10G]  06_MAG_binning
│   ├── [ 181K]  all_GenomeInfo_auto.tsv
│   ├── [ 2.2G]  backmapping
│   │   └── [  18M]  M1.5
│   │       ├── [ 304K]  ${SAMPLE}_mapped_to_${SAMPLE}.depth
│   │       ├── [ 304K]  ${OTHER_SAMPLE}_mapped_to_${SAMPLE}.depth
│   │       └── [ 3.1M]  ${SAMPLE}_merged.depth
│   ├── [ 1.6G]  bins
│   │   └── [  22M]  ${SAMPLE}-metabat2
│   │       ├── [ 2.6M]  MAG.10.fa
│   │       └── [ 858K]  MAG.9.fa
│   ├── [ 1.6G]  bins_renamed
│   │   └── [  22M]  ${SAMPLE}
│   │       ├── [ 2.6M]  MAG_${SAMPLE}_10.fa
│   │       └── [ 858K]  MAG_${SAMPLE}_9.fa
│   ├── [ 2.5G]  contig_fates
│   │   ├── [ 1.9G]  backmapping_coverages
│   │   │   └── [  15M]  ${SAMPLE}_contig_coverages.csv
│   │   └── [ 2.1M]  ${SAMPLE}_contig_fates.csv
│   ├── [  30M]  drep_results
│   │   ├── [ 8.0K]  data
│   │   │   └── [ 4.0K]  Clustering_files
│   │   ├── [  30M]  data_tables
│   │   │   ├── [ 102K]  Bdb.csv
│   │   │   ├── [  33K]  genomeInfo.csv
│   │   │   ├── [  29M]  Mdb.csv
│   │   │   ├── [ 599K]  Ndb.csv
│   │   │   ├── [  28K]  Sdb.csv
│   │   │   └── [  20K]  Widb.csv
│   │   ├── [ 4.0K]  dereplicated_genomes
│   │   ├── [ 4.0K]  figures
│   │   └── [  92K]  log
│   │       └── [  88K]  logger.log
│   ├── [ 2.4G]  evaluate_bins
│   │   ├── [  38M]  ${SAMPLE}
│   │   │   ├── [  35M]  bins
│   │   │   │   ├── [ 3.8M]  MAG_${SAMPLE}_11
│   │   │   │   │   ├── [ 888K]  genes.faa
│   │   │   │   │   ├── [ 535K]  genes.gff
│   │   │   │   │   ├── [ 2.3M]  hmmer.analyze.txt
│   │   │   │   │   └── [  15K]  hmmer.tree.txt
│   │   │   ├── [ 2.7K]  checkm.log
│   │   │   ├── [1000K]  lineage.ms
│   │   │   └── [ 1.7M]  storage
│   │   │       ├── [  72K]  aai_qa
│   │   │       │   └── [ 4.7K]  MAG_${SAMPLE}_9
│   │   │       │       └── [  695]  PF01121.15.masked.faa
│   │   │       ├── [ 6.7K]  bin_stats.analyze.tsv
│   │   │       ├── [  54K]  bin_stats_ext.tsv
│   │   │       ├── [ 6.7K]  bin_stats.tree.tsv
│   │   │       ├── [ 621K]  checkm_hmm_info.pkl.gz
│   │   │       ├── [ 184K]  marker_gene_stats.tsv
│   │   │       ├── [ 4.6K]  phylo_hmm_info.pkl.gz
│   │   │       └── [ 761K]  tree
│   │   │           ├── [  96K]  concatenated.fasta
│   │   │           ├── [ 327K]  concatenated.pplacer.json
│   │   │           ├── [ 265K]  concatenated.tre
│   │   │           ├── [  699]  PF00164.20.masked.faa
│   │   │           ├── [ 1.0K]  pplacer.out
│   │   └── [ 1.3K]  ${SAMPLE}_checkm.summary
│   ├── [ 126K]  ForTree_GenomeInfo_auto.tsv
│   └── [ 2.0M]  gtdbtk_out_dir
│       └── [ 2.0M]  classify
│           └── [ 2.0M]  gtdbtk.bac120.summary.tsv
# examples of ${GENOME}: ESL0304.fa, MAG_M1.3_4.fa
├── [  24G]  07_Phylogenies
│   ├── [ 1.4G]  00_genomes
│   │   ├── [ 2.4M]  ${GENOME}.fa
│   ├── [  15G]  01_prokka
│   │   ├── [  26M]  ${GENOME}
│   │   │   ├── [ 203K]  ${GENOME}.err
│   │   │   ├── [ 766K]  ${GENOME}.faa
│   │   │   ├── [ 2.1M]  ${GENOME}.ffn
│   │   │   ├── [ 2.4M]  ${GENOME}.fna
│   │   │   ├── [ 2.4M]  ${GENOME}.fsa
│   │   │   ├── [ 5.1M]  ${GENOME}.gbk
│   │   │   ├── [ 3.2M]  ${GENOME}.gff
│   │   │   ├── [  32K]  ${GENOME}.log
│   │   │   ├── [ 8.7M]  ${GENOME}.sqn
│   │   │   ├── [ 621K]  ${GENOME}.tbl
│   │   │   ├── [ 200K]  ${GENOME}.tsv
│   │   │   └── [  105]  ${GENOME}.txt
# examples of ${GROUP}: g__Lactobacillus, g__Apibacter
│   ├── [ 6.9G]  02_orthofinder
│   │   ├── [ 249M]  ${GROUP}
│   │   │   ├── [ 805K]  ${GENOME}.faa
│   │   │   ├── [ 177M]  OrthoFinder
│   │   │   │   ├── [ 215K]  ${GROUP}_single_ortho_MAGs.txt
│   │   │   │   ├── [ 587K]  ${GROUP}_single_ortho.txt
│   │   │   │   └── [ 176M]  Results_${GROUP}
│   │   │   │       ├── [ 2.5K]  Citation.txt
│   │   │   │       ├── [  25K]  Comparative_Genomics_Statistics
│   │   │   │       │   ├── [ 6.0K]  Orthogroups_SpeciesOverlaps.tsv
│   │   │   │       │   ├── [ 1.6K]  Statistics_Overall.tsv
│   │   │   │       │   └── [  14K]  Statistics_PerSpecies.tsv
│   │   │   │       ├── [ 1.2K]  Log.txt
│   │   │   │       ├── [ 2.6M]  Orthogroups
│   │   │   │       │   ├── [ 261K]  Orthogroups.GeneCount.tsv
│   │   │   │       │   ├── [ 3.2K]  Orthogroups_SingleCopyOrthologues.txt
│   │   │   │       │   ├── [ 1.1M]  Orthogroups.tsv
│   │   │   │       │   ├── [ 1.1M]  Orthogroups.txt
│   │   │   │       │   └── [  79K]  Orthogroups_UnassignedGenes.tsv
│   │   │   │       └── [  23M]  Orthogroup_Sequences
│   │   │   │           └── [  13K]  OG0000184.fa
│   │   │   └── [    0]  single_ortholog_sequences.done
│   │   └── [ 147K]  ${GROUP}_Orthogroups_summary.csv
│   ├── [ 506M]  03_aligned_orthogroups
│   │   ├── [  21M]  ${GROUP}
│   │   │   ├── [  10M]  CoreGeneAlignment.fasta
│   │   │   ├── [  31K]  OG0000184_aligned_pruned.fasta
│   └── [ 1.8M]  05_IQTree
│       ├── [ 325K]  ${GROUP}
│       │   ├── [ 4.1K]  ${GROUP}_Phylogeny.bionj
│       │   ├── [  92K]  ${GROUP}_Phylogeny.ckp.gz
│       │   ├── [ 4.3K]  ${GROUP}_Phylogeny.contree
│       │   ├── [  48K]  ${GROUP}_Phylogeny.iqtree
│       │   ├── [  30K]  ${GROUP}_Phylogeny.log
│       │   ├── [ 122K]  ${GROUP}_Phylogeny.mldist
│       │   ├── [ 5.2K]  ${GROUP}_Phylogeny.model.gz
│       │   ├── [  11K]  ${GROUP}_Phylogeny.splits.nex
│       │   └── [ 4.3K]  ${GROUP}_Phylogeny.treefile
├── [  23M]  08_motus_profile
│   └── [  23M]  samples_merged.motus
├── [  59G]  09_MapToAssembly
│   └── [ 1.0G]  ${SAMPLE}.bam
├── [ 172K]  211018_Medgenome_india_samples_Report.Rmd
├── [ 419M]  assembly_summary_genbank.txt
├── [ 155K]  config
│   ├── [ 4.8K]  config.yaml
│   ├── [  34K]  initialize_project-211018_Medgenome_india_samples.sh
│   ├── [  32K]  IsolateGenomeInfo.csv
│   ├── [  31K]  Metadata_211018_Medgenome_india_samples.csv
│   └── [  51K]  reinitialize_project-211018_Medgenome_india_samples.sh
├── [ 3.5G]  database
│   ├── [ 869M]  4_host_db
│   ├── [ 1.1M]  4_host_db.amb
│   ├── [  298]  4_host_db.ann
│   ├── [ 869M]  4_host_db.bwt
│   ├── [  163]  4_host_db_list.txt
│   ├── [ 217M]  4_host_db.pac
│   ├── [ 435M]  4_host_db.sa
│   ├── [ 448M]  genome_db_220315_AP
│   ├── [10.0K]  genome_db_220315_AP.amb
│   ├── [ 7.0K]  genome_db_220315_AP.ann
│   ├── [ 441M]  genome_db_220315_AP.bwt
│   ├── [ 6.0K]  genome_db_220315_AP.fai
│   ├── [ 1.4K]  genome_db_220315_AP_list.txt
│   ├── [ 8.8K]  genome_db_220315_AP_metafile.txt
│   ├── [ 110M]  genome_db_220315_AP.pac
│   ├── [ 220M]  genome_db_220315_AP.sa
├── [ 6.2K]  envs
│   ├── [  218]  core-cov-env.yaml
│   ├── [   95]  gtdb-env.yaml
│   ├── [  135]  mags-env.yaml
│   ├── [  249]  mapping-env.yaml
│   ├── [   68]  motus-env.yaml
│   ├── [  179]  phylogenies-env.yaml
│   ├── [  274]  popcogent-env.yaml
│   ├── [  477]  rmd-env.yaml
│   ├── [  153]  scripts-env.yaml
│   ├── [  127]  snv-env.yaml
│   ├── [  162]  spades-env.yaml
│   └── [  129]  trim-qc-env.yaml
├── [ 275M]  fastqc
│   ├── [ 130M]  raw
│   │   ├── [ 562K]  ${SAMPLE}_R1_fastqc.html
│   │   ├── [ 765K]  ${SAMPLE}_R1_fastqc.zip
│   │   ├── [ 566K]  ${SAMPLE}_R2_fastqc.html
│   │   └── [ 772K]  ${SAMPLE}_R2_fastqc.zip
│   └── [ 145M]  trim
│       ├── [ 602K]  ${SAMPLE}_R1_trim_fastqc.html
│       ├── [ 753K]  ${SAMPLE}_R1_trim_fastqc.zip
│       ├── [ 606K]  ${SAMPLE}_R2_trim_fastqc.html
│       └── [ 758K]  ${SAMPLE}_R2_trim_fastqc.zip
├── [  15M]  Figures
│   ├── [ 6.0K]  00a-Total_reads_trimming.pdf
│   ├── [ 5.6K]  00b-Total_reads_species_after_trimming.pdf
│   ├── [ 6.6K]  00c-Mapped_Unmapped_reads_prop.pdf
├── [ 620M]  logs
├── [ 172K]  README.md
├── [ 190K]  scripts
│   ├── [ 1.7K]  aln_aa_to_dna_phylo.py
│   ├── [ 1.7K]  aln_aa_to_dna.py
│   ├── [  795]  aln_calc.sh
│   ├── [  878]  calc_assembly_length.pl
│   ├── [ 1.9K]  calc_filt_core_length.py
│   ├── [ 3.7K]  calc_jaccard.pl
│   ├── [ 4.3K]  calc_perc_id_orthologs_phylo.py
│   ├── [ 3.9K]  calc_shared_SNV_fraction.pl
│   ├── [  26K]  copy_faa_files.sh
│   ├── [ 1.3K]  copy_raw_data.sh
│   ├── [ 9.5K]  core_cov.py
│   ├── [ 8.1K]  core_cov.R
│   ├── [ 2.0K]  corecov_split_by_sdp.py
│   ├── [  749]  csv_to_tsv.py
│   ├── [ 3.7K]  cum_curve_SNVs_host.pl
│   ├── [ 8.7K]  cum_curve_SNVs_host.py
│   ├── [ 1.3K]  distance_matrix.R
│   ├── [ 4.5K]  download_data.py
│   ├── [ 3.9K]  extract_orthologs_phylo.py
│   ├── [ 2.7K]  fasta_generate_regions.py
│   ├── [ 3.0K]  filt_core_bed.py
│   ├── [ 2.1K]  filter_bam.py
│   ├── [ 2.5K]  filter_orthologs_phylo.py
│   ├── [ 1.8K]  filter_orthologs.py
│   ├── [  965]  filter_sam_aln_length.pl
│   ├── [ 1.2K]  filter_sam_aln_length_unmapped.pl
│   ├── [ 3.8K]  filter_snvs.pl
│   ├── [ 4.3K]  filter_snvs.py
│   ├── [ 3.0K]  filt_vcf_samples.pl
│   ├── [ 1.4K]  filt_vcf_samples.py
│   ├── [ 1.4K]  filt_vcf_samples.py
│   ├── [ 3.6K]  get_single_ortho_phylo.py
│   ├── [ 1.5K]  get_single_ortho.py
│   ├── [ 2.1K]  KEs_trim_aln.py
│   ├── [ 4.0K]  Make_metadata.R
│   ├── [  10K]  make_phylo_table.R
│   ├── [ 1.4K]  merge_depths.pl
│   ├── [ 8.2K]  my_core_cov.R
│   ├── [ 2.8K]  parse_core_cov.py
│   ├── [ 1.1K]  parse_spades_metagenome.pl
│   ├── [ 1.6K]  parse_spades_metagenome.py
│   ├── [ 7.5K]  prune_and_concat_alns.py
│   ├── [  976]  rearange_faa.py
│   ├── [ 3.8K]  subset_metagenome_db.py
│   ├── [ 3.3K]  subset_orthofile.py
│   ├── [  973]  summarise_filter_snvs.py
│   ├── [ 6.6K]  summarise_orthogroups.py
│   ├── [ 4.7K]  summarise_snps.pl
│   ├── [ 2.1K]  trim_aln_phylo.py
│   ├── [ 2.1K]  trim_aln.py
│   ├── [ 1.2K]  vcf_split_by_sdp.py
│   └── [ 1011]  write_adapters.py
├── [  75K]  Snakefile
```

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

The raw reads and trimmed reads were checked for quality using the tool [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). The report (html format) also summarises basic statistics including number of reads. The QC results can be found in their respective folders at `./fastqc/raw/{SAMPLE}_R*_fastqc.html` and `./fastqc/trim/{SAMPLE}_R*_trim_fastqc.html`.

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
  20160415_OBIWAN225	20160415	Kirsten_Ellegaard	12	GTF	Illumina	100	PE	HiSeq 2500	Genomic diversity landscape of the honey bee gut microbiota (2019, NatCom)	Foragers/Winterbees, Year 1, Les Droites \
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

conda environments are all specified in `envs/` and built by snakemake under various names in `/work/FAC/FBM/DMF/pengel/spirit/aprasad/Miniconda3/spirit_envs`

Run the pipeline in the conda environment called `snakmake_with_samtools` in the cluster. It is a clone of the snakemake environment made as recommended by Snakemake [docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba) followed by `conda install biopython` and later `conda install samtools` in it. This is so that Kirsten's core_cov script works (specific conda environments can only be specified for rules using bash).

## Description of directory structure

Directory names are largely self-explanatory.

>`./00_rawdata`, `./01_Trimmed`, `./02_HostMapping`, `./03_MicrobiomeMapping`
>`database` contains databases to be used for mapping. It also contains `./Orthofinder` files. These are described later in the sections describing associated rules.
>`./envs` contains all yaml files required for this pipeline. They contain a list of packages needed to specify the conda environment for various rules to work within.
>`./logs` contains log files
>`./scripts` contains all scripts needed for the snakemake pipeline. Many of these scripts are adapted from Kirsten's scripts from the zenodo directories, github or from the lab_resources directories.
The **results** of the core coverage estimation are stored in,
> `./04_CoreCov_211018_Medgenome_india_samples`
> `./07_SNVProfiling` is not fully implemented (yet) for these samples as it is not relevant at this time.
>`./fastqc` contains fastqc **results** for trimmed and raw files
+ bamfile_list_red.txt - required by KE's core coverage pipeline
+ bamfile_list.txt - required by KE's core coverage pipeline
+ Adapters-PE.fa - is generated based on index sequences by the script `./scripts/write_adapters.py` (was deleted earlier as it was on scratch. Needs to be re-written.)
+ config.yaml - comprises information including list of samples
+ index_table.csv - used by the script `./scripts/write_adapters.py` to make indexed adapters
+ Mapping_summary.csv - result from the rule summarize_mapping
+ rulegraph.pdf - summary DAG of rules in the pipeline (made using `snakemake --forceall --rulegraph | dot -Tpdf > Figuers/rulegraph.pdf`)
+ Report.Rmd - this report !
+ Report.html - this report compiled !
+ Snakefile - the pipipeline !!!

## Rules

![DAG of all rules in the pipeline](Figures/rulegraph.pdf){width=65%}

+ `rule raw_qc`
  - This rule runs [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on raw files and saves the output in `./fastqc/raw`.
+ `rule make_adapters`
  - This rule uses the script `_./scripts/write_adapters.py`_ which was deleted earlier.
  - It uses the index_table.csv files to make the Adapters-PE.fa file containing indexed adapters.
+ `rule trim`
  - This rules trims reads using [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).
  - The _Adapters-PE.fa_ files is used.
  - The trimming parameters are.
+ `rule trim_qc`
  - This rule runs [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on trimmed files and saves the output in `./fastqc/trim`.
+ `rule index_bwa`
  - Indexes genomes in `./database` for use by [bwa](http://bio-bwa.sourceforge.net/) using [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml#3).
+ `rule index_samtools`
  - Indexes genomes in `./database` for use by [samtools](http://www.htslib.org/doc/#manual-pages).
+ `rule make_genome_list`
  - Creats a text file corresponding to each set of genomes in `./database` to be used when we need to know which genomes are present in given genome database.
+ `rule host_mapping`
  - Uses [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml#:~:text=BWA%20is%20a%20software%20package,such%20as%20the%20human%20genome.&text=BWA%2DMEM%20also%20has%20better,genome%20(the%20index%20command)) to map reads for each sample to a database containing host genomes, `./database/4_host_db`.
  - Unmapped reads identified by samtools with the option `-f4` are stored in a seperate bam file.
  - The bam file with all alignments is used later by the counting rule and then deleted after counting.
+ `rule host_mapping_extract_reads`
  - Reads that did not map to the host database are extracted and then mapped to the microbiome database.
  - They are extracted using [picard](https://broadinstitute.github.io/picard/).
  - The option `-Xmx8g` ensures that java is given 8 GB memory. If suffecient memory is not allocated, the job will fail.
+ `rule host_mapping_count`
  - Counts the number of mapped, properly and unmapped reads from host mapping.
  - It uses the following flags to identify each kind of read:
    - **count number of properly mapped reads: `-f 67 -F 2304`**
        - 67 (include -f) flags
            + read paired (0x1)
            + read mapped in proper pair (0x2)
            + first in pair (0x40)
        - 2308 (exclude -F) flags
            + read unmapped (0x4)
            + supplementary alignment (0x800)
            + not primary (0x100)
    - **count number of mapped reads: `-f 67 -F 2304`**
        - 65 (include -f) flags
            + read paired (0x1)
            + first in pair (0x40)
        - 2308 (exclude -F) flags
            + read unmapped (0x4)
            + supplementary alignment (0x800)
            + not primary (0x100)
    - **count number of unmapped reads: `-f 67 -F 2304`**
        - 69 (include -f) flags
            + read paired (0x1)
            + read unmapped (0x4)
            + first in pair (0x40)
        - 2304 (exclude -F) flags
            + supplementary alignment (0x800)
            + not primary (0x100)
  - After counting, the bam file is deleted and an empty file is touched to mark that counting is complete for said file.
+ `rule microbiomedb_mapping`
  - The host unmapped reads extracted earlier are mapped to the microbiome database.
  - Mapped reads are extracted using a perls script as follows. First, unmapped reads are excluded using `-F4` and then supplementary reads are excluded `-F 0x800`. Finally, the remaining reads are sent through `./scripts/filter_sam_aln_length.pl`. The script filters away reads that have less than 50bps matching in the alignment.
+ `rule microbiomedb_extract_reads`
  - Extracts mapped reads identified as mentioned in the previous rule and saves them as fastq files.
+ `rule microbiome_mapping_count`
  - Counts reads as explained in the other counting rule, _host_mapping_count_.
+ `rule cat_and_clean_counts`
  - Compiles all the counts into 1 file for easier parsing by the summarize rule.
+ `rule summarize_mapping`
  - Summarizes counts in a csv file using the results of earlier rules and by counting fastq files.
+ `rule run_orthofinder`
  - Runs [Orthofinder](https://github.com/davidemms/OrthoFinder) for each phylotype.
  - **Before** running this, group genomes by phylotype in directories for Orthofinder to be able to get which groups to consider together. When the genomes for the database are downloaded at `./database/faa_files/{genome}`, they are all in one directory. Grouping was done using `./scripts/rearange_faa.py`. As written, it is to be run from the scripts directory in which it resides (!! it uses relative paths !!).
  - faa files for each genome comes from the respective databese (NCBI for example)
  - When orthofinder finishes, the following file will be generated and used for the following steps, `./database/faa_files/{phylotype}/OrthoFinder/Results_dir/Orthogroups/Orthogroups.txt`.
  - The file _Orthogroups.txt_ contains a list of orthogroups. Eg, each line would look like
    - **OG0000003**: C4S76_01365 C4S76_01370 C4S76_01375 C4S77_06100 C4S77_06130 C4S77_06135 C4S77_06775 C4S77_06780 C4S77_06785 C4S77_06790 C4S77_06795 C4S77_06800 C4S77_06805 C4S77_06810 C4S77_09595 C4S77_09600 C4S77_09605 C4S77_09610 **C4S77_09615 C4S77_09620 C4S77_10540 Ga0307799_111506**
    - where, **OG0000003** is an orthogroup for this group of genomes (phylotype) and **C4S77**, **Ga0307799** etc. are genomes that belong to that group. **09615, 09620, 10540** are genes from **C4S77** and **111506** from **Ga0307799** that belong to orthogroup OG0000003.
+ `rule get_single_ortho`
  - The files _Orthogroups.txt_ is parsed by `./scripts/get_single_ortho.py` and single-copy orthologs are written to `./database/Orthofinder/{phylotype}_single_ortho.txt`
  - The script reads each orthogroup and counts the number of genomes present in genes of that orthogroup. If the number of genes in the orthogroup and the number of genomes in the orthogroup are the same as the total number of genomes in the database for said phylotype, the genes in the group are considered single-copy core genes and included for core coverage estimation.
+ `rule extract_orthologs`
  - This rule prepares files with sequences of orthologs in order to calculate percentage identity (perc_id).
  - First, it reads the file `./database/Orthofinder/{phylotype}_single_ortho.txt` and gets all the genome-ids present in the ortholog-file, and all the gene-ids associated with each gene-family. Using this list it extracts and stores the sequences of each of the genes of an orthogroup in an faa file and ffn file corresponding to each group in the directory`./database/Orthofinder/{phylotype}_ortho_sequences/`.
  EXAMPLE:
  - `cat ./database/Orthofinder/firm5_ortho_sequences/OG0001034.faa` \
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
  - `cat ./database/Orthofinder/firm5_ortho_sequences/OG0001034.ffa` \
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
  - this rule relies on various tools and scripts tied together by `./scripts/aln_calc.sh`. The scripts are:
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
      + Uses as input, trimmed aligned sequences and a metafile (`./database/genome_db_210402_metafile.txt`) which is a tab-delim file with genome-id in tab1 and SDP-affiliation in tab 3
      + First, it checks the number of SDPs contained within the alignment. If more than one, it continues by calculating alignment percentage identity stats across SDPs. If only one SDP, exits script.
      + Next, it Compares the genomes in each SDP to all other genomes in the alignment: calculates percentage identity for all pairwise combinations. Calculates the max, min, and mean values, prints to file `./database/Orthofinder/{phylotype}_perc_id.txt` showing one orthogroup per line.
      EXAMPLE:
      + `cat ./database/Orthofinder/firm5_perc_id.txt` \
        ... \
        OG0001034	0.674	0.586	0.972 \
        ... \

+ `rule filter_orthogroups`
  - orthogroups are filtered based on:
    - Minimum gene-length 300bp (applied to all members of each gene-family)
    - Inter-SDP max alignment identity 95% (only if the phylotype contain multiple SDPs)
  - Short genes are filtered off, because they are likley to be less reliable for accurate coverage estimates. Similarly, the inter-SDP similarity threshold is used to ensure that there is enough divergence between the SDPs for reliable mapping (at least as estimated from the currently availabe genomes). It is worthwhile to check the number of gene-families before/after filtering. If a lot of gene-families were filtered off, this could be an indication that the SDPs are not properly discrete.
  - This finally results in the single-copy core genes that have been filtered to be used for core coverage estimation. `./database/Orthofinder/{phylotype}_single_ortho_filt.txt`.
+ `rule core_cov`
  - takes as input bam files with the alignments for each sample to be considered (as a text file containing a list of these files) and the _\_single_ortho_filt_ file. Outputs are written to `./04_CoreCov_"+ProjectIdentifier+"/{phylotype}_corecov.txt`.
  - The script reads the filtered orthofile `./database/Orthofinder/{phylotype}_single_ortho_filt.txt` and gets the gene-famililes and genome-ids for each SDP.
  - Then, from bed files, it finds the start and end positions of each of the genes in an orthogroup for each of the genomes of the orthogroup. It writes these to the file, `./04_CoreCov_*/{phylotype}_corecov.txt` each SDP in the phylotype, start position for each gene family in the genome marked reference for that SDP.
  - The coverage is also written to this file.
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
  - SDP abundance is estimated based on mapped read coverage of core genes. It sums up gene coverages of all the genes og OG families associated with said SDP across genomes belogining to the SDP.
  - It also reports PTR (Peak-Trough Ratio).
  - Most species in the database are represented by multiple genomes (< 98.5% gANI between genomes). Core genes are inferred at the phylotype. More accurate estimates can be obtained by using a large number (+700) of core genes.
+ `rule core_cov_plots`
  - This R-script will estimate the coverage at the terminus, using the summed core gene family coverages. If the cov-ter cannot be properly estimated (fx. due to draft genome status or lack of replication), an estimate will be generated using the median coverage across core gene families, and the PTR is set to NA. If more than 20\% of the core gene families have no coverage, the abundance will be set to zero. As output, a tabular file is generated (including the cov-ter/median cov, and PTR), and a pdf-file with plots for visual validation.
  - First, filter for samples with coverage of at least 1 on > 80% of the core genes. Next, values that are deviating no more than 2 times the median are kept others are discarded as outliers.
  - Next, gets fitted coordinates and append values to coord-table. It does this by using the segmented package. As explained below,

```{r eval=FALSE, echo=TRUE}
    x <- data_filt$Ref_pos
 	  y <- data_filt$Coverage
 	  psi_est <- max(x)/2
    lin.mod = y ~ x
    segmented(lin.mod, seg.Z=~x, psi=psi_est)
```
where, `lin.mod` is a simple linear model that was made by base R. The R package `Segmented` supports breakpoint analysis. The methods used by this package are applicable when segments are (nearly continuous) so this means that for the regression to make sense the core gene families selected should cover the reference genome well and without too many huge gaps. `psi`, is a starting value of the breakpoint. Example of a model fit using segemented,
```{r eval=FALSE, echo=TRUE}
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
```{r eval=FALSE, echo=TRUE}
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
```{r eval=FALSE, echo=TRUE}
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

  - the PTR was set to `NA`, the median will be plotted and used for quantification. Else, the segmented regression line is plotted, and the terminus coverage is used for quantification.

+ rule assemble_host_unmapped
  - Takes as input the R1 and R2 reads that were not mapped to the host and assembles them using [spades](https://cab.spbu.ru/software/meta-spades/) with the `--meta` tag and default parameters.
  - Memory allocation is not obvious. More documentation on this soon.
+ rule map_to_assembly
  - Map reads that were assembled against the contigs that they were assembled into using [bwa mem](http://bio-bwa.sourceforge.net).
+ rule cat_and_clean_counts_assembly
  - Compiles counts into one file for summarizing
+ rule summarize_mapping_assembly
  - similar to earlier rule "summarize_mapping"
+ rule backmapping
  - NxN mapping for
+ rule merge_depths
+ rule binning
+ rule process_metabat2
+ rule checkm_evaluation
+ rule prepare_info_for_drep
+ rule drep
+ rule gtdb_annotate
+ rule compile_report
+ rule backup
<!-- + `rule assemble_host_unmapped`
+ `rule mapping_red_db`
+ `rule subset_ortho_and_meta`
+ `rule core_cov_red`
+ `rule core_cov_red_plots`
+ `rule parse_core_cov_red`
+ `rule de_duplicate`
+ `rule freebayes_profiling`
+ `rule vcf_summary_stats`
+ `rule vcf_filtering1`
+ `rule vcf_filtering2`
+ `rule vcf_filtering3` -->
+ `rule onsuccess`

## Scripts

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
+ `./scripts/write_adapters.py`

## Envs


`core-cov-env.yaml`

```{yaml echo=TRUE}
name: core-cov-env

channels:

  - bioconda
  - conda-forge

dependencies:
  - python=3.*
  - samtools
  - bwa
  - mafft
  - orthofinder
  - biopython
  - r-base=3.5.1
  - r-plyr=1.8.6
  - r-segmented=1.1_0
  - r-cairo
```

`mapping-env.yaml`

```{yaml echo=TRUE}
name: mapping-env

channels:
  - bioconda
  - conda-forge
  - hcc

dependencies:
  - python=3.9.7
  - openjdk=11.0.9.1
  - perl=5.32.1=0_h7f98852_perl5
  - samtools=1.13=h8c37831_0
  - picard=2.26.2=hdfd78af_0
  - bwa=0.7.17=h5bf99c6_8
```

`rmd-env.yaml`

```{yaml echo=TRUE}
name: rmd-env

channels:
  - bioconda
  - conda-forge

dependencies:
  - r-base
  - r-rmarkdown
  - r-ggplot2
  - r-kableExtra
  - r-codetools
  - r-tidyverse
  - r-prettydoc
  - r-viridis
  - r-hrbrthemes
  - r-ggthemes
  - r-RColorBrewer
  - r-scales
  - r-segmented
  - r-shiny
  - r-dplyr
  - r-xlsx
  - r-DT
  - r-leaflet
```

`snv-env.yaml`

```{yaml echo=TRUE}
name: snv-env

channels:
  - bioconda
  - conda-forge
  - hcc

dependencies:
  - freebayes
  - vcftools
  - vcflib
```

`trim-qc-env.yaml`

```{yaml echo=TRUE}
name: trim-qc-env

channels:
  - bioconda
  - conda-forge
  - hcc

dependencies:
  - python
  - trimmomatic
  - fastqc
  - quast
```

# Next steps

It is clear that the database is not best suited for some SDPs found especially in host species other than _Apis mellifera_. So, the next step would be to implement a MAG based analysis to compare these samples. However, as the database was already shown to be well-suited for _Apis mellifera_ and _Apis cerana_, another set of analysis would compare these samples from the [publication](https://www.sciencedirect.com/science/article/pii/S0960982220305868) ([zenodo](https://zenodo.org/record/3747314#.YcGkvRPMK3I)) with the samples from India.

## Choosing representative genomes

Information in [drep](https://drep.readthedocs.io/en/latest/choosing_parameters.html) documentation.

$$A * Completeness - B * Contamination + C * (Contamination * frac{strainheterogeneity}{100}) + D * log(N50) + E * log(size) + F * (centrality - S_{a}ni) $$s
