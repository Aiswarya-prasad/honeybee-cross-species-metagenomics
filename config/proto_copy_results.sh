# prefix="aprasad@curnagl.dcsr.unil.ch:/scratch/aprasad/211018_Medgenome_india_samples"
prefix="/scratch/aprasad/211018_Medgenome_india_samples"
# outdir="/Users/aiswarya/One*/Work/2-2-Projects/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/211018_Medgenome_india_samples"
outdir="/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/211018_Medgenome_india_samples_host_filtered"
# config / input files
mkdir -p "$outdir"
mkdir -p "$outdir"/config
rsync -r --progress -v $prefix/config/* "$outdir"/config
rsync -r --progress -v $prefix/*_Report.Rmd "$outdir"/
# results
mkdir -p "$outdir"/fastqc
mkdir -p "$outdir"/fastqc/raw
rsync -r --progress -v $prefix/fastqc/raw/*_*_fastqc.html "$outdir"/fastqc/raw
mkdir -p "$outdir"/fastqc/trim
rsync -r --progress -v $prefix/fastqc/trim/*_*_trim_fastqc.html "$outdir"/fastqc/trim
mkdir -p "$outdir"/02_HostMapping
rsync -r --progress -v $prefix/02_HostMapping/*_flagstat.tsv "$outdir"/02_HostMapping
rsync -r --progress -v $prefix/02_HostMapping/*_coverage.tsv "$outdir"/02_HostMapping
rsync -r --progress -v $prefix/02_HostMapping/*_coverage_histogram.txt "$outdir"/02_HostMapping
mkdir -p "$outdir"/03_MicrobiomeMapping
rsync -r --progress -v $prefix/03_MicrobiomeMapping/*_flagstat.tsv "$outdir"/03_MicrobiomeMapping
mkdir -p "$outdir"/04_MicrobiomeMappingDirect
rsync -r --progress -v $prefix/04_MicrobiomeMappingDirect/*_flagstat.tsv "$outdir"/04_MicrobiomeMappingDirect
mkdir -p "$outdir"/02_motus_profile
rsync -r --progress -v $prefix/02_motus_profile/samples_merged.motus "$outdir"/02_motus_profile
mkdir -p "$outdir"/05_Assembly/
mkdir -p "$outdir"/05_Assembly/MapToAssembly
rsync -r --progress -v $prefix/05_Assembly/MapToAssembly/Assembly_mapping_summary.csv "$outdir"/05_Assembly/MapToAssembly
mkdir -p "$outdir"/06_MAG_binning
mkdir -p "$outdir"/06_MAG_binning/contig_fates
# we want some depth files but not all! i.e. we only want sample mapped to sample
# TODO: parse sample list from config file
for sample in M1.1  M1.2  M1.3  M1.4  M1.5 C1.1  C1.2  C1.3  C1.4  C1.5 C2.1  C2.2  C2.3  C2.4  C2.5 C3.1  C3.2  C3.3  C3.4  C3.5 D1.1 D1.2 D1.3 D1.4 D1.5 D2.1 D2.2 D2.3 D2.4 D2.5 D3.1 D3.2 D3.3 D3.4 D3.5 F1.1 F1.2 F1.3 F1.4 F1.5 F2.1 F2.2 F2.3 F2.4 F2.5 F3.1 F3.2 F3.3 F3.4 F3.5;
do
  mkdir -p "$outdir"/06_MAG_binning/backmapping/${sample}
  rsync -r --progress -v $prefix/06_MAG_binning/backmapping/${sample}/${sample}_mapped_to_${sample}.depth "$outdir"/06_MAG_binning/backmapping/${sample}
done
rsync -r --progress -v $prefix/06_MAG_binning/contig_fates/*_contig_fates.csv "$outdir"/06_MAG_binning/contig_fates/
mkdir -p "$outdir"/06_MAG_binning/
rsync -r --progress -v $prefix/06_MAG_binning/*_auto.tsv "$outdir"/06_MAG_binning
mkdir -p "$outdir"/06_MAG_binning/contig_fates
mkdir -p "$outdir"/06_MAG_binning/contig_fates/backmapping_coverages/
rsync -r --progress -v $prefix/06_MAG_binning/contig_fates/backmapping_coverages/*_contig_coverages.csv "$outdir"/06_MAG_binning/contig_fates/backmapping_coverages/
# since trees are only made for specific "groups"
# TODO: parse "group" list from config file
for group in "g__Bombilactobacillus" "g__Lactobacillus" "g__Lactobacillus" "g__Bifidobacterium" "g__Gilliamella" "g__Gilliamella" "g__Frischella" "g__Snodgrassella" "g__Bartonella" "g__Enterobacter" "g__Pectinatus" "g__Apibacter" "g__Dysgonomonas" "g__Spiroplasma" "g__Entomomonas" "g__Saezia" "g__Parolsenella" "g__WRHT01" "g__Commensalibacter" "g__Apilactobacillus" "g__Bombella"
do
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Comparative_Genomics_Statistics
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Orthogroups
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies/05_IQTree/
  mkdir -p "$outdir"/07_AnnotationAndPhylogenies/05_IQTree/${group}/
  rsync -r --progress -v $prefix/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Orthogroups/Orthogroups.txt "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Orthogroups/Orthogroups.txt
  rsync -r --progress -v $prefix/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Orthogroups/Orthogroups.GeneCount.tsv "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Orthogroups/Orthogroups.GeneCount.tsv
  rsync -r --progress -v $prefix/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Comparative_Genomics_Statistics/Statistics_Overall.tsv "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}/OrthoFinder/Results_${group}/Comparative_Genomics_Statistics/Statistics_Overall.tsv
  rsync -r --progress -v $prefix/07_AnnotationAndPhylogenies/02_orthofinder/${group}_Orthogroups_summary.csv "$outdir"/07_AnnotationAndPhylogenies/02_orthofinder/${group}_Orthogroups_summary.csv
  rsync -r --progress -v $prefix/07_AnnotationAndPhylogenies/05_IQTree/${group}/${group}_Phylogeny.treefile "$outdir"/07_AnnotationAndPhylogenies/05_IQTree/${group}
  rsync -r --progress -v $prefix/07_AnnotationAndPhylogenies/05_IQTree/${group}/${group}_Phylogeny.contree "$outdir"/07_AnnotationAndPhylogenies/05_IQTree/${group}
done
mkdir -p "$outdir"/09_MagDatabaseProfiling
mkdir -p "$outdir"/09_MagDatabaseProfiling/MAGsDatabaseMapping
rsync -r --progress -v $prefix/09_MagDatabaseProfiling/MAGsDatabaseMapping/*_flagstat.tsv "$outdir"/09_MagDatabaseProfiling/MAGsDatabaseMapping
rsync -r --progress -v $prefix/09_MagDatabaseProfiling/MAGsDatabaseMapping/*_host-filtered_flagstat.tsv "$outdir"/09_MagDatabaseProfiling/MAGsDatabaseMapping
mkdir -p "$outdir"/database
mkdir -p "$outdir"/database/MAGs_database_Orthofinder
rsync -r --progress -v $prefix/database/MAGs_database_Orthofinder/*_Orthogroups_summary.csv "$outdir"/database/MAGs_database_Orthofinder
rsync -r --progress -v $prefix/database/MAGs_database_Orthofinder/*_Orthogroups_filtered_summary.csv "$outdir"/database/MAGs_database_Orthofinder
# use magic to get list of names of groups
for group in g__113_1 g__172_1 g__73_1 g__Apibacter g__Bartonella g__Bifidobacterium g__Bombella g__Bombilactobacillus g__Commensalibacter g__Dysgonomonas g__Entomomonas g__Frischella g__Gilliamella g__Lactobacillus g__Parolsenella g__Pectinatus g__Saezia g__Snodgrassella g__Spiroplasma g__WRHT01 g__Zymobacter;
do
  mkdir -p "$outdir"/database/MAGs_database_Orthofinder/${group}
  mkdir -p "$outdir"/database/MAGs_database_Orthofinder/${group}/OrthoFinder
  rsync -r --progress -v $prefix/database/MAGs_database_Orthofinder/${group}/OrthoFinder/${group}_single_ortho.txt "$outdir"/database/MAGs_database_Orthofinder/${group}/OrthoFinder
  rsync -r --progress -v $prefix/database/MAGs_database_Orthofinder/${group}/OrthoFinder/${group}_single_ortho_filt.txt "$outdir"/database/MAGs_database_Orthofinder/${group}/OrthoFinder
done
mkdir -p "$outdir"/09_MagDatabaseProfiling/
mkdir -p "$outdir"/09_MagDatabaseProfiling/CoverageEstimation
mkdir -p "$outdir"/09_MagDatabaseProfiling/CoverageEstimation/Merged
rsync -r --progress -v $prefix/09_MagDatabaseProfiling/CoverageEstimation/Merged/*_coord.pdf "$outdir"/09_MagDatabaseProfiling/CoverageEstimation/Merged
rsync -r --progress -v $prefix/09_MagDatabaseProfiling/CoverageEstimation/Merged/*.txt "$outdir"/09_MagDatabaseProfiling/CoverageEstimation/Merged
mkdir -p "$outdir"/10_instrain
for sample in M1.1  M1.2  M1.3  M1.4  M1.5 C1.1  C1.2  C1.3  C1.4  C1.5 C2.1  C2.2  C2.3  C2.4  C2.5 C3.1  C3.2  C3.3  C3.4  C3.5 D1.1 D1.2 D1.3 D1.4 D1.5 D2.1 D2.2 D2.3 D2.4 D2.5 D3.1 D3.2 D3.3 D3.4 D3.5 F1.1 F1.2 F1.3 F1.4 F1.5 F2.1 F2.2 F2.3 F2.4 F2.5 F3.1 F3.2 F3.3 F3.4 F3.5;
do
  mkdir -p $prefix/10_instrain/${sample}_profile.IS
  mkdir -p $prefix/10_instrain/${sample}_profile.IS/output
  mkdir -p $prefix/10_instrain/${sample}_profile.IS/log
  mkdir -p $prefix/10_instrain/${sample}_profile.IS/figures
  rm -rf $prefix/10_instrain/${sample}_profile.IS/raw_data
  rsync -r --progress -v $prefix/10_instrain/${sample}_profile.IS/output/* "$outdir"/10_instrain/${sample}_profile.IS/output
  rsync -r --progress -v $prefix/10_instrain/${sample}_profile.IS/log/* "$outdir"/10_instrain/${sample}_profile.IS/log
  rsync -r --progress -v $prefix/10_instrain/${sample}_profile.IS/figures/* "$outdir"/10_instrain/${sample}_profile.IS/figures
done
rsync -r --progress -v $prefix/10_instrain/rep_mags.IS.COMPARE "$outdir"/10_instrain
rsync -r --progress -v $prefix/*_Report.html "$outdir"
# other results not (currently) explicitly mentioned in snakefile ...
mkdir -p "$outdir"/06_MAG_binning/drep_results/
mkdir -p "$outdir"/06_MAG_binning/drep_results/data_tables
rsync -r --progress -v $prefix/06_MAG_binning/drep_results/data_tables/*.csv "$outdir"/06_MAG_binning/drep_results/data_tables
mkdir -p "$outdir"/06_MAG_binning/drep_results/figures
rsync -r --progress -v $prefix/06_MAG_binning/drep_results/figures/* "$outdir"/06_MAG_binning/drep_results/figures
mkdir -p "$outdir"/06_MAG_binning/gtdbtk_out_dir/
mkdir -p "$outdir"/06_MAG_binning/gtdbtk_out_dir/classify/
rsync -r --progress -v $prefix/06_MAG_binning/gtdbtk_out_dir/classify/*.summary.tsv "$outdir"/06_MAG_binning/gtdbtk_out_dir/classify/
mkdir -p "$outdir"/06_MAG_binning/evaluate_bins
rsync -r --progress -v $prefix/06_MAG_binning/evaluate_bins/*.summary "$outdir"/06_MAG_binning/evaluate_bins
# find a way to get "all" groups
mkdir -p 12_species_validation
rsync -r --progress $prefix/12_species_validation/*/*perc_id* "$outdir"/12_species_validation
rsync -r --progress $prefix/Figures/* "$outdir"/Figures

# parsed outputs - this script does not do it but this can be run on remote and the result can be copied
echo "benchmarks" > benchmarks_cat.txt; for file in logs/*.benchmark; do echo $file >> benchmarks_cat.txt; cat $file | cut -f2,3,4,5 >> benchmarks_cat.txt; done
mkdir -p "$outdir"/logs
rsync -r --progress -v $prefix/benchmarks_cat.txt "$outdir"/logs/benchmarks_cat.txt
# logs/backup.log - output from backup script can be copied if interested=
