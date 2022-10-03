#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(dplyr)

script_usage <- function() {
    cat("\nUsage: Rscript make_phylo_table.R /directory/xxxx/ \n", fill=TRUE)
    cat("No argument provided for file prefix or path provided is not the absolute path. Script will exit.", fill=TRUE)
}

Args <- commandArgs(trailingOnly=TRUE)
if (length(Args) == 0 || Args[1] == "--help") {
  script_usage()
  quit()
}

if (length(Args) == 1) {
  project_dir_prefix <- Args[1]
  if ( !startsWith(project_dir_prefix, "/") ) {
    script_usage()
    quit()
  }
}

get_sample_name <- function(magname){
  paste0(head(strsplit(strsplit(magname, "MAG_")[[1]][2], "_")[[1]], -1), collapse="_")
}

get_host_name <- function(magname){
  sample_name = strsplit(magname, "MAG_")[[1]][2]
  if (grepl("Dr|Gr", sample_name)) {
    return("Apis mellifera")
  }
  if (grepl("Am", sample_name)) {
    return("Apis mellifera")
  }
  if (grepl("Ac", sample_name)) {
    return("Apis cerana")
  }
  if (grepl("M1.|M2.|M3.", sample_name)) {
    return("Apis mellifera")
  }
  if (grepl("C1.|C2.|C3.", sample_name)) {
    return("Apis cerana")
  }
  if (grepl("D1.|D2.|D3.", sample_name)) {
    return("Apis dorsata")
  }
  if (grepl("F1.|F2.|F3.", sample_name)) {
    return("Apis florea")
  }
  return("EMPTY")
}

get_origin_name <- function(magname){
  sample_name = strsplit(magname, "MAG_")[[1]][2]
  if (grepl("Dr|Gr", sample_name)) {
    return("Switzerland, Engel apiary")
  }
  if (grepl("Am", sample_name)) {
    return("Japan")
  }
  if (grepl("Ac", sample_name)) {
    return("Japan")
  }
  if (grepl("M1.|M2.|M3.", sample_name)) {
    return("India")
  }
  if (grepl("C1.|C2.|C3.", sample_name)) {
    return("India")
  }
  if (grepl("D1.|D2.|D3.", sample_name)) {
    return("India")
  }
  if (grepl("F1.|F2.|F3.", sample_name)) {
    return("India")
  }
  return("EMPTY")
}

phy_group_dict = c("firm4" = "g__Bombilactobacillus",
            "g__Bombilactobacillus_outgroup" = "g__Bombilactobacillus",
            "firm5" = "g__Lactobacillus",
            "lacto" = "g__Lactobacillus",
            "g__Lactobacillus_outgroup" = "g__Lactobacillus",
            "bifido" = "g__Bifidobacterium",
            "g__Bifidobacterium_outgroup" = "g__Bifidobacterium",
            "gilli" = "g__Gilliamella",
            "entero" = "g__Gilliamella",
            "g__Gilliamella_outgroup" = "g__Gilliamella",
            "fper" = "g__Frischella",
            "g__Frischella_outgroup" = "g__Frischella",
            "snod" = "g__Snodgrassella",
            "g__Snodgrassella_outgroup" = "g__Snodgrassella",
            "bapis" = "g__Bartonella",
            "g__Bartonella_outgroup" = "g__Bartonella",
            # "" = "g__Enterobacter",
            "g__Enterobacter_outgroup" = "g__Enterobacter",
            # "" = "g__",
            # "" = "g__Pectinatus",
            "g__Pectinatus_outgroup" = "g__Pectinatus",
            "api" = "g__Apibacter",
            "g__Apibacter_outgroup" = "g__Apibacter",
            # "" = "g__Dysgonomonas",
            "g__Dysgonomonas_outgroup" = "g__Dysgonomonas",
            # "" = "g__Spiroplasma",
            "g__Spiroplasma_outgroup" = "g__Spiroplasma",
            # "" = "g__Zymobacter",
            "g__Zymobacter_outgroup" = "g__Zymobacter",
            # "" = "g__Entomomonas",
            "g__Entomomonas_outgroup" = "g__Entomomonas",
            # "" = "g__Saezia",
            "g__Saezia_outgroup" = "g__Saezia",
            # "" = "g__Parolsenella",
            "g__Parolsenella_outgroup" = "g__Parolsenella",
            # "" = "g__WRHT01",
            "g__WRHT01_outgroup" = "g__WRHT01",
            "com" = "g__Commensalibacter",
            "g__Commensalibacter_outgroup" = "g__Commensalibacter",
            "lkun" = "g__Apilactobacillus",
            "g__Apilactobacillus_outgroup" = "g__Apilactobacillus",
            "bom" = "g__Bombella",
            "g__Bombella_outgroup" = "g__Bombella"
          )

get_group_phy <- function(phy){
  return(c(phy_group_dict[phy][[1]]))
}
# setwd("/Volumes/Storage/Work/Temp-From-NAS/cross-species-analysis-India")
MAG_taxonomy_info <- read.csv(paste0(project_dir_prefix, "/06_MAG_binning/gtdbtk_out_dir/classify/gtdbtk.bac120.summary.tsv"), sep = "\t")
MAG_info <- read.csv(paste0(project_dir_prefix, "/06_MAG_binning/all_genomes.csv"), sep = "\t")

isolates <- read.csv(paste0(project_dir_prefix, "config/IsolateGenomeInfo.csv"))

MAG_clusters_info <- read.csv(paste0(project_dir_prefix, "/06_MAG_binning/drep_results/data_tables/Cdb.csv"))
format_name <- function(genome){
  MAG = strsplit(genome, ".fa")[[1]]
  return(MAG)
}
MAG_clusters_info <- MAG_clusters_info %>% mutate(genome = Vectorize(format_name)(genome))
clusters <- MAG_clusters_info %>%
              select(genome, secondary_cluster)
genomes_info <- MAG_info %>%
                  select(Bin.Id, Completeness, Contamination)
colnames(genomes_info) = c("Genome", "Completeness", "Contamination")
taxonomy <- MAG_taxonomy_info %>%
              select(user_genome, classification) %>%
                separate(classification, c("domain","phylum","class","order","family","genus","species"), sep = ";")

MAGs_collated <- left_join(clusters, taxonomy, by = c("genome"="user_genome"))
MAGs_collated <- left_join(genomes_info, MAGs_collated, by = c("Genome"="genome"))

MAGs_collated_for_tree <- MAGs_collated %>%
                        filter(Completeness > 50, Contamination < 5)

final_genome_info <- data.frame("ID" = isolates$ID,
                                 "Accession" = isolates$Accession,
                                 "Locus_tag" = isolates$Locus_tag_KE_DB,
                                 "Strain_name" = isolates$Strain_name,
                                 "Phylotype" = isolates$Phylotype,
                                 "SDP" = isolates$SDP,
                                 "Species" = isolates$Species,
                                 "Host" = isolates$Host,
                                 "Study" = isolates$Study,
                                 "Origin" = isolates$Origin,
                                 "Source_database" = isolates$Source_database,
                                 "Cluster" = isolates$SDP,
                                 "Group" = isolates$Group,
                                 "Genus" = rep(NA, length(isolates$ID)),
                                 "Family" = rep(NA, length(isolates$ID)),
                                 "Order" = rep(NA, length(isolates$ID)),
                                 "Class" = rep(NA, length(isolates$ID)),
                                 "Sample" = rep(NA, length(isolates$ID)),
                                 "Group_auto" = rep("EMPTY", length(isolates$ID)))

final_genome_info <- final_genome_info %>%
                      mutate("Group_auto" = ifelse(is.na(Phylotype), Vectorize(get_group_phy)(Group), Vectorize(get_group_phy)(Phylotype)))
final_genome_info$Group_auto <- as.character(final_genome_info$Group_auto)

MAGs_collated_for_tree <- data.frame("ID" = MAGs_collated_for_tree$Genome,
                                 "Accession" = rep(NA, length(MAGs_collated_for_tree$Genome)),
                                 "Locus_tag" = rep(NA, length(MAGs_collated_for_tree$Genome)),
                                 "Strain_name" = MAGs_collated_for_tree$Genome,
                                 "Phylotype" = rep(NA, length(MAGs_collated_for_tree$Genome)),
                                 "SDP" = rep(NA, length(MAGs_collated_for_tree$Genome)),
                                 "Species" = MAGs_collated_for_tree$species,
                                 "Host" = rep("EMPTY", length(MAGs_collated_for_tree$Genome)),
                                 "Study" = rep(NA, length(MAGs_collated_for_tree$Genome)),
                                 "Origin" = rep("EMPTY", length(MAGs_collated_for_tree$Genome)),
                                 "Source_database" = rep("MAGs", length(MAGs_collated_for_tree$Genome)),
                                 "Cluster" = MAGs_collated_for_tree$secondary_cluster,
                                 "Group" = rep("EMPTY", length(MAGs_collated_for_tree$Genome)),
                                 "Genus" = MAGs_collated_for_tree$genus,
                                 "Family" = MAGs_collated_for_tree$family,
                                 "Order" = MAGs_collated_for_tree$order,
                                 "Class" = MAGs_collated_for_tree$class,
                                 "Sample" = rep("EMPTY", length(MAGs_collated_for_tree$Genome)))

MAGs_collated <- data.frame("ID" = MAGs_collated$Genome,
                                 "Accession" = rep(NA, length(MAGs_collated$Genome)),
                                 "Locus_tag" = rep(NA, length(MAGs_collated$Genome)),
                                 "Strain_name" = MAGs_collated$Genome,
                                 "Phylotype" = rep(NA, length(MAGs_collated$Genome)),
                                 "SDP" = rep(NA, length(MAGs_collated$Genome)),
                                 "Species" = MAGs_collated$species,
                                 "Host" = rep("EMPTY", length(MAGs_collated$Genome)),
                                 "Study" = rep(NA, length(MAGs_collated$Genome)),
                                 "Origin" = rep("EMPTY", length(MAGs_collated$Genome)),
                                 "Source_database" = rep("MAGs", length(MAGs_collated$Genome)),
                                 "Cluster" = MAGs_collated$secondary_cluster,
                                 "Group" = rep("EMPTY", length(MAGs_collated$Genome)),
                                 "Genus" = MAGs_collated$genus,
                                 "Family" = MAGs_collated$family,
                                 "Order" = MAGs_collated$order,
                                 "Class" = MAGs_collated$class,
                                 "Sample" = rep("EMPTY", length(MAGs_collated$Genome)))

MAGs_collated <- MAGs_collated %>%
  mutate(Group_auto = Genus) %>%
    mutate(Host = Vectorize(get_host_name)(ID)) %>%
      mutate(Origin = Vectorize(get_origin_name)(ID)) %>%
        mutate(Sample = Vectorize(get_sample_name)(ID))

MAGs_collated_for_tree <- MAGs_collated_for_tree %>%
  mutate("Group_auto" = Genus) %>%
    mutate("Host" = Vectorize(get_host_name)(ID)) %>%
      mutate(Origin = Vectorize(get_origin_name)(ID)) %>%
        mutate(Sample = Vectorize(get_sample_name)(ID))

write.csv(rbind(final_genome_info, MAGs_collated_for_tree), paste0(project_dir_prefix, "06_MAG_binning/ForTree_GenomeInfo_auto.csv"), row.names = F, quote = T)
write.csv(rbind(final_genome_info, MAGs_collated), paste0(project_dir_prefix, "06_MAG_binning/all_GenomeInfo_auto.csv"), row.names = F, quote = T)
write.csv(rbind(MAGs_collated), paste0(project_dir_prefix, "06_MAG_binning/MAGs_GenomeInfo_auto.csv"), row.names = F, quote = T)
