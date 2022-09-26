library(tidyverse)

MAG_taxonomy_info <- read.csv("ELLEGARD_2018_2020-gtdbtk.bac120.summary.tsv", sep = "\t")
MAG_info <- read.csv("ELLEGARD_2018_2020-evaluate_summary.tsv", sep = "\t")
isolates <- read.csv("all_combined.csv")
MAG_clusters_info <- read.csv("Cdb.csv")
format_name <- function(genome){
  MAG = strsplit(genome, ".fa")[[1]]
  return(MAG)
}
MAG_clusters_info <- MAG_clusters_info %>% mutate(genome = Vectorize(format_name)(genome))
clusters <- MAG_clusters_info %>%
              select(genome, secondary_cluster)
genomes_info <- MAG_info %>%
                  select(Genome, Mean.Completeness, Mean.Contamination)
taxonomy <- MAG_taxonomy_info %>%
              select(user_genome, classification) %>%
                separate(classification, c("domain","phylum","class","order","family","genus","species"), sep = ";")

MAGs_collated <- left_join(clusters, taxonomy, by = c("genome"="user_genome"))
MAGs_collated <- left_join(genomes_info, MAGs_collated, by = c("Genome"="genome"))
write.csv(MAGs_collated, "all_info_MAGs.csv", row.names = F)

MAGs_collated_high <- MAGs_collated %>%
                        filter(Mean.Completeness > 95, Mean.Contamination < 5)
write.csv(MAGs_collated_high, "High_quality_all_info_MAGs.csv", row.names = F)

isolates <- isolates %>%
              filter(DB == "KE_2021_expandedGBR")

final_genome_info <- data.frame("ID" = isolates$Locus_tag,
                                 "Accession" = isolates$Genbank_acc,
                                 "Locus_tag" = isolates$Locus_tag,
                                 "Strain_name" = isolates$Strain_name,
                                 "Phylotype" = isolates$Phylotype,
                                 "SDP" = isolates$SDP,
                                 "Species" = isolates$Species,
                                 "Host" = isolates$Host,
                                 "Study" = isolates$Study,
                                 "Origin" = isolates$Origin,
                                 "Source_database" = isolates$DB,
                                 "Cluster" = isolates$SDP,
                                 "Group" = rep("EMPTY", length(isolates$Locus_tag)),
                                 "Genus" = rep(NA, length(isolates$Locus_tag)),
                                 "Family" = rep(NA, length(isolates$Locus_tag)),
                                 "Order" = rep(NA, length(isolates$Locus_tag)),
                                 "Class" = rep(NA, length(isolates$Locus_tag)))
MAGs_collated_high <- data.frame("ID" = MAGs_collated_high$Genome,
                                 "Accession" = rep(NA, length(MAGs_collated_high$Genome)),
                                 "Locus_tag" = rep(NA, length(MAGs_collated_high$Genome)),
                                 "Strain_name" = MAGs_collated_high$Genome,
                                 "Phylotype" = rep(NA, length(MAGs_collated_high$Genome)),
                                 "SDP" = rep(NA, length(MAGs_collated_high$Genome)),
                                 "Species" = MAGs_collated_high$species,
                                 "Host" = rep("EMPTY", length(MAGs_collated_high$Genome)),
                                 "Study" = rep(NA, length(MAGs_collated_high$Genome)),
                                 "Origin" = rep("EMPTY", length(MAGs_collated_high$Genome)),
                                 "Source_database" = rep("Shini_lab_MAGs", length(MAGs_collated_high$Genome)),
                                 "Cluster" = MAGs_collated_high$secondary_cluster,
                                 "Group" = rep("EMPTY", length(MAGs_collated_high$Genome)),
                                 "Genus" = MAGs_collated_high$genus,
                                 "Family" = MAGs_collated_high$family,
                                 "Order" = MAGs_collated_high$order,
                                 "Class" = MAGs_collated_high$class)
write.csv(rbind(final_genome_info, MAGs_collated_high), "all_GenomeInfo.csv", row.names = F)

# Host, Origin and Group are to be filled in manually!
