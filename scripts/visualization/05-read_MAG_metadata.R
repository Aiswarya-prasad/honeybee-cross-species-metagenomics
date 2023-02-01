##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

##############
# files to be read
##############

df_gtdbk_bac <- read.csv("06_MAG_binning/gtdbtk_out_dir/classify/gtdbtk.bac120.summary.tsv", sep = "\t")
drep_Bdb <- read.csv("06_MAG_binning/drep_results/data_tables/Bdb.csv", sep = ",")
drep_gInformation <- read.csv("06_MAG_binning/drep_results/data_tables/genomeInfo.csv", sep = ",")
drep_Sdb <- read.csv("06_MAG_binning/drep_results/data_tables/Sdb.csv", sep = ",")
drep_Widb <- read.csv("06_MAG_binning/drep_results/data_tables/Widb.csv", sep = ",")
MAGs_collated <- read.csv("06_MAG_binning/all_GenomeInfo_auto.tsv", sep = "\t")

##############
# analyse data and plot
##############

cluster_info <- drep_Widb %>%
                  select(genome, score, cluster, cluster_members, furthest_cluster_member) %>%
                    rename(score_representative = score) %>%
                    mutate(representative_genome = Vectorize(format_genome_name)(genome)) %>%
                      select(!genome)
MAGs_collated_info <- left_join(MAGs_collated, mutate(drep_gInformation, genome = Vectorize(format_genome_name)(genome)), by = c("ID"="genome"))
# Number of MAGs coming from each host species
vis_magOTUs_df_all <- MAGs_collated_info %>%
                        group_by(Sample) %>%
                            mutate(all_quality = ifelse(completeness > 50 & contamination < 5 & N50 > 10000, "Pass", "Fail")) %>%
                            mutate(Completeness_quality = cut(completeness ,breaks=c(-1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100, 110),
                                                            labels = c("<10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-95", "95-100", "100"))
                                                          ) %>%
                              mutate(Contamination_quality = cut(contamination ,breaks=c(-1, 0, 5, 10, 100),
                                                              labels = c("0", "0-5", "5-10", ">10"))
                                                          ) %>%
                              mutate(N50_quality = cut(N50 ,breaks=c(0, 10000, 50000, 100000, 200000, 500000, 1000000, 2000000, Inf), labels = c("<10Kb", "10-50Kb", "50-100Kb", "100-200Kb", "200-500Kb", "0.5-1Mb", "1-2Mb", ">2Mb"))) %>%
                                mutate(sample = Vectorize(get_sample_name)(ID)) %>%
                                filter(grepl("MAG_", ID)) %>%
                                  mutate(Host = Vectorize(get_host_name)(ID)) # %>%
                                    # mutate(num_contigs = Vectorize(get_num_contigs_per_bin)(ID))
vis_magOTUs_df_all <- vis_magOTUs_df_all %>%
                          group_by(Cluster) %>%
                            mutate(Num_mags = n()) %>%
                              mutate(Prevalence_overall = Num_mags/length(samples))

vis_magOTUs_df_all <- vis_magOTUs_df_all %>%
                          group_by(Cluster, Host) %>%
                            mutate(Present = n()) %>%
                            mutate(Prevalence = ifelse(Host=="Apis mellifera", Present/length(samples_am), NA)) %>%
                            mutate(Prevalence = ifelse(Host=="Apis cerana", Present/length(samples_ac), Prevalence)) %>%
                            mutate(Prevalence = ifelse(Host=="Apis dorsata", Present/length(samples_ad), Prevalence)) %>%
                            mutate(Prevalence = ifelse(Host=="Apis florea", Present/length(samples_af), Prevalence)) %>%
                              # left_join(summarise(group_by(contigs_depths_df, bin), mean_coverage = mean(depth), .groups = "drop"), by = c("ID" = "bin")) %>%
                                arrange(Genus)
vis_magOTUs_df_all <- left_join(vis_magOTUs_df_all, cluster_info, by = c("Cluster" = "cluster")) %>%
                        left_join(drep_Sdb %>%
                                    mutate(ID = Vectorize(format_genome_name)(genome)) %>%
                                    select(!genome)
                        )
vis_magOTUs_df <- vis_magOTUs_df_all %>%
                    filter(all_quality == "Pass")
vis_magOTUs_df$Host <- as.factor(recode(vis_magOTUs_df$Host, Am="Apis mellifera", Ac="Apis cerana", Ad="Apis dorsata", Af="Apis florea"))
vis_magOTUs_df$Cluster <- as.factor(vis_magOTUs_df$Cluster)
vis_magOTUs_df$sample <- as.factor(vis_magOTUs_df$sample)
vis_magOTUs_df$Family <- as.factor(vis_magOTUs_df$Family)
vis_magOTUs_df$Genus <- as.factor(vis_magOTUs_df$Genus)
vis_magOTUs_df <- vis_magOTUs_df[order(vis_magOTUs_df$Genus), ]


write.csv(vis_magOTUs_df_all, row.names = F, file = "Figures/exported-vis_magOTUs_df_all.csv", quote = F)
write.csv(vis_magOTUs_df, row.names = F, file = "Figures/exported-vis_magOTUs_df.csv", quote = F)
