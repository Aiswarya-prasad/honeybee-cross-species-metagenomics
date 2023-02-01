#*********Work in progress*********#
# conda clean --lock
##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

##############
# files to be read
##############

source('scripts/visualization/05-read_MAG_metadata.R', chdir = TRUE)

genera_list <- c()

for (file_name in list.files("database/MAGs_database_Orthofinder/")) {
  if (!grepl(".csv", file_name)) {
    genus_name <- file_name
  }
  if (!(genus_name %in% genera_list)) {
    genera_list <- c(genera_list, genus_name)
  }
}
all_Isolates <- c()
all_MAGs <- c()
Number_genomes_df <- data.frame()
for (group in Groups) {
  OG_numbers_df <- read.csv(paste0("07_AnnotationAndPhylogenies/02_orthofinder/", group, "/OrthoFinder/Results_", group, "/Orthogroups/Orthogroups.GeneCount.tsv"), sep = "\t")
  genomes <- head(colnames(OG_numbers_df)[-1], -1)
  MAGs <- genomes[which(grepl("MAG", genomes, fixed=TRUE))]
  all_MAGs <- c(all_MAGs, MAGs)
  Isolates <- genomes[which(!grepl("MAG", genomes, fixed=TRUE))]
  all_Isolates <- c(all_Isolates, Isolates)
  Number_genomes_df <- rbind(Number_genomes_df,
                             cbind(
                                MAGs = length(MAGs),
                                Isolates = length(Isolates),
                                genomes = length(genomes),
                                Group = group
                              )
                            )
}
Number_genomes_df <- pivot_longer(Number_genomes_df, !Group, names_to = "Type", values_to = "Number") %>%
                      mutate(Number = as.numeric(Number))
ortho_mags_df <- data.frame()
ortho_mags_filt_df <- data.frame()
for (genus_name in genera_list) {
  orthofinder_filtered_summary <- read.csv(paste0("database/MAGs_database_Orthofinder/", genus_name,"_Orthogroups_filtered_summary.csv"))
  orthofinder_summary <- read.csv(paste0("database/MAGs_database_Orthofinder/", genus_name,"_Orthogroups_summary.csv"))
  ortho_mags_df <- rbind(ortho_mags_df, orthofinder_summary)
  ortho_mags_filt_df <- rbind(ortho_mags_filt_df, orthofinder_filtered_summary)
}
ortho_w_isolates_df <- data.frame()
for (group in Groups) {
  orthofinder_summary_w_isolates <- read.csv(paste0("07_AnnotationAndPhylogenies/02_orthofinder/", group, "_Orthogroups_summary.csv"))
  ortho_w_isolates_df <- rbind(ortho_w_isolates_df, orthofinder_summary_w_isolates)
}

Orthogroups_summary_df <- data.frame()
for (group in Groups) {
  OG_summary <- read.csv(paste0("07_AnnotationAndPhylogenies/02_orthofinder/", group, "_Orthogroups_summary.csv"), sep = ",")
  Orthogroups_summary_df <- rbind(Orthogroups_summary_df,
                                  OG_summary
                            ) %>%
                              filter(present_in != "None")
}

##############
# analyse data and plot
##############

# core genes used for trees
# ortho_w_isolates_df %>% pull(status) %>% unique()
ortho_w_isolates_df %>%
  filter(status %in% c("Isolates:ScpCore ; MAGs:ScpCore", "Isolates:ScpCore ; MAGs:ScpCore0.5")) %>%
    group_by(group) %>%
      summarise(genes_used = n())

plot_number_OGs <- ortho_w_isolates_df %>%
  mutate(Type = ifelse(status %in% c("Isolates:ScpCore ; MAGs:ScpCore0.5"), "Present in at least 50% MAGs", NA)) %>%
  mutate(Type = ifelse(status %in% c("Isolates:ScpCore ; MAGs:ScpCore"), "Present in all MAGs", Type)) %>%
  filter(!is.na(Type)) %>%
    group_by(group, Type) %>%
     summarise(Number = n())

# core genes used for core_cov
ortho_mags_df %>% pull(status) %>% unique()
ortho_mags_df %>%
  filter(status %in% c("; MAGs:ScpCore", "; MAGs:ScpCore0.5")) %>%
    group_by(group) %>%
      summarise(genes_used = n()) %>%
        print(n = 50)

plot_number_OGs_mags <- ortho_mags_df %>%
  mutate(Type = ifelse(status %in% c("; MAGs:ScpCore0.5"), "Present in at least 50% MAGs", NA)) %>%
  mutate(Type = ifelse(status %in% c("; MAGs:ScpCore"), "Present in all MAGs", Type)) %>%
  filter(!is.na(Type)) %>%
    group_by(group, Type) %>%
    filter(group %in% Groups) %>%
     summarise(Number = n())

ggplot() +
  geom_bar(data = plot_number_OGs_mags,
           aes(
            x = Number,
            y = group,
            fill = Type), stat = "identity") +
  labs(x = "Number of orthogroups", y = "Genus", fill = "Type") +
  xlim(c(0, 3500)) +
  scale_fill_manual(values = c(brewer.pal(9, "Pastel1")[1],
                               brewer.pal(9, "Pastel1")[2]
                              )
  ) +
  new_scale_fill() +
  geom_rect(data =plot_number_OGs_mags %>%
                    select(group) %>% unique %>%
                    left_join(Number_genomes_df %>%
                                pivot_wider(names_from = Type, values_from = Number), by = c("group" = "Group")) %>%
                        ungroup() %>%
                        filter(group %in% Groups) %>%
                           mutate(SampleID = row_number()), 
            aes(ymin=SampleID - 0.5,
                ymax=SampleID + 0.5,
                xmax = 3000,
                xmin = 3200,
                fill = MAGs)
          ) +
  labs(fill = "Number of MAGs") +
  scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Blues"), guide = "colourbar") +
  new_scale_fill() +
  geom_rect(data = plot_number_OGs_mags %>%
                    select(group) %>% unique %>%
                    left_join(Number_genomes_df %>%
                                pivot_wider(names_from = Type, values_from = Number), by = c("group" = "Group")) %>%
                        ungroup() %>%
                        filter(group %in% Groups) %>%
                           mutate(SampleID = row_number()), 
            aes(ymin=SampleID - 0.5,
                ymax=SampleID + 0.5,
                xmax = 3200,
                xmin = 3400,
                fill = Isolates)
          ) +
  labs(fill = "Number of Isolates") +
  scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar") +
  make_theme(leg_pos="right", setFill = F, modify_guide = F
  )

freq_MAGs_df_all <- data.frame()
freq_Isolates_df_all <- data.frame()
for (group in Groups) {
  OG_numbers_df <- read.csv(paste0("07_AnnotationAndPhylogenies/02_orthofinder/", group, "/OrthoFinder/Results_", group, "/Orthogroups/Orthogroups.GeneCount.tsv"), sep = "\t")
  OG_numbers_df_MAGs <- OG_numbers_df %>% select(Orthogroup | Total | starts_with("MAG"))
  OG_numbers_df_Isolates <- OG_numbers_df %>% select(Orthogroup | !starts_with("MAG"))
  freq_MAGs_df <- OG_numbers_df_MAGs %>%
                # single-copy and half core
                filter(Total <= sum(colnames(.) %in% all_MAGs) - 2) %>%
                filter(Total > (sum(colnames(.) %in% all_MAGs) - 2)/2 ) %>%
                                  select(Total) %>%
                                    mutate(Type = "MAGs") %>%
                                      mutate(Total_prop = Total/max(Total))
  freq_Isolates_df <- OG_numbers_df_Isolates %>%
                    # single-copy and half core
                    filter(Total <= sum(colnames(.) %in% all_Isolates) - 2) %>%
                    filter(Total == (sum(colnames(.) %in% all_Isolates) - 2)) %>%
                      select(Total) %>%
                        mutate(Type = "Isolates") %>%
                          mutate(Total_prop = Total/max(Total))
  freq_MAGs_df_all <- rbind(freq_MAGs_df_all, cbind(freq_MAGs_df, "Group" = rep(group, dim(freq_MAGs_df)[[1]])))
  freq_Isolates_df_all <- rbind(freq_Isolates_df_all, cbind(freq_Isolates_df, "Group" = rep(group, dim(freq_Isolates_df)[[1]])))
}

ggplot(freq_MAGs_df_all, aes(x = Total_prop, fill = Type)) +
    geom_histogram(binwidth=0.05) +
    labs(x = "Proportion of genomes containing orthogroup", y = "Frequency") +
    facet_wrap(~ Group) +
      make_theme(palettefill="Pastel1")

orthogroups_present_in <- Orthogroups_summary_df %>%
  group_by(group, present_in) %>%
    summarise(Number = n(), .groups="drop")

num_genomes <- ggplot(Number_genomes_df %>% filter(Type != "genomes"),
                      aes(x = Number, y = factor(Group, genera), fill = Type)) +
                geom_bar(stat = "identity") +
                labs(x = "Number of genomes", y = "Genus group") +
                make_theme(leg_pos="right",
                           x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                           y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                ggsave("Figures/07a-Number_of_genomes_per_group.pdf")


plot_present_in <- ggplot(orthogroups_present_in,
                          aes(x = Number, y = factor(group, genera), fill = present_in)) +
                    geom_bar(stat = "identity") +
                    labs(x = "Number of orthogroups", y = "Genus group", fill = "Present in") +
                    make_theme(leg_pos="right",
                               x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                               y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                               ggsave("Figures/07b-Number_of_orthogroups_per_group.pdf")

ScpCore <- c("Isolates:ScpCore ; MAGs:ScpCore"
           )
ScpCoreMAGs <- c("Isolates:-Core ; MAGs:ScpCore",
                 "Isolates:Scp- ; MAGs:ScpCore",
                 "Isolates:-- ; MAGs:ScpCore"
           )
ScpCoreIsolates <- c("Isolates:ScpCore ; MAGs:-Core",
                     "Isolates:ScpCore ; MAGs:Scp-",
                     "Isolates:ScpCore ; MAGs:--",
                     "Isolates:ScpCore ; MAGs:-Core0.5"
           )
ScpCorehalfMAGs <- c("Isolates:ScpCore ; MAGs:ScpCore0.5"
           )

orthogroups_status <- Orthogroups_summary_df %>%
                            group_by(group, status) %>%
                            mutate(status = ifelse(status %in% c(ScpCore, ScpCorehalfMAGs, ScpCoreMAGs, ScpCoreIsolates), status, "Others")) %>%
                            mutate(status = ifelse(status %in% ScpCore, "Single-copy core orthogroups", status)) %>%
                            mutate(status = ifelse(status %in% ScpCorehalfMAGs, "Single-copy core in half of the MAGs and single-copy core in isolates", status)) %>%
                            mutate(status = ifelse(status %in% ScpCoreMAGs, "Single-copy core only in MAGs", status)) %>%
                            mutate(status = ifelse(status %in% ScpCoreIsolates, "Single-copy core only in Isolates", status)) %>%
                              summarise(Number = n(), .groups="drop")
status <- ggplot(orthogroups_status,
                  aes(x = Number, y = factor(group, genera), fill = status)) +
                geom_bar(stat = "identity", position = "stack") +
                labs(x = "Number of orthogroups", y = "Genus group", fill = "Core in") +
                scale_x_continuous(limits = c(0, 1500)) +
                make_theme(leg_pos="bottom", leg_size = 8, guide_nrow = 5,
                           max_colors = length(unique(orthogroups_status$status)),
                           palettefill = "Spectral",
                           x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                           y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
score_hist <- ggplot(vis_magOTUs_df_all, aes(x = score, fill = factor(Host, host_order))) +
              geom_histogram(binwidth = 1) +
              labs(x = "Score", y = "Number of MAGs", fill = "Host") +
              geom_vline(xintercept = 15) +
                make_theme(setF = F,
                           leg_pos="right", leg_size = 7,
                           guide_nrow=4
                ) +
                  scale_fill_manual(values=host_order_color)
                  ggsave("Figures/11a-mag_drep_scores_histogram.pdf")
cluster_size_histogram <- ggplot(vis_magOTUs_df_all, aes(x = cluster_members, fill = factor(Genus, genera))) +
  geom_histogram(binwidth = 0.5) +
  labs(x = "Number of members", y = "Number of MAGs", fill = "Genus") +
    make_theme(setF = F,
               leg_pos="right", leg_size = 7,
               guide_nrow= 21
    ) +
      scale_fill_manual(values=genusColors)
      ggsave("Figures/11b-cluster_size_histogram.pdf")

get_num_clusters_and_MAGs_CMT <- function(df_info, cut_off_value){
  df_info %>%
    ungroup() %>%
      mutate(CMT = completeness - 5*contamination) %>%
        filter(CMT > cut_off_value) %>%
          summarise(MAGs = n_distinct(ID), Clusters = n_distinct(Cluster)) %>%
            as.list()
}
get_num_clusters_and_MAGs_CMT(vis_magOTUs_df_all, 0)[["MAGs"]]
get_num_clusters_and_MAGs_CMT(vis_magOTUs_df_all, 0)[["Clusters"]]

ggplot(vis_magOTUs_df_all, aes(y = Cluster, fill = factor(Genus, genera))) +
  geom_bar() +
  labs(x = "Number of MAGs", y = "Cluster", fill = "Genus") +
    make_theme(setF = F,
               leg_pos="right", leg_size = 7,
               guide_nrow= 21, y_size = 3
    ) +
    geom_vline(xintercept = 1) +
    facet_grid(~Host) +
      scale_fill_manual(values=genusColors)
cluster_completeness_number <- ggplot(vis_magOTUs_df_all, aes(y = Completeness_quality, x=cluster_members, fill = factor(Genus, genera))) +
  geom_bar(stat = "identity") +
  labs(x = "Completeness", y = "Cluster", fill = "Genus") +
    make_theme(setF = F,
               leg_pos="right", leg_size = 7,
               guide_nrow= 21, y_size = 3
    ) +
    geom_vline(xintercept = 1) +
    facet_wrap_paginate(~Cluster, ncol = 5, nrow = 7,
                      scales = "fixed", strip.position = "top", page = 5) +
      scale_fill_manual(values=genusColors)

paginate_save(cluster_completeness_number, "Cluster", "Figures/11x-clusters_completeness_numbers.pdf",
              pass_nrow = 7)