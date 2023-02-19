##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

condense_motu <- function(motu){
  motu_condensed = strsplit(motu, " ")[[1]][1]
  return(motu_condensed)
}

extend_colors_motus <- function(names_vec){
  final_list <- list()
  for (a_name in names_vec) {
    if (a_name %in% names(motuColors)) {
      final_list[a_name] = motuColors[a_name]
    } else {
      final_list[a_name] = "grey"
    }
  }
  return(final_list)
}

##############
# files to be read
##############

df_motus_raw <- read.csv("02_motus_profile/samples_merged.motus", sep = "\t", skip = 2, stringsAsFactors=FALSE)

##############
# analyse data and plot
##############

colnames(df_motus_raw)[1] <- "taxonomy"
df_motus_raw <- df_motus_raw %>%
              mutate(across(!c("taxonomy"), as.numeric)) %>%
              mutate(sum_ab = rowSums(across(where(is.numeric)))) %>%
                filter(sum_ab > 0) %>%
                  select(!sum_ab) %>%
                    column_to_rownames("taxonomy")

df_motus <- as.data.frame(t(df_motus_raw)) %>%
              rownames_to_column("Sample") %>%
              pivot_longer(!Sample, names_to = "motu", values_to = "rel_ab") %>%
                group_by(Sample, motu) %>%
                  mutate(Present = ifelse(rel_ab > 0, 1, 0)) %>%
                    group_by(motu) %>%
                     mutate(Prevalence_num = sum(Present)) %>%
                      mutate(Prevalence = mean(Present))

df_motus_combined <- left_join(df_motus, df_meta)  %>%
                      group_by(SpeciesID, motu) %>%
                        mutate(Present_host = ifelse(rel_ab > 0, 1, 0)) %>%
                          group_by(SpeciesID, motu) %>%
                            mutate(Prevalence_num_host = sum(Present_host)) %>%
                            mutate(Prevalence_host = mean(Present_host)) %>%
                              mutate(mean_rel_ab_host = mean(rel_ab)) %>%
                              group_by(Sample) %>%
                              mutate(mean_rel_ab = mean(rel_ab)) %>%
                                mutate(motu_condensed = Vectorize(condense_motu)(motu))

motu_list <- df_motus_combined %>%
              pull(motu) %>% unique

high_motu_list <- df_motus_combined %>%
                        filter(Prevalence_host >= 0.5) %>%
                          pull(motu_condensed) %>% unique

df_assigned_plot <- df_motus_combined %>%
                      group_by(Sample) %>%
                        mutate(sum = sum(rel_ab)) %>%
                          mutate(unassigned = rel_ab[motu == "unassigned"]) %>%
                            mutate("assigned" = sum - unassigned) %>%
                              select(Sample, assigned, unassigned) %>%
                                pivot_longer(!Sample, values_to = "proportion", names_to = "Type")

assigned_plot <- ggplot(df_assigned_plot, aes(y = factor(Sample, samples), x = proportion, fill = factor(Type, levels = c("unassigned", "assigned")))) +
                  geom_bar(position = "stack", stat = "identity") +
                  labs(fill = "Type", x = "Proportion", y = "Sample") +
                    make_theme(leg_pos = "right", guide_nrow = 20, leg_size = 7,
                               y_size = 8, y_hj = 1, y_vj = 0.5
                    )
                    # ggsave("Figures/02-motus_unassigned.pdf")

plot_motus_high_prev <- ggplot(df_motus_combined, aes(x = factor(Sample, samples), y = rel_ab, fill = factor(motu_condensed, motus))) +
                geom_bar(position = "stack", stat = "identity") +
                labs(fill = "mOTU", x = "mOTUs2 Relative abundance", y = "Sample") +
                  make_theme(max_colors = length(unique(motus)),
                             setFill = F,
                             leg_pos = "bottom",
                             guide_nrow = 3,
                             leg_size = 8,
                             x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                             y_angle=0 ,y_vj=0, y_hj=0, y_size=12
                  ) +
                  scale_fill_manual(values=motuColors)
                    ggsave("Figures/02-motus.pdf")

df_motus_filt <- filter(df_motus_combined, Prevalence_num >= 1 & rel_ab > 0.01)
plot_motus_filt <- ggplot(df_motus_filt, aes(y = factor(Sample, samples), x = rel_ab, fill = factor(motu_condensed))) +
                geom_bar(position = "stack", stat = "identity") +
                labs(fill = "mOTU", x = "mOTUs2 Relative abundance", y = "Sample") +
                  make_theme(max_colors = length(unique(df_motus_filt$motu_condensed)),
                             palettefill = "Set1",
                             leg_pos = "right",
                             guide_nrow = 20,
                             leg_size = 7,
                             y_size = 8, y_hj = 1,
                  )
                  ggsave("Figures/02-motus_filt.pdf")

ggplot(df_motus_combined %>% filter(rel_ab > 0.01), aes(y = factor(Sample, samples), x = rel_ab, fill = factor(motu))) +
                geom_bar(position = "stack", stat = "identity") +
                labs(fill = "mOTU", x = "mOTUs2 Relative abundance", y = "Sample") +
                  make_theme(max_colors = length(unique(df_motus_combined$motu)),
                             leg_pos = "none",
                             leg_size = 7,
                             y_size = 8, y_hj = 1
                  )
ggsave("Figures/02c-motus_all.pdf")
grid.arrange(get_only_legend(ggplot(df_motus_combined  %>% filter(rel_ab > 0.01), aes(y = factor(Sample, samples), x = rel_ab, fill = factor(motu))) +
                geom_bar(position = "stack", stat = "identity") +
                labs(fill = "mOTU", x = "mOTUs2 Relative abundance", y = "Sample") +
                  make_theme(max_colors = length(unique(df_motus_combined$motu)),
                             leg_size = 7,
                             guide_nrow = 40,
                             y_size = 8, y_hj = 1
                  ) +
                  theme(legend.key=element_blank(), legend.key.size=unit(5,"point"))
                )
              )
ggsave("Figures/02-motus_all_legend.pdf")

plot_motus_filt
plot_motus_high_prev