#*********Work in progress*********#

##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

##############
# files to be read
##############

source('scripts/visualization/05-read_MAG_metadata.R', chdir = TRUE)
# system("bash scripts/calc_genome_sizes.sh")
genome_lengths_df <- read.csv("/scratch/aprasad/211018_Medgenome_india_samples/Figures/Genome_sizes.csv")
reference_genomes_info <- read.csv("/scratch/aprasad/211018_Medgenome_india_samples/06_MAG_binning/all_GenomeInfo_auto.csv")
drep_N <- read.csv("06_MAG_binning/drep_results/data_tables/Ndb.csv") %>%
            mutate(reference = Vectorize(format_genome_name)(reference)) %>%
              mutate(querry = Vectorize(format_genome_name)(querry))
##############
# analyse data and plot
##############

vis_magOTUs_df_all_means <- vis_magOTUs_df_all %>%
                        group_by(Sample) %>%
                          summarise(Sample, completeness, contamination, N50, Host) %>%
                            mutate(Completeness_mean = mean(completeness)) %>%
                              mutate(Contamination_mean = mean(contamination)) %>%
                              mutate(N50_mean = mean(N50))
Completeness_hist <- ggplot(vis_magOTUs_df_all, aes(x = completeness, fill = Host)) +
    geom_histogram(binwidth=2) +
    geom_vline(xintercept = 50) +
      make_theme(palettefill="Spectral")
N50_hist <- ggplot(vis_magOTUs_df_all, aes(x = N50, fill = Host)) +
    geom_histogram(bins = 150) +
      scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
      geom_vline(xintercept = 10000) +
      make_theme(palettefill="Spectral")
Contamination_hist <- ggplot(vis_magOTUs_df_all, aes(x = contamination, fill = Host)) +
    geom_histogram(binwidth=2) +
    geom_vline(xintercept = 5) +
      make_theme(palettefill="Spectral")
completeness_per_sample <- ggplot(vis_magOTUs_df_all, aes(y = factor(Sample, levels = samples), fill = Completeness_quality)) +
    geom_bar(position = "stack") +
    labs(fill = "Quality", y = "Sample") +
      make_theme(palettefill="RdYlGn", max_colors = length(levels(vis_magOTUs_df_all$Completeness_quality)))
N50_per_sample <- ggplot(vis_magOTUs_df_all, aes(y = factor(Sample, levels = samples), fill = N50_quality)) +
        geom_bar(position = "stack") +
        labs(fill = "Quality", y = "Sample") +
        make_theme(palettefill="RdYlBu", max_colors = length(levels(vis_magOTUs_df_all$N50_quality)))
contamination_per_sample <- ggplot(vis_magOTUs_df_all, aes(y = factor(Sample, levels = samples), fill = Contamination_quality)) +
    geom_bar(position = "stack") +
    labs(fill = "Quality", y = "Sample") +
      make_theme(palettefill="RdYlBu", max_colors = length(levels(vis_magOTUs_df_all$Contamination_quality)))
MAG_quality_per_sample <- ggplot(vis_magOTUs_df_all, aes(y = factor(Sample, levels = samples), fill = all_quality)) +
    geom_bar(position = "stack") +
    labs(fill = "Quality", y = "Sample") +
      make_theme(palettefill="Set1",)

prev_vs_abud_all <- ggplot(vis_magOTUs_df_all, aes(x = mean_coverage, y = Prevalence, size = completeness, color = Genus, alpha = 0.5)) +
                    geom_point(position = position_jitter(w = 0, h = 0.05)) +
                      make_theme(setFill = F, setCol = F,
                        leg_pos = "bottom",
                        guide_nrow = 8,
                        leg_size = 12
                      ) +
                      theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
                      scale_color_manual(values=genusColors) +
                      scale_x_continuous(trans = "log10") +
                        scale_alpha(guide = "none")

prev_vs_abud <- ggplot(vis_magOTUs_df, aes(x = mean_coverage, y = Prevalence, size = completeness, color = Genus, alpha = 0.5)) +
  geom_point(position = position_jitter(w = 0, h = 0.05)) +
    make_theme(setFill = F, setCol = F,
      leg_pos = "bottom",
      guide_nrow = 8,
      leg_size = 12
    ) +
    theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
    scale_color_manual(values=genusColors) +
    scale_x_continuous(trans = "log10") +
    # facet_wrap(~ factor(Host, host_order)) +
      scale_alpha(guide = "none")

prev_overall_vs_abud_all <- ggplot(vis_magOTUs_df_all, aes(x = mean_coverage, y = Prevalence_overall, size = completeness, color = Genus, alpha = 0.5)) +
                    geom_point(position = position_jitter(w = 0, h = 0.05)) +
                      make_theme(setFill = F, setCol = F,
                        leg_pos = "bottom",
                        guide_nrow = 8,
                        leg_size = 12
                      ) +
                      theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
                      scale_color_manual(values=genusColors) +
                      scale_x_continuous(trans = "log10") +
                        scale_alpha(guide = "none")

prev_overall_vs_abud <- ggplot(vis_magOTUs_df, aes(x = mean_coverage, y = Prevalence_overall, size = completeness, color = Genus, alpha = 0.5)) +
  geom_point(position = position_jitter(w = 0, h = 0.05)) +
    make_theme(setFill = F, setCol = F,
      leg_pos = "bottom",
      guide_nrow = 8,
      leg_size = 12
    ) +
    theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
    scale_color_manual(values=genusColors) +
    scale_x_continuous(trans = "log10") +
    # facet_wrap(~ factor(Host, host_order)) +
      scale_alpha(guide = "none")

prev_vs_abud_all_host <- ggplot(vis_magOTUs_df_all, aes(x = mean_coverage, y = Prevalence, size = completeness, color = Genus, alpha = 0.5)) +
  geom_point(position = position_jitter(w = 0, h = 0.05)) +
    make_theme(setFill = F, setCol = F,
      leg_pos = "none",
      guide_nrow = 8,
      leg_size = 12
    ) +
    ggtitle("MAGs with > 70% completeness and < 5% contamination") +
    theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
    scale_color_manual(values=genusColors) +
    scale_x_continuous(trans = "log10") +
    facet_wrap(~ factor(Host, host_order)) +
      scale_alpha(guide = "none")
      ggsave("Figures/06-prev_vs_coverage_all_MAGs_genus_by_host.pdf")

prev_vs_abud_host <- ggplot(vis_magOTUs_df, aes(x = mean_coverage, y = Prevalence, size = completeness, color = Genus, alpha = 0.5)) +
  geom_point(position = position_jitter(w = 0, h = 0.05)) +
    make_theme(setFill = F, setCol = F,
      leg_pos = "none",
      guide_nrow = 8,
      leg_size = 12
    ) +
    ggtitle("MAGs with > 70% completeness and < 5% contamination") +
    theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
    scale_color_manual(values=genusColors) +
    scale_x_continuous(trans = "log10") +
    facet_wrap(~ factor(Host, host_order)) +
      scale_alpha(guide = "none")
      ggsave("Figures/06-prev_vs_coverage_filtered_MAGs_genus_by_host.pdf")

prev_vs_abud_all_host <- ggplot(vis_magOTUs_df_all, aes(x = mean_coverage, y = Prevalence, color = Genus)) +
  geom_point(position = position_jitter(w = 0, h = 0.05)) +
    make_theme(setFill = F, setCol = F,
      leg_pos = "none",
      guide_nrow = 8,
      leg_size = 12
    ) +
    ggtitle("MAGs with > 70% completeness and < 5% contamination") +
    theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
    scale_color_manual(values=genusColors) +
    scale_x_continuous(trans = "log10") +
    facet_wrap(~ factor(Host, host_order)) +
      scale_alpha(guide = "none")
      ggsave("Figures/06-prev_vs_coverage_all_MAGs_genus_by_host_no_size.pdf")

prev_vs_abud_host <- ggplot(vis_magOTUs_df, aes(x = mean_coverage, y = Prevalence, color = Genus)) +
  geom_point(position = position_jitter(w = 0, h = 0.05)) +
    make_theme(setFill = F, setCol = F,
      leg_pos = "none",
      guide_nrow = 8,
      leg_size = 12
    ) +
    ggtitle("MAGs with > 70% completeness and < 5% contamination") +
    theme(legend.margin=margin(-1,-1,-1,-1), legend.box="vertical") +
    scale_color_manual(values=genusColors) +
    scale_x_continuous(trans = "log10") +
    facet_wrap(~ factor(Host, host_order)) +
      scale_alpha(guide = "none")
      ggsave("Figures/06-prev_vs_coverage_filtered_MAGs_genus_by_host_no_size.pdf")

legend_hist_host <- get_only_legend(Completeness_hist)
g  <- arrangeGrob(
  arrangeGrob(
      Completeness_hist + make_theme(setFill = F, setCol = F, leg_pos = "none"),
      Contamination_hist + make_theme(setFill = F, setCol = F, leg_pos = "none"),
      N50_hist + make_theme(setFill = F, setCol = F, leg_pos = "none"),
      nrow = 2,
      layout_matrix = rbind(c(1,2), c(3,3))
    ),
    legend_hist_host,
    heights = c(10, 1)
  )
  ggsave("Figures/06-QC_MAG_histogram.pdf", g)
g <- grid.arrange(
    MAG_quality_per_sample + make_theme(setFill = F, setCol = F, leg_size = 7, y_size = 5, leg_pos = "right", guide_nrow = 2),
    N50_per_sample + make_theme(setFill = F, setCol = F, leg_size = 7, y_size = 5, leg_pos = "right", guide_nrow = 7),
    contamination_per_sample + make_theme(setFill = F, setCol = F, leg_size = 7, y_size = 5, leg_pos = "right", guide_nrow = 4),
    completeness_per_sample + make_theme(setFill = F, setCol = F, leg_size = 7, y_size = 5, leg_pos = "right", guide_nrow = 10)
  )
  ggsave("Figures/06-QC_MAG_per_sample.pdf", g)
g <- grid.arrange(prev_vs_abud_all + make_theme(setFill = F, setCol = F, leg_pos = "none"),
             prev_vs_abud + make_theme(setFill = F, setCol = F, leg_pos = "none") + ggtitle("MAGs with > 70% completeness and < 5% contamination"),
             genus_legend,
             heights = c(3,3,2)
           )
    ggsave("Figures/06-prev_vs_coverage_MAGs_genus.pdf", g)
g <- grid.arrange(prev_overall_vs_abud_all + make_theme(setFill = F, setCol = F, leg_pos = "none"),
             prev_overall_vs_abud + make_theme(setFill = F, setCol = F, leg_pos = "none") + ggtitle("MAGs with > 70% completeness and < 5% contamination"),
             genus_legend,
             heights = c(3,3,2)
           )
    ggsave("Figures/06-prev_overall_vs_coverage_MAGs_genus.pdf", g)
vis_magOTUs_df %>%
  group_by(Host) %>%
  # filter(mean_coverage > 0) %>%
    mutate(N_MAGs = n_distinct(ID)) %>%
    mutate(N_magOTUs = n_distinct(Cluster)) %>%
    mutate(N_Genera = n_distinct(Genus)) %>%
      summarise(N_MAGs, N_magOTUs, N_Genera) %>%
        unique()

genus_MAG_quality_host_all <- ggplot(vis_magOTUs_df_all, aes(y = Genus, fill = Host)) +
        geom_bar(position = "stack") +
        labs(fill = "Host", y = "Genus", x = "Number of MAGs") +
        make_theme(palettefill = "Spectral")
      ggsave("Figures/06-QC_per_Genus_per_host_all.pdf")

genus_MAG_quality_host <- ggplot(filter(vis_magOTUs_df_all, all_quality=="Pass"), aes(y = Genus, fill = Host)) +
        geom_bar(position = "stack") +
        labs(fill = "Host", y = "Genus", x = "Number of MAGs") +
        make_theme(palettefill = "Spectral")
      ggsave("Figures/06-QC_per_Genus_per_host_passed.pdf")

line_list <- c()
for (num in 1:94){
  add_line <- geom_vline(xintercept=num+0.5, size=0.1, color="black")
  # add_line <- geom_vline(xintercept=num+0.5, size=0.1, alpha=0.5, color="grey")
  line_list <- c(line_list, add_line)
}

magOTUs_per_sample <- ggplot(vis_magOTUs_df_all, aes(y = factor(Cluster), x = factor(sample, samples), fill = Host)) +
                            geom_tile() +
                              labs(x = "Sample", y = "Cluster")+
                              make_theme(setFill=F,
                              # make_theme(palettefill="Spectral", max_colors = length(unique(vis_magOTUs_df$Cluster)),
                              leg_pos="none", guide_nrow=6,
                              y_hj=1, y_size=7, leg_size=8, y_vj=0.5,
                              x_vj=0, x_hj=1, x_size=6, x_angle=90) +
                              scale_fill_manual(values=host_order_color) +
                              line_list
                                    ggsave("Figures/06-magOTUs_per_sample.pdf")

magOTUs_per_sample_genus <- ggplot(vis_magOTUs_df_all, aes(y = factor(Cluster), x = factor(sample, samples), fill = factor(Genus, genera))) +
                            geom_tile() +
                              labs(x = "Sample", y = "Cluster")+
                              make_theme(setFill=F,
                              # make_theme(palettefill="Spectral", max_colors = length(unique(vis_magOTUs_df$Cluster)),
                              leg_pos="none", guide_nrow=6,
                              y_hj=1, y_size=7, leg_size=8, y_vj=0.5,
                              x_vj=0, x_hj=1, x_size=6, x_angle=90) +
                              scale_fill_manual(values=genusColors, guide = F) +
                              line_list
                                    ggsave("Figures/06-magOTU_per_sample_genus.pdf")

magOTUs_per_sample_by_host_genus <- ggplot(vis_magOTUs_df_all, aes(y = Cluster, x = sample, fill = factor(Genus, genera))) +
        geom_tile() +
        labs(y = "Cluster", x = "Prevalence", fill = "Genus") +
        make_theme(setFill = F, setCol = F,
                   y_size = 2, y_hj = 1.5, y_vj = 0.5,
                   x_size = 7, x_angle = 40, x_hj = 1, x_vj = 1,
                   leg_size = 5, leg_pos = "none") +
        scale_fill_manual(values=genusColors) +
          facet_wrap(~ factor(Host, host_order), scales = "free")
      ggsave("Figures/06-magOTU_by_host_genus.pdf")
# contigs_depths_df %>%
#   filter(bin == "MAG_C2.4_8")
#
# vis_magOTUs_df_all %>%
#   filter(Cluster == "116_1") %>%
#     select(Class, Genus, ID, completeness, length)

magOTUs_per_sample_by_host_coverage <- ggplot(vis_magOTUs_df_all, aes(y = Cluster, x = sample, fill = mean_coverage)) +
        geom_tile() +
        labs(y = "Cluster", x = "Prevalence", fill = "Log of mean of contig mean coverage") +
        make_theme(setFill = F, setCol = F,
                   y_size = 7, y_hj = 1, y_vj = 0.5,
                   x_size = 7, x_angle = 40, x_hj = 1, x_vj = 1,
                   guide_nrow = 1,
                   leg_pos = "bottom"
                 ) +
          guides(fill = guide_colorbar(barhwight = 1, barwidth = 10)) +
          scale_fill_gradientn(colors=brewer.pal(5, "RdYlGn"), na.value = "transparent",
                              trans = "log10") +
          facet_wrap(~ factor(Host, host_order), scales = "free")
          ggsave("Figures/06-magOTU_by_host_MAGs_coverage.pdf")

magOTUs_per_sample_by_host_prevalence <- ggplot(vis_magOTUs_df_all, aes(y = Cluster, x = sample, fill = Prevalence)) +
        geom_tile() +
        labs(y = "Cluster", x = "Prevalence", fill = "Prevalence within host") +
        make_theme(setFill = F, setCol = F,
                   y_size = 7, y_hj = 1, y_vj = 0.5,
                   x_size = 7, x_angle = 40, x_hj = 1, x_vj = 1,
                   guide_nrow = 1,
                   leg_pos = "bottom"
                 ) +
          guides(fill = guide_colorbar(barhwight = 1, barwidth = 10)) +
          scale_fill_gradientn(colors=brewer.pal(5, "RdYlGn"), na.value = "transparent") +
          facet_wrap(~ factor(Host, host_order), scales = "free")
          ggsave("Figures/06-magOTU_by_host_MAGs_prevalence.pdf")

vis_magOTUs_df_all_shared_cluster <- vis_magOTUs_df_all %>%
                              group_by(Cluster) %>%
                                mutate(Num_hosts = n_distinct(Host)) %>%
                                  filter(Num_hosts > 1)

magOTUs_shared_per_sample_genus <- ggplot(vis_magOTUs_df_all_shared_cluster, aes(y = factor(Cluster), x = factor(sample, samples), fill = factor(Genus, genera))) +
                            geom_tile() +
                              labs(x = "Sample", y = "Cluster")+
                              make_theme(setFill=F,
                              leg_pos="none", guide_nrow=6,
                              y_hj=1, y_size=7, leg_size=8, y_vj=0.5,
                              x_vj=0, x_hj=1, x_size=6, x_angle=90) +
                              scale_fill_manual(values=genusColors) +
                              line_list
                                    ggsave("Figures/06-magOTU_shared_per_sample_genus.pdf")

magOTUs_shared_per_sample_prev_abund <- ggplot(vis_magOTUs_df_all_shared_cluster, aes(y = factor(Cluster), x = factor(sample, samples), fill = Prevalence)) +
                            geom_tile() +
                            geom_text(aes(label = round(mean_coverage, 2)), size = 1) +
                              labs(x = "Sample", y = "Cluster", fill = "Prevalence within host")+
                              make_theme(setFill=F,
                              y_hj=1, y_size=7, leg_size=8, y_vj=0.5,
                              x_vj=0, x_hj=1, x_size=6, x_angle=90) +
                              guides(fill = guide_colorbar(barhwight = 1, barwidth = 10)) +
                              scale_fill_gradientn(colors=brewer.pal(5, "RdYlGn"), na.value = "transparent", trans = "log10") +
                              line_list
                                    ggsave("Figures/06-magOTUs_shared_per_sample_prev_abund.pdf")

magOTUs_per_sample_by_host_completeness_genus <- ggplot(vis_magOTUs_df_all, aes(y = Cluster, x = Num_mags, size = factor(Completeness_quality), color = Genus, alpha = 0.5)) +
        geom_point() +
        labs(y = "Cluster", x = "Number of MAGs", size = "Completeness") +
        make_theme(setFill = F, setCol = F,
                   y_size = 3, y_hj = 1, y_vj = 0.5,
                   leg_size = 5, leg_pos = "right") +
        scale_color_manual(values=genusColors) +
          facet_wrap(~ factor(Host, host_order)) +
            guides(color = "none", alpha = "none")
      ggsave("Figures/06-magOTUs_by_host_completeness.pdf")

prev_abd_by_genus_completeness <- ggplot(vis_magOTUs_df_all, aes(y = Prevalence, x = mean_coverage, color = factor(Completeness_quality), alpha = 0.5)) +
        geom_point(position = position_jitter(w = 0, h = 0.05)) +
        labs(y = "Prevalence within host", x = "Mean of mean contig coverage", color = "Completeness") +
        make_theme(setFill = F, setCol = T,
                   palettecolor = "RdYlGn",
                   # y_size = 3, y_hj = 1, y_vj = 0.5,
                   x_angle = 30, x_hj = 1, x_vj = 1,
                   leg_size = 8, leg_pos = "right",
                   guide_nrow = 11
                 ) +
        scale_x_continuous(trans="log10") +
          facet_wrap(~ Genus) +
            guides(alpha = "none", color = "none")
            ggsave("Figures/06-mags_prev_vs_abd_by_genus_completeness.pdf")

prev_abd_by_genus_host <- ggplot(vis_magOTUs_df_all, aes(y = Prevalence, x = mean_coverage, color = Host, alpha = 0.5)) +
        geom_point(position = position_jitter(w = 0, h = 0.05)) +
        labs(y = "Prevalence within host", x = "Mean of mean contig coverage", color = "Host species") +
        make_theme(setFill = F, setCol = F,
                   # y_size = 3, y_hj = 1, y_vj = 0.5,
                   x_angle = 30, x_hj = 1, x_vj = 1,
                   leg_size = 5, leg_pos = "right") +
        scale_color_manual(values=host_order_color) +
        scale_x_continuous(trans="log10") +
          facet_wrap(~ Genus) +
            guides(color = "none", alpha = "none")
            ggsave("Figures/06-mags_prev_vs_abd_by_genus_host.pdf")

magOTUs_per_sample
magOTUs_per_sample_genus
magOTUs_per_sample_by_host_genus
magOTUs_per_sample_by_host_coverage
magOTUs_per_sample_by_host_prevalence
magOTUs_shared_per_sample_genus
magOTUs_shared_per_sample_prev_abund
magOTUs_per_sample_by_host_completeness_genus
prev_abd_by_genus_completeness
prev_abd_by_genus_host

reference_genomes <- reference_genomes_info %>%
                            filter(Source_database != "MAGs" & !grepl("MAG_", ID)) %>%
                            select(ID,Group_auto,Host,Species,Phylotype) %>%
                              rename(Genus = Group_auto) %>%
                                left_join(genome_lengths_df) %>%
                                  group_by(Genus) %>%
                                  mutate(length_med = median(length))

ggplot() +
  geom_boxplot(data = vis_magOTUs_df %>% filter(!is.na(Genus)),
               aes(y = Genus,
                   x = length,
                  #  fill = Genus
               ),
               outlier.shape = NA
               ) +
    scale_fill_manual(values=genusColors) +
    geom_jitter(data = vis_magOTUs_df %>% filter(!is.na(Genus)),
                aes(y = Genus,
                    x = length,
                    color = factor(Ref_status)
                )
               ) +
    scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
    make_theme(palettefill = "Spectral", setCol = F, setFill = F,
               guide_nrow = 1,
               x_size = 14,
               y_size = 14
              ) +
    labs(x = "MAG size / Genome size", y = "Genus", color = "MAG type") +
    geom_point(data = reference_genomes,
                aes(y = Genus,
                    x = length
                ), size = 2,
                # width = 0.01,
                shape = 3, color = "red"
               ) +
      scale_color_manual(values = c("1" = "blue", "0" = "black"),
                         label = c("1" = "Best scoring MAG", "0" = "Other MAGs")
      ) +
    labs(x = "MAG size / Genome size", y = "Genus", color = "MAG type") +
    guides(fill = FALSE) +
    geom_text(aes(x = 3.5e6, y = 1, label = "+ "), color = "red") +
    geom_text(aes(x = 4.5e6, y = 1, label = "represents reference isolate genomes"))
    ggsave("Figures/06-Genome_sizes_by_genus.pdf")


# maybe do this genus-wise
matrix <- drep_N %>%
                  select(reference, querry, ani) %>%
                  filter(reference %in% vis_magOTUs_df$ID) %>%
                    right_join(vis_magOTUs_df %>% ungroup() %>% select(ID, Cluster, Genus), by = c("reference" = "ID")) %>%
                      # filter(Genus %in% c("g__Bifidobacterium")) %>%
                      mutate(reference = paste0(Cluster, "_", Genus, "_", reference)) %>%
                      mutate(querry = paste0(Cluster, "_", Genus, "_", querry)) %>%
                      select(!c(Genus, Cluster)) %>%
                        pivot_wider(names_from=querry, values_from=ani) %>%
                          replace(is.na(.), 0.7) %>%
                            column_to_rownames("reference") %>%
                            select(rownames(.)) %>%
                                as.matrix
drep_dist_ordered <- matrix %>%
                      as.data.frame() %>%
                        rownames_to_column(var = "reference") %>%
                          pivot_longer(!reference, values_to = "ani", names_to = "querry") %>%
                            mutate(ani = as.numeric(ani))


drep_dist <- drep_N %>%
                  select(reference, querry, ani) %>%
                  filter(reference %in% vis_magOTUs_df$ID) %>%
                  left_join(vis_magOTUs_df %>% ungroup() %>% select(ID, Cluster, Genus), by = c("reference" = "ID"))



ggplot(drep_dist %>% filter(!is.na(Genus)), aes(x = reference, y = querry)) +
  geom_tile(aes(fill = ani)) +
    make_theme(setFill = F, modify_guide = F, leg_pos = "right",
               x_size = 0, y_size = 0
  ) +
  # facet_wrap(~ factor(Genus, genera), scales="free") +
    scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#66bd63", midpoint = 0.85, guide = "colourbar", limits=c(0.7,1), na.value = "#d73027")


column_ha = HeatmapAnnotation(Host = rownames(abundance_df_matrix_pa[samples_IN,]) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color_dark)
                          )
column_ba = HeatmapAnnotation(Location = rownames(abundance_df_matrix_pa[samples_IN,]) %>% as.data.frame() %>% mutate(Location = Vectorize(get_location_from_sample_name)(.)) %>% pull(Location),
                              col = list(Location = location_country_colors)
                          )
row_ha = rowAnnotation(Genus = colnames(abundance_df_matrix_pa[samples_IN,]) %>% as.data.frame() %>% left_join(cluster_tax_info %>% ungroup() %>% select(Genus, cluster_name), by = c(. = "cluster_name")) %>% pull(Genus),
                       col = list(Genus = unlist(genusColors))
                    )
# pdf("Figures/09-presence_absence_heatmap_IN.pdf")

make_col_list <- function(my_vec, my_palette = "Spectral") {
  uniq_names <- unique(my_vec)
  num_ele <- length(uniq_names)
  if (num_ele <= 3) {
    cols_used <- brewer.pal(5, my_palette)[1:num_ele]
  }
  if (num_ele <= 9 & num_ele > 3) {
    cols_used <- brewer.pal(num_ele, my_palette)
  } else {
    if (my_palette == "Pastel2") {
      cols_used <- colorRampPalette(brewer.pal(8, my_palette))(num_ele)
    } else {
      cols_used <- colorRampPalette(brewer.pal(9, my_palette))(num_ele)
    }
  }
  col_list <- c(cols_used)
  names(col_list) = uniq_names
  return(col_list)
}

plot_ani_heatmap <- function(matrix, row_split = "magotus", add_values = F, value_size = 10, limit_l = 0.95, limit_h = 1, col_fun = colorRamp2(c(limit_l, limit_h), c("#9ecae1", "#08306b"))) {
  host_names <- rownames(matrix) %>% as.data.frame() %>%
                  as.data.frame() %>%
                        left_join(vis_magOTUs_df %>%
                          ungroup() %>%
                            select(ID, Host), by = c(. = "ID")) %>% pull(Host)
  magOTU_names <- rownames(matrix) %>% 
                      as.data.frame() %>%
                        left_join(vis_magOTUs_df %>%
                          ungroup() %>%
                            select(ID, Cluster), by = c(. = "ID")) %>% pull(Cluster)
  Genus_names <- rownames(matrix) %>% 
                    as.data.frame() %>%
                      left_join(vis_magOTUs_df %>%
                        ungroup() %>%
                          select(ID, Genus), by = c(. = "ID")) %>% pull(Genus)
  species_names <- rownames(matrix) %>% 
                    as.data.frame() %>%
                      left_join(vis_magOTUs_df %>%
                        ungroup() %>%
                          select(ID, Species), by = c(. = "ID")) %>% pull(Species)
  anno_host_col = HeatmapAnnotation(Host = host_names,
                                      col = list(Host = host_order_color_dark),
                                      border = TRUE
                               )
  anno_host_row = rowAnnotation(Host = host_names,
                                      col = list(Host = host_order_color_dark),
                                      border = TRUE
                                                  )
  anno_gtdb_col = HeatmapAnnotation(`GTDB Species` = species_names,
                                    col = list(`GTDB Species` = make_col_list(species_names, "Paired")),
                                    border = TRUE
                          )
  anno_gtdb_row = rowAnnotation(`GTDB Species` = species_names,
                                    col = list(`GTDB Species` = make_col_list(species_names, "Paired")),
                                    border = TRUE
                          )
  anno_magotu_col = HeatmapAnnotation(`magOTU` = magOTU_names,
                                col = list(`magOTU` = make_col_list(magOTU_names, "Spectral")
                                          ),
                                          border = TRUE
                          )
  anno_magotu_row = rowAnnotation(`magOTU` = magOTU_names,
                                col = list(`magOTU` = make_col_list(magOTU_names, "Spectral")
                                          ),
                                          border = TRUE
                          )
  anno_tax_col = HeatmapAnnotation(`GTDB Species` = species_names,
                                   `magOTU` = magOTU_names,
                                   col = list(`magOTU` = make_col_list(magOTU_names, "Spectral"),
                                              `GTDB Species` = make_col_list(species_names, "Paired")
                                             ),
                                             border = TRUE
                                   )
  anno_tax_row = rowAnnotation(`GTDB Species` = species_names,
                                   `magOTU` = magOTU_names,
                                   col = list(`magOTU` = make_col_list(magOTU_names, "Spectral"),
                                              `GTDB Species` = make_col_list(species_names, "Paired")
                                             ),
                                             border = TRUE
                                   )
  if (row_split == "magotus") {
    row_split = as.character(as.factor(magOTU_names))
    # column_split = as.numeric(as.factor(host_names))
  }
  if (add_values) {
    cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(ifelse(matrix[i,j] > 0.7, sprintf("%.2f", matrix[i, j]), "-"), x, y, gp = gpar(fontsize = value_size))
          }
    heatmap_obj = Heatmap(matrix, 
            col = col_fun,
            # col = brewer.pal(2, "Pastel1"),
            # clustering_method_columns = 'ward.D2',
            top_annotation = anno_host_col,
            right_annotation = anno_gtdb_row,
            left_annotation = anno_magotu_row,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 0),
            heatmap_legend_param = list(title = "ANI", color_bar = "Continuous"),
            row_split = row_split,
            row_title_rot = 0,
            cell_fun = cell_fun,
            cluster_columns = F,
            border = TRUE,
            # row_dend_reorder = T,
            cluster_rows = F
            )
  } else {
    heatmap_obj = Heatmap(matrix, 
            col = col_fun,
            # col = brewer.pal(2, "Pastel1"),
            # clustering_method_columns = 'ward.D2',
            top_annotation = anno_host_col,
            right_annotation = anno_gtdb_row,
            left_annotation = anno_magotu_row,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 0),
            heatmap_legend_param = list(title = "ANI", color_bar = "Continuous"),
            row_split = row_split,
            row_title_rot = 0,
            # cell_fun = cell_fun,
            cluster_columns = F,
            border = TRUE,
            # row_dend_reorder = T,
            cluster_rows = F
            )
  }
  return(heatmap_obj)
  # draw(heatmap_obj, merge_legend = TRUE)
}



# #d53e4f, #f46d43, #fee08b, #66c2a5, #3288bd
my_col_fun = colorRamp2(c(0.7, 0.8, 0.9, 0.93, 0.94, 0.95, 1), c("#c6dbef", "#9ecae1", "#2171b5", "#4d9221", "#fddbc7", "#b2182b", "#67000d"))

get_ani_matrix <- function(genera_list = c(), mag_names = c(), magOTUs = c()) {
  if (length(genera_list) > 0) {
    matrix <- drep_N %>%
                  select(reference, querry, ani) %>%
                  filter(reference %in% vis_magOTUs_df$ID) %>%
                    right_join(vis_magOTUs_df %>% ungroup() %>% select(ID, Cluster, Host, Species, Genus), by = c("reference" = "ID")) %>%
                      filter(Genus %in% genera_list) %>%
                      # mutate(reference = paste0(Cluster, "_", Genus, "_", reference)) %>%
                      # mutate(querry = paste0(Cluster, "_", Genus, "_", querry)) %>%
                      arrange(Genus, Cluster, Host, Species) %>%
                      select(!c(Genus, Cluster)) %>%
                        pivot_wider(names_from=querry, values_from=ani) %>%
                          replace(is.na(.), 0.7) %>%
                            column_to_rownames("reference") %>%
                            select(rownames(.)) %>%
                                as.matrix
  }
  if (length(mag_names) > 0) {
    matrix <- drep_N %>%
                  select(reference, querry, ani) %>%
                  filter(reference %in% vis_magOTUs_df$ID) %>%
                    right_join(vis_magOTUs_df %>% ungroup() %>% select(ID, Cluster, Host, Species, Genus), by = c("reference" = "ID")) %>%
                      filter(reference %in% mag_names) %>%
                      # filter(reference %in% mag_names | querry %in% mag_names) %>%
                      # mutate(reference = paste0(Cluster, "_", Genus, "_", reference)) %>%
                      # mutate(querry = paste0(Cluster, "_", Genus, "_", querry)) %>%
                      arrange(Genus, Cluster, Host, Species) %>%
                      select(!c(Genus, Cluster)) %>%
                        pivot_wider(names_from=querry, values_from=ani) %>%
                          replace(is.na(.), 0.7) %>%
                            column_to_rownames("reference") %>%
                            select(rownames(.)) %>%
                                as.matrix  
  }
  if (length(magOTUs) > 0) {
    matrix <- drep_N %>%
                  select(reference, querry, ani) %>%
                  filter(reference %in% vis_magOTUs_df$ID) %>%
                    right_join(vis_magOTUs_df %>% ungroup() %>% select(ID, Cluster, Host, Genus), by = c("reference" = "ID")) %>%
                      filter(Cluster %in% genera_list) %>%
                      # mutate(reference = paste0(Cluster, "_", Genus, "_", reference)) %>%
                      # mutate(querry = paste0(Cluster, "_", Genus, "_", querry)) %>%
                      arrange(Genus, Cluster, Host, Species) %>%
                      select(!c(Genus, Cluster)) %>%
                        pivot_wider(names_from=querry, values_from=ani) %>%
                          replace(is.na(.), 0.7) %>%
                            column_to_rownames("reference") %>%
                            select(rownames(.)) %>%
                                as.matrix  
  }
  return(matrix)
}

# example of magOTU seperation, true species separation?
draw(plot_ani_heatmap(get_ani_matrix(mag_names = c("MAG_C3.3_9", "MAG_C1.1_13", "MAG_C1.5_6", 
                                              "MAG_C3.5_1", "MAG_C2.4_10", "MAG_C3.1_2", "MAG_C3.2_7")), col_fun = my_col_fun, add_values = T),
                                              merge_legend = TRUE)

draw(plot_ani_heatmap(get_ani_matrix(genera_list = c("g__Bifidobacterium")), col_fun = my_col_fun),
                                              merge_legend = TRUE)
draw(plot_ani_heatmap(get_ani_matrix(genera_list = c("g__Lactobacillus")), col_fun = my_col_fun),
                                              merge_legend = TRUE)

all_genera <- vis_magOTUs_df %>%
  pull(Genus) %>% unique()

for (genus_iter in all_genera) {
  pdf(paste0("Figures/06-ANI_heatmaps/06-", genus_iter, "_ANI_heatmap.pdf"))
  draw(plot_ani_heatmap(get_ani_matrix(genera_list = c(genus_iter)), col_fun = my_col_fun),
      merge_legend = T)
  dev.off()
  pdf(paste0("Figures/06-ANI_heatmaps/06-", genus_iter, "_ANI_heatmap_values.pdf"))
  draw(plot_ani_heatmap(get_ani_matrix(genera_list = c(genus_iter)), col_fun = my_col_fun, add_values = T, value_size = 2),
      merge_legend = T)
  dev.off()
}


vis_magOTUs_df %>%
  filter(Cluster =="33_1") %>%
    pull(N50)
