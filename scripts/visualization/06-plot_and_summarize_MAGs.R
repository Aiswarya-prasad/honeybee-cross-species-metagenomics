#*********Work in progress*********#

##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)
source('scripts/visualization/01-plot_mapping_output.R')
source('scripts/visualization/05-read_MAG_metadata.R')

##############
# files to be read
##############


##############
# analyse data and plot
##############

num_mags_df <- df_reads %>% select(Sample, Trimmed, MAGs_DB, Host_mapped) %>%
                left_join(vis_magOTUs_df %>% ungroup() %>% group_by(Sample) %>% summarise(Host, Number_genera = n_distinct(Genus), Number_of_mags = n())) %>%
                left_join(vis_magOTUs_df %>% ungroup() %>% group_by(Sample) %>% summarise(Host, Number_clusters = n_distinct(Cluster))) %>%
                filter(Sample %in% samples_IN)


ggplot() +
  geom_point(data = num_mags_df,
             aes(x = MAGs_DB,
                 y = Number_of_mags,
                 size = Number_genera,
                 color = Host)) +
                 labs(x = "Reads mapped to MAG DB",
                      y = "Number of medium quality MAGs",
                      size = "Number of Genera",
                      color = "Host species"
                 ) +
                 make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_color_manual(values = host_order_color_dark)
                                ggsave("Figures/06-number_of_MAGs_per_depth_mapped.pdf")
ggplot() +
  geom_point(data = num_mags_df,
             aes(x = Trimmed,
                 y = Number_of_mags,
                 size = Number_genera,
                 color = Host)) +
                 labs(x = "Reads after trimming",
                      y = "Number of medium quality MAGs",
                      size = "Number of Genera",
                      color = "Host species") +
                 make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_color_manual(values = host_order_color_dark)
                                ggsave("Figures/06-number_of_MAGs_per_depth.pdf")
ggplot() +
  geom_point(data = num_mags_df,
             aes(x = Trimmed,
                 y = Number_of_mags,
                 size = Number_clusters,
                 color = Host)) +
                 labs(x = "Reads after trimming",
                      y = "Number of medium quality MAGs",
                      size = "Number of Genera",
                      color = "Host species") +
                 make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 2) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_color_manual(values = host_order_color_dark)

ggplot() +
  geom_point(data = num_mags_df,
             aes(x = Host_mapped,
                 y = Number_of_mags,
                 color = Host), size = 4) +
                 make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_color_manual(values = host_order_color)

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

