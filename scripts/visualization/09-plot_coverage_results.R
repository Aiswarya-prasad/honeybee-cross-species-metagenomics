#*********Work in progress*********#

##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

shorten_list <- function(covs_list_full, start = 1, end = 2) {
  covs_list_short = strsplit(covs_list_full, ",")[[1]][start:end]
  return(paste0(covs_list_short, collapse = ","))
}

trim_side_mean <- function(x, trim, type="both") {
    if (type == "both") {
        mean(x,trim)}
    else if (type == "right") {
        x <- sort(x)
        mean(x[1:floor(length(x)*(1-trim))])}
    else if (type == "left"){
        x <- sort(x)
        mean(x[max(1,floor(length(x)*trim)):length(x)])}
}

make_cum_curve <- function(samples_vector, pa_df, iterations, name = NA) {
  num_clusters_matrix <- matrix(nrow = length(samples_vector), ncol = iterations)
  for (iter in 1:iterations) {
    clusters_found_cumulative <- c()
    for (num_samples in 1:length(samples_vector)) {
      num_new_clusters = 0
      selected_samples <- sample(samples_vector, num_samples)
      clusters_found <- colnames(pa_df[selected_samples, which(colSums(pa_df[selected_samples, ]) > 1)])
      for (cluster in clusters_found) {
        if (cluster %in% clusters_found_cumulative) {
          invisible()
        } else {
          num_new_clusters <- num_new_clusters + 1
          clusters_found_cumulative <- c(clusters_found_cumulative, cluster)
        }
      }
      num_clusters_matrix[num_samples, iter] = length(clusters_found_cumulative)
    }
  }
  num_clusters_df <- as.data.frame(num_clusters_matrix)
  colnames(num_clusters_df) <- do.call(function(x) paste0("curve_", x), list(c(1:iterations)))
  num_clusters_df <- cbind(sample_size = c(1:length(samples_vector)), num_clusters_df)
  plot_cum_curve_df <- pivot_longer(num_clusters_df, !sample_size, values_to = "number_of_clusters", names_to = "curve")
  plot_cum_curve_df <- cbind("name" = name, plot_cum_curve_df)
  return(plot_cum_curve_df)
}
pcoa_plot <- function(df_pcoa, metadata, variable, color_add=F, color_list, colname_in_metadata = "ID") {
          matrix <- as.matrix(df_pcoa)
          dist <- as.dist(matrix)
          res_pcoa <- pcoa(dist)
          ev1 <- res_pcoa$vectors[,1]
          ev2 <- res_pcoa$vectors[,2]
          df_pcoa_new <- data.frame(cbind(ev1,ev2))
          df_pcoa_new$Sample <- rownames(df_pcoa_new)
          rownames(df_pcoa_new) <- NULL
          df_pcoa_new <- left_join(df_pcoa_new, metadata, by = c("Sample" = colname_in_metadata))
          perc_axis <- round(((res_pcoa$values$Relative_eig[c(1,2)])*100), digits=1)
          axis_x_title <- paste0("PCo1 (",perc_axis[1],"%)")
          axis_y_title <- paste0("PCo2 (",perc_axis[2],"%)")
          if(color_add) {
            p <- ggplot(df_pcoa_new, aes(x = ev1,
                                       y = ev2,
                                       colour = get(variable)))+
                geom_point(stat="identity", size=4, shape=19) +
                  labs(x=axis_x_title, y = axis_y_title, color = variable) +
                    make_theme(setFill = F, setCol = F, guide_nrow = 4, leg_size = 10 ) +
                      scale_color_manual(values=color_list)
          } else {
            p <- ggplot(df_pcoa_new, aes(x = ev1,
                                       y = ev2,
                                       color = get(variable)))+
                geom_point(stat="identity", size=4, shape=19) +
                  labs(x=axis_x_title, y = axis_y_title, color = variable) +
                    make_theme( max_colors = length(unique(df_pcoa_new[, variable])), guide_nrow = 4, leg_size = 10 ) 
          }
          return(p)
}

##############
# files to be read
##############

source('scripts/visualization/01-plot_mapping_output.R', chdir = TRUE)
source('scripts/visualization/05-read_MAG_metadata.R', chdir = TRUE)
corecov_data_dir <- "09_MagDatabaseProfiling/CoverageEstimation/Merged"
coords_df <- data.frame()
for (file in list.files(corecov_data_dir)) {
  if (endsWith(file, "_coord.txt")) {
    df_temp <- read.csv(paste0(corecov_data_dir, "/", file), sep = "\t", header = TRUE)
    coords_df <- rbind(coords_df, df_temp)
  }
}

##############
# analyse data and plot
##############

cluster_tax_info <- vis_magOTUs_df %>%
                      ungroup() %>%
                      select(Genus, Species, Cluster, Num_mags) %>%
                      unique() %>%
                        mutate(cluster_name = paste0(Genus, "-", Cluster))

coverage_df <- coords_df %>%
  rename(Cluster = magOTU) %>%
    left_join(cluster_tax_info) %>%
      mutate(Host = Vectorize(get_host_name)(Sample)) %>%
      mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
        mutate(Detected = ifelse(Cov > 0, 1, 0)) %>%
          mutate(Detected = ifelse(is.na(Detected), 0, Detected)) %>%
          mutate(Cov = ifelse(is.na(Cov), 0, Cov)) %>%
          group_by(Host, Cluster) %>%
            mutate(Present = sum(Detected)) %>%
            mutate(Mean = mean(Cov)) %>%
            mutate(Prevalence_host = ifelse(Host=="Apis mellifera", Present/length(samples_am), NA)) %>%
            mutate(Prevalence_host = ifelse(Host=="Apis cerana", Present/length(samples_ac), Prevalence_host)) %>%
            mutate(Prevalence_host = ifelse(Host=="Apis dorsata", Present/length(samples_ad), Prevalence_host)) %>%
            mutate(Prevalence_host = ifelse(Host=="Apis florea", Present/length(samples_af), Prevalence_host)) %>%
            mutate(Prevalence_host = ifelse(Host=="Apis andreniformis", Present/length(samples_aa), Prevalence_host))

coverage_df_genus <- coords_df %>%
                rename(Cluster = magOTU) %>%
                  left_join(cluster_tax_info) %>%
                      mutate(Detected = ifelse(Cov > 1, 1, 0)) %>%
                        mutate(Detected = ifelse(is.na(Detected), 0, Detected)) %>%
                        mutate(Cov = ifelse(is.na(Cov), 0, Cov)) %>%
                        ungroup() %>%
                          group_by(Genus, Sample) %>%
                              summarise(Genus_Present = ifelse(any(Detected == 1), 1, 0), Num_mags = sum(Num_mags), Genus_cov = mean(Cov)) %>%
                               mutate(Host = Vectorize(get_host_name)(Sample)) %>%
                                group_by(Host, Genus) %>%
                                  mutate(Mean = round(trim_side_mean(Genus_cov, trim = 0.01)), 2) %>%
                                  mutate(Median = round(mean(Genus_cov)), 2) %>%
                                  mutate(covs_list = paste0(round(Genus_cov, 0), collapse = ",")) %>%
                                  mutate(covs_list = Vectorize(shorten_list)(covs_list, 1, 3)) %>%
                                  mutate(Present = sum(Genus_Present)) %>%
                                  mutate(Prevalence_host = ifelse(Host=="Apis mellifera", Present/length(samples_am), NA)) %>%
                                  mutate(Prevalence_host = ifelse(Host=="Apis cerana",Present/length(samples_ac), Prevalence_host)) %>%
                                  mutate(Prevalence_host = ifelse(Host=="Apis dorsata",Present/length(samples_ad), Prevalence_host)) %>%
                                  mutate(Prevalence_host = ifelse(Host=="Apis florea",Present/length(samples_af), Prevalence_host)) %>%
                                  mutate(Prevalence_host = ifelse(Host=="Apis andreniformis",Present/length(samples_aa), Prevalence_host))

ggplot(coverage_df, aes(y = factor(cluster_name), x = factor(Sample, samples_IN_MY), fill = Cov)) +
                            geom_tile() +
                              labs(x = "Sample", y = "Cluster")+
                              scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#1a88c9", guide = "colourbar", trans = "log10") +
                              # scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#ffc20e", guide = "colourbar", trans = "log10") +
                                make_theme(setFill = F, modify_guide = F, y_size = 8,
                                  x_size = 0, x_angle = 30, x_hj = 1 , x_vj = 1,
                                  y_hj =1, leg_pos = "right"
                                )
ggplot(coverage_df, aes(y = factor(cluster_name), x = factor(Host, host_order), fill = Prevalence_host)) +
                            geom_tile() +
                              labs(x = "Sample", y = "magOTU", fill = "Prevalence\n(Present = \nCov > 1)\n")+
                              scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#1a88c9", guide = "colourbar") +
                              # scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#ffc20e", guide = "colourbar") +
                                coord_fixed(0.1) +
                                make_theme(setFill = F, modify_guide = F, y_size = 8,
                                  y_hj =1, leg_pos = "right"
                                )
                                ggsave("Figures/09-Prevalence_magOTU_by_host.pdf")

ggplot(coverage_df_genus %>% filter(Genus != "g__"), aes(y = factor(Genus, genera), x = factor(Host, host_order), fill = Prevalence_host)) +
                            geom_tile() +
                            geom_text(aes(label = ifelse(Mean > 0, Mean, ""))) +
                              labs(x = "Sample", y = "Genus", fill = "Prevalence\n(Present = \nCov > 1)\nTruncated mean\n(0.01)")+
                              # scale_fill_gradient(na.value = "transparent", low = "#fff5f0" , high = "#99000d", guide = "colourbar") +
                              scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#1a88c9", guide = "colourbar") +
                                guides(fill = guide_colourbar(nbin = 8)) +
                                make_theme(setFill = F, modify_guide = F, y_size = 10,
                                  x_hj= 0.5, x_size = 10,
                                  y_hj =1, leg_pos = "right", leg_size = 10
                                )

ggplot(coverage_df %>% filter(Genus != "g__") %>% mutate(Cov = round(Cov, 4)),
        aes(y = factor(Genus, genera), x = Cov,, color =  factor(Genus, genera))) +
                            geom_jitter() +
                            # geom_text(aes(label = Mean)) +
                              labs(x = "Coverage", y = "Cluster")+
                              # scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#ffc20e", guide = "colourbar") +
                              scale_color_manual(values = genusColors) +
                                # guides(fill = guide_colourbar(nbin = 8)) +
                                facet_wrap(~ factor(Host, host_order)) +
                                # scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
                                make_theme(setCol = F, modify_guide = F, y_size = 10,
                                  x_hj= 0.5, x_size = 10,
                                  y_hj =1, leg_pos = "right", leg_size = 10
                                )

ggplot(coverage_df %>% filter(Genus != "g__") %>% mutate(Cov = round(Cov, 4)),
        aes(y = factor(Genus, genera), x = Cov,, color =  factor(Genus, genera))) +
                            geom_point() +
                            # geom_text(aes(label = Mean)) +
                              labs(x = "Coverage", y = "Cluster")+
                              # scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#ffc20e", guide = "colourbar") +
                              scale_color_manual(values = genusColors) +
                                # guides(fill = guide_colourbar(nbin = 8)) +
                                facet_wrap(~ factor(Host, host_order)) +
                                scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
                                make_theme(setCol = F, modify_guide = F, y_size = 10,
                                  x_hj= 0.5, x_size = 10,
                                  y_hj =1, leg_pos = "right", leg_size = 10
                                )
ggplot(coverage_df %>% filter(Genus != "g__") %>% filter(Sample %in% samples_IN_MY) %>% mutate(Cov = round(Cov, 4)) %>% filter(Cov != 0),
        aes(y = factor(Genus, genera), x = Cov, shape = Location, color =  Location)) +
        # aes(y = factor(Genus, genera), x = Cov, shape = Location, color =  factor(Genus, genera))) +
                            geom_jitter() +
                            # geom_text(aes(label = Mean)) +
                              labs(x = "Coverage", y = "Cluster")+
                              # scale_fill_gradient(na.value = "transparent", low = "#ffffff" , high = "#ffc20e", guide = "colourbar") +
                            #   scale_color_manual(values = genusColors) +
                                # guides(fill = guide_colourbar(nbin = 8)) +
                                facet_wrap(~ factor(Host, host_order)) +
                                scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
                                make_theme(setCol = F, modify_guide = F, y_size = 10,
                                  x_hj= 0.5, x_size = 10,
                                  y_hj =1, leg_pos = "right", leg_size = 10
                                )
                                ggsave(paste0("Figures/", "09-Coverage_scatter_location_compared.pdf"))

abundance_df <- coords_df %>%
  rename(Cluster = magOTU) %>%
    left_join(cluster_tax_info %>% filter(!is.na(cluster_name))) %>%
      mutate(Host = Vectorize(get_host_name)(Sample)) %>%
        mutate(Cov = ifelse(is.na(Cov), 0, Cov)) %>%
          mutate(MedianCov = ifelse(is.na(MedianCov), 0, MedianCov)) %>%
            select(Sample, Cov, cluster_name, Host) %>%
              left_join(select(df_reads, Trimmed, Sample), by = "Sample") %>%
                  rename(NumReads = Trimmed) %>%
                    mutate(CovNorm = Cov/NumReads)

abundance_df_matrix <-  abundance_df %>% filter(!is.na(cluster_name)) %>%
                          mutate(PresAb = ifelse(CovNorm > 0, 1, 0)) %>%
                          select(Sample, cluster_name, PresAb) %>%
                            pivot_wider(names_from=cluster_name, values_from=PresAb, values_fn = as.vector) %>% remove_rownames %>%
                            column_to_rownames(var="Sample")

abundance_df_matrix_pa <-  abundance_df %>% filter(!is.na(cluster_name)) %>%
                          mutate(PresAb = ifelse(Cov > 1, 1, 0)) %>%
                          select(Sample, cluster_name, PresAb) %>%
                            pivot_wider(names_from=cluster_name, values_from=PresAb, values_fn = as.vector) %>% remove_rownames %>%
                            column_to_rownames(var="Sample")

abundance_df_matrix_Cov_rel <-  abundance_df %>% filter(!is.na(cluster_name)) %>%
                          select(Sample, cluster_name, Cov) %>%
                            pivot_wider(names_from=cluster_name, values_from=Cov, values_fn = as.vector) %>% remove_rownames %>%
                            column_to_rownames(var="Sample") %>%
                              mutate(across(-1)*100/rowSums(across(-1)))

abundance_df_matrix_Cov <-  abundance_df %>% filter(!is.na(cluster_name)) %>%
                          select(Sample, cluster_name, Cov) %>%
                          filter(Sample %in% samples_IN) %>%
                            pivot_wider(names_from=cluster_name, values_from=Cov, values_fn = as.vector) %>% remove_rownames %>%
                            column_to_rownames(var="Sample")

column_ha = HeatmapAnnotation(Host = rownames(abundance_df_matrix_pa[samples_IN,]) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color)
                          )
column_ba = HeatmapAnnotation(Location = rownames(abundance_df_matrix_pa[samples_IN,]) %>% as.data.frame() %>% mutate(Location = Vectorize(get_location_from_sample_name)(.)) %>% pull(Location),
                              col = list(Location = location_country_colors)
                          )
row_ha = rowAnnotation(Genus = colnames(abundance_df_matrix_pa[samples_IN,]) %>% as.data.frame() %>% left_join(cluster_tax_info %>% ungroup() %>% select(Genus, cluster_name), by = c(. = "cluster_name")) %>% pull(Genus),
                       col = list(Genus = unlist(genusColors))
                    )
pdf("Figures/09-presence_absence_heatmap_IN.pdf")

Heatmap(t(abundance_df_matrix_pa[samples_IN,]), 
        col = c("#ffffff", "#1a88c9"),
        # col = brewer.pal(2, "Pastel1"),
        clustering_method_columns = 'ward.D2',
        top_annotation = column_ha, right_annotation = row_ha,
        bottom_annotation = column_ba,
        column_names_gp = grid::gpar(fontsize = 0),
        heatmap_legend_param = list(title = "Present / Absent", color_bar = "discrete"),
          cluster_rows = F)
dev.off()
column_ha = HeatmapAnnotation(Host = rownames(abundance_df_matrix_pa[samples_IN_MY,]) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color)
                          )
column_ba = HeatmapAnnotation(Location = rownames(abundance_df_matrix_pa[samples_IN_MY,]) %>% as.data.frame() %>% mutate(Location = Vectorize(get_location_from_sample_name)(.)) %>% pull(Location),
                              col = list(Location = location_country_colors)
                          )
row_ha = rowAnnotation(Genus = colnames(abundance_df_matrix_pa[samples_IN_MY,]) %>% as.data.frame() %>% left_join(cluster_tax_info %>% ungroup() %>% select(Genus, cluster_name), by = c(. = "cluster_name")) %>% pull(Genus),
                       col = list(Genus = unlist(genusColors))
                    )
pdf("Figures/09-presence_absence_heatmap_IN.pdf")

Heatmap(t(abundance_df_matrix_pa[samples_IN_MY,]), 
        col = c("#ffffff", "#1a88c9"),
        # col = brewer.pal(2, "Pastel1"),
        clustering_method_columns = 'ward.D2',
        top_annotation = column_ha, right_annotation = row_ha,
        bottom_annotation = column_ba,
        column_names_gp = grid::gpar(fontsize = 0),
        heatmap_legend_param = list(title = "Present / Absent", color_bar = "discrete"),
          cluster_rows = F)
dev.off()

column_ha = HeatmapAnnotation(Host = rownames(abundance_df_matrix_pa) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color)
                          )
column_ba = HeatmapAnnotation(Location = rownames(abundance_df_matrix_pa) %>% as.data.frame() %>% mutate(Location = Vectorize(get_location_from_sample_name)(.)) %>% pull(Location),
                              col = list(Location = location_country_colors)
                          )
row_ha = rowAnnotation(Genus = colnames(abundance_df_matrix_pa) %>% as.data.frame() %>% left_join(cluster_tax_info %>% ungroup() %>% select(Genus, cluster_name), by = c(. = "cluster_name")) %>% pull(Genus),
                       col = list(Genus = unlist(genusColors))
                    )
pdf("Figures/09-presence_absence_heatmap.pdf")

Heatmap(t(abundance_df_matrix_pa), 
        col = c("#ffffff", "#1a88c9"),
        # col = brewer.pal(2, "Pastel1"),
        clustering_method_columns = 'ward.D2',
        top_annotation = column_ha, right_annotation = row_ha,
        bottom_annotation = column_ba,
        column_names_gp = grid::gpar(fontsize = 0),
        heatmap_legend_param = list(title = "Present / Absent", color_bar = "discrete"),
          cluster_rows = F)
dev.off()
# pdf("Figures/09-relative_abundance_heatmap.pdf")
Heatmap(t(abundance_df_matrix_Cov_rel),
        col = colorRamp2(c(0, 10, 100), colors = c("#ffffff", "#1a88c9", "#1a3869")),
        # col = colorRamp2(c(0, 10, 100), colors = c(brewer.pal(9, "BrBG")[9], brewer.pal(9, "BrBG")[5], brewer.pal(9, "BrBG")[1])),
        clustering_method_columns = 'ward.D2',
        top_annotation = column_ha, right_annotation = row_ha,
          cluster_rows = F)
# dev.off()
column_ha = HeatmapAnnotation(Host = rownames(abundance_df_matrix_Cov) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color)
                          )
row_ha = rowAnnotation(Genus = colnames(abundance_df_matrix_Cov) %>% as.data.frame() %>% left_join(cluster_tax_info %>% ungroup() %>% select(Genus, cluster_name), by = c(. = "cluster_name")) %>% pull(Genus),
                       col = list(Genus = unlist(genusColors))
                    )
Heatmap(t(abundance_df_matrix_Cov),
        col = colorRamp2(c(0, 100, 1000), colors = c("#ffffff", "#1a88c9", "#08306b")),
        # col = colorRamp2(c(0, 10, 100), colors = c(brewer.pal(9, "BrBG")[9], brewer.pal(9, "BrBG")[5], brewer.pal(9, "BrBG")[1])),
        clustering_method_columns = 'ward.D2',
        top_annotation = column_ha, right_annotation = row_ha,
          cluster_rows = F)

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
        filter(Location %in% c("Malaysia", "India")) %>%
          left_join(as.data.frame(rowSums(data_otu)) %>% setNames("cov_sums") %>% rownames_to_column("Sample"), by = "Sample"),
       aes(x = MAGs_DB, y = cov_sums, 
           color = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads", y = "Sum of all coverages", color = "Host") +
    geom_point(size = 3) +
      make_theme(setCol = F) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))

ggplot(df_reads %>% filter((Sample %in% samples)) %>% 
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample))
        ) +
    geom_point(aes(x = factor(Sample, samples),
                   y = MAGs_DB,
                   color = Location,
                   shape = factor(Host, host_order))) +
    geom_hline(yintercept = 1e+06) +
    # scale_color_manual(values=host_order_color_dark) +
    scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
    make_theme(setCol = F)

ggplot(df_reads %>% filter((Sample %in% samples)) %>% 
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample))
        ) +
    geom_point(aes(x = factor(Sample, samples),
                   y = MAGs_DB,
                   shape = Location,
                   color = factor(Host, host_order))) +
    geom_hline(yintercept = 1e+06) +
    scale_color_manual(values=host_order_color_dark) +
    scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
    make_theme(setCol = F)


ggplot(df_reads %>%
        select(Sample, MAGs_DB) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% ungroup() %>% select(Sample, Cov, Genus, cluster_name) %>% filter(Cov > 0) %>% filter(!is.na(cluster_name)) %>% group_by(Sample, Genus) %>%  mutate(Sample, Genus, genus_Cov = sum(Cov))),
       aes(x = MAGs_DB, y = genus_Cov, 
        #    color = Location)
        #    color = factor(Genus, genera), shape = Location)
           color = factor(Host, host_order), shape = Location)
    ) +
    # facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads", y = "Sum of all coverages", color = "Host") +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(setCol = F, guide_nrow = 5,
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Lactobacillus")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Lactobacillus coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Lactobacillus.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Bombilactobacillus")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Bombilactobacillus coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Bombilactobacillus.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Bifidobacterium")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Bifidobacterium coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Bifidobacterium.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Gilliamella")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Gilliamella coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Gilliamella.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Snodgrassella")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Snodgrassella coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Snodgrassella.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Apibacter")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Apibacter coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Apibacter.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Frischella")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Frischella coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Frischella.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Dysgonomonas")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Dysgonomonas coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          scale_size(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Dysgonomonas.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Commensalibacter")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Commensalibacter coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          scale_size(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Commensalibacter.pdf")


ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Pectinatus")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Pectinatus coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          scale_size(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Pectinatus.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Bartonella")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Bartonella coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          scale_size(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Bartonella.pdf")
# coverage_df %>% filter(grepl("kun", Species)) %>% pull(Species)

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__Bombella")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__Bombella coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          scale_size(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__Bombella.pdf")

ggplot(df_reads %>%
        select(Sample, MAGs_DB, Trimmed) %>%
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
          left_join(coverage_df %>% 
                    ungroup() %>%
                    select(Sample, Cov, Genus, cluster_name) %>%
                    filter(Cov > 0)%>%
                    filter(!is.na(cluster_name)) %>%
                    group_by(Sample, Genus) %>% 
                    mutate(Sample, Genus, genus_Cov = sum(Cov)) %>%
                    filter(Genus == "g__WRHT01")
                    ),
       aes(x = MAGs_DB, y = genus_Cov, size = Trimmed,
           color = Location)
        #    color = factor(Genus, genera), shape = Location)
        #    color = factor(Host, host_order), shape = Location)
    ) +
    facet_wrap(~factor(Host, host_order)) +
    labs(x = "Numer of MAGs DB mapped reads",
         y = "Sum of all coverages",
         color = "Host",
         size = "Number of trimmed reads",
         title = "g__WRHT01 coverage summed") +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = 5, linetype = "dashed") +
    geom_hline(yintercept = 10, linetype = "dotted") +
      make_theme(guide_nrow = 3, palettecolor = "Set1",
                 x_angle = 0, x_size = 7) +
        # scale_color_manual(values=genusColors) +
        # scale_color_manual(values=host_order_color_dark) +
          scale_y_continuous(trans = "log10") +
          scale_size(labels=unit_format(unit = "M", scale = 1e-6)) +
          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
          scale_size(labels=unit_format(unit = "M", scale = 1e-6))
          ggsave("Figures/09-Coverage_vs_mapped_reads_g__WRHT01.pdf")

deep_samples <- df_reads %>% filter(Trimmed > 2e+6) %>% pull(Sample)
length(deep_samples)
df_magOTUs_vegan <- abundance_df_matrix_pa %>%
                      filter(rowSums(.) > 0) %>%
                        filter(row.names(.) %in% deep_samples)
dist_matrix <- as.matrix(vegdist(df_magOTUs_vegan, "bray"))

pcoa_magotus <- pcoa_plot(dist_matrix, df_meta_complete, "Species", color_add=T, color_list=host_order_color_dark, colname_in_metadata = "ID")
          # ggsave("Figures/09d-magOTUs_pcoa.pdf")
pcoa_magotus
# code for betapart functions from
# https://github.com/cran/betapart

# Latest commit a222c6d on Feb 24, 2012
betapart.core <- function(x){
	
	# test for a binary numeric matrix or data.frame
	if(! is.matrix(x)){
		x<-as.matrix(x)
  	}
	
	if(! is.numeric(x))
    	stop("The data in x is not numeric.",call.=TRUE)
	
	# simple test for binary data
	xvals <-  unique(as.vector(x))
	if (any(!is.element(xvals, c(0, 1)))) 
        stop("The table contains values other than 0 and 1: data should be presence/absence.", call. = TRUE)

	shared <- x %*% t(x)
    not.shared <-  abs(sweep(shared, 2, diag(shared)))
		
	sumSi <- sum(diag(shared)) # species by site richness
    St <- sum(colSums(x) > 0)  # regional species richness
    a <- sumSi - St            # multi site shared species term

    sum.not.shared <- not.shared + t(not.shared)
    max.not.shared <- pmax(not.shared, t(not.shared))
    min.not.shared <- pmin(not.shared, t(not.shared))

	computations<-list(data=x, sumSi=sumSi, St=St, a=a, shared=shared, not.shared=not.shared, 
	                   sum.not.shared=sum.not.shared, max.not.shared=max.not.shared, 
	                   min.not.shared=min.not.shared)
    class(computations)<-"betapart"

	return(computations)
} 
# Latest commit 4b59d9c on Feb 24, 2012
beta.pair <- function(x, index.family="sorensen"){
	
	# test for a valid index
	index.family <- match.arg(index.family, c('jaccard','sorensen'))
	
	# test for pre-existing betapart objects
	if (! inherits(x, "betapart")){	
		x <- betapart.core(x)
	}
	
	# run the analysis given the index
	switch(index.family,
		sorensen = {
			beta.sim <- x$min.not.shared / (x$min.not.shared + x$shared)
			beta.sne <- ((x$max.not.shared - x$min.not.shared) / ((2 * x$shared) + x$sum.not.shared)) * (x$shared / (x$min.not.shared + x$shared))
            beta.sor <- x$sum.not.shared / (2 * x$shared + x$sum.not.shared)
                
            pairwise <- list(beta.sim=as.dist(beta.sim), beta.sne=as.dist(beta.sne),beta.sor=as.dist(beta.sor))},
		jaccard = {
			beta.jtu <- (2*x$min.not.shared) / ((2*x$min.not.shared) + x$shared)
	        beta.jne <- ((x$max.not.shared - x$min.not.shared) / (x$shared + x$sum.not.shared)) * (x$shared / ((2*x$min.not.shared) + x$shared))
	        beta.jac <- x$sum.not.shared / (x$shared + x$sum.not.shared)

	        pairwise <- list(beta.jtu=as.dist(beta.jtu), beta.jne=as.dist(beta.jne),beta.jac=as.dist(beta.jac))})

	return(pairwise)
}
beta.multi <- function(x, index.family="sorensen"){

	# test for a valid index
	index.family <- match.arg(index.family, c('jaccard','sorensen'))
	
	# test for pre-existing betapart objects
	if (! inherits(x, "betapart")){	
		x <- betapart.core(x)
	}

	maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
    minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])
	
	# run the analysis given the index
	switch(index.family,
		sorensen = {
            beta.sim <- minbibj / (minbibj + x$a)
            beta.sne <- (x$a / (minbibj + x$a)) * ((maxbibj - minbibj) / ((2 * x$a) + maxbibj + minbibj))
            beta.sor <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * x$a))

           	multi <- list(beta.SIM=beta.sim, beta.SNE=beta.sne,beta.SOR=beta.sor)},
		jaccard = {
            beta.jtu <- (2*minbibj) / ((2*minbibj) + x$a)
            beta.jne <- (x$a / ((2*minbibj) + x$a)) * ((maxbibj - minbibj) / ((x$a) + maxbibj + minbibj))
            beta.jac <- (minbibj + maxbibj) / (minbibj + maxbibj + x$a)

           	multi <- list(beta.JTU=beta.jtu, beta.JNE=beta.jne, beta.JAC=beta.jac)})

	return(multi)

}
str()
dist_matrix_betapart <- beta.pair(df_magOTUs_vegan, "jaccard")
dist_matrix_betapart <- beta.pair(df_magOTUs_vegan, "sorensen")
str(dist_matrix_betapart)
pcoa_plot(dist_matrix_betapart$beta.sor, df_meta_complete, "Species", color_add=T, color_list=host_order_color_dark)
pcoa_plot(dist_matrix_betapart$beta.sim, df_meta_complete, "Species", color_add=T, color_list=host_order_color_dark)

pcoa_plot(dist_matrix_betapart$beta.jac, df_meta_complete, "Species", color_add=T, color_list=host_order_color_dark)
pcoa_plot(dist_matrix_betapart$beta.jtu, df_meta_complete, "Species", color_add=T, color_list=host_order_color_dark)

pcoa_plot(dist_matrix_betapart$beta.jac, df_meta_complete, "Country", color_add=F, color_list=host_order_color_dark)
pcoa_plot(dist_matrix_betapart$beta.jtu, df_meta_complete, "Country", color_add=F, color_list=host_order_color_dark)




df_plot_cum_curve_m <- make_cum_curve(samples_am, df_magOTUs_vegan, 100, "Apis mellifera")
df_plot_cum_curve_c <- make_cum_curve(samples_ac, df_magOTUs_vegan, 100, "Apis cerana")
df_plot_cum_curve_d <- make_cum_curve(samples_ad, df_magOTUs_vegan, 100, "Apis dorsata")
df_plot_cum_curve_f <- make_cum_curve(samples_af, df_magOTUs_vegan, 100, "Apis florea")
df_plot_cum_curve_a <- make_cum_curve(samples_aa, df_magOTUs_vegan, 100, "Apis andreniformis")
df_plot_cum_curve <- rbind(df_plot_cum_curve_m, df_plot_cum_curve_c, df_plot_cum_curve_d, df_plot_cum_curve_f, df_plot_cum_curve_a)
ggplot(data = df_plot_cum_curve, aes(x = sample_size, y = number_of_clusters, color = factor(name, host_order))) +
                      geom_jitter(position = position_dodge(width=0.7)) +
                        geom_smooth(se = FALSE) +
                          labs(color = "Host species", x = "# Bees", y = "Number of magOTUs") +
                          scale_color_manual(values=host_order_color) +
                            make_theme(leg_pos = "bottom", setCol = F, guide_nrow = 1)
