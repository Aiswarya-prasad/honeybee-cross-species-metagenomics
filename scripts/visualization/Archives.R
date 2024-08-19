# cazy info from Cayman profiling output

```{r}
# visualize cazy
# not sure how this was parsed out - probably somewherte in python, worst case, interactively...
cazy_df <- read.csv("results/figures/cazy_df.tsv", sep = "\t") %>%
            filter(sample %in% samples) %>%
            mutate(Host = Vectorize(get_host_from_sample_name)(sample))
host_count <- cazy_df %>% group_by(Host) %>% summarize(count = n_distinct(sample))
get_prev <- function(present, host) {
  total <- pull(filter(host_count, Host == host), count)
  return(present / total)
}

cazy_df <- read.csv("results/figures/cazy_df.tsv", sep = "\t") %>%
            mutate(Host = Vectorize(get_host_from_sample_name)(sample)) %>%
              mutate(present = ifelse(uniq_rpkm > 0, 1, 0)) %>%
              group_by(feature, Host) %>%
              mutate(count = sum(present)) %>%
              ungroup() %>%
              mutate(Prevalence_host = Vectorize(get_prev)(count, Host))

# # count is number of samples of host in which it is present
# # present if uniq_rpkm > 0
# order_feature_count <- cazy_df %>%
#             group_by(feature) %>%
#             mutate(sum = sum(count)) %>%
#             arrange(sum) %>%
#             pull(feature) %>% unique

# order_feature_prev <- cazy_df %>% 
#             reframe(feature, Host, present, Prevalence_host) %>%
#             group_by(feature) %>%
#             mutate(sum = sum(Prevalence_host)) %>%
#             arrange(sum) %>%
#             pull(feature) %>% unique

# ggplot(cazy_df %>% reframe(feature, Host, count, Prevalence_host) %>% unique) +
#   geom_tile(aes(x = Host, y = feature, fill = Prevalence_host)) +
#   labs(y = "CAZy GH family", x = "Host", fill = "Prevalence") +
#   scale_y_discrete(limits = order_feature_prev) +
#   coord_fixed(ratio = 0.05) +
#   make_theme(setFill = FALSE, leg_pos = "right", 
#              x_angle = 45, x_size = 10,
#              y_size = 6, y_hj = 1)


# completed_substrate_annotations <- read_xlsx("data/cayman_gene_db/20230607_glycan_annotations_cleaned.xlsx")
completed_substrate_annotations <- read.csv("data/cayman_gene_db/20230607_glycan_annotations_cleaned_manually.tsv", header = T, sep = "\t")
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily,ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
glycan_annotations_final_cleaned_long <- glycan_annotations_final_cleaned %>% separate_rows(ORIGIN,sep=",") %>% separate_rows(FUNCTION_IN_ORIGIN,sep=",") %>% separate_rows(FUNCTION_AT_DESTINATION_1,sep=",") %>% separate_rows(FUNCTION_AT_DESTINATION_2,sep=",") %>% separate_rows(FUNCTION_AT_DESTINATION_3,sep=",")
glycan_annotations_cleaned <- glycan_annotations_final_cleaned

glycan_annotations_cleaned_filt <- glycan_annotations_cleaned %>% filter(!is.na(FUNCTION_AT_DESTINATION_3))

# and abbreviate long glycan names
cazy_df_annotated <- cazy_df %>%
                      left_join(glycan_annotations_cleaned, by = c("feature" = "Subfamily")) %>%
                      mutate(mod_names_pec = ifelse(grepl("pectin", FUNCTION_AT_DESTINATION_3), 1, 0)) %>%
                      mutate(short_names = abbreviate(FUNCTION_AT_DESTINATION_3, minlength = 50))

# # heatmap of Glycan_annotatoin vs host showing prevalence as size of dot
# ggplot(cazy_df_annotated %>% filter(sample %in% samples_IN_MY) %>% filter(uniq_rpkm > 0) %>% filter(grepl("ectin", FUNCTION_AT_DESTINATION_3))) +
#   geom_tile(aes(x = sample, y = feature, fill = uniq_rpkm)) +
#   labs(y = "Cazyme family", x = "sample", fill = "uniq_rpkm") +
#   scale_fill_viridis_c(trans = "log10") +
#   new_scale_fill() +
#   geom_tile(aes(x = sample, y = 0, fill = Host)) +
#   scale_fill_manual(values = host_order_color) +
#   make_theme(setFill = FALSE, leg_pos = "right", guide_nrow = 10, modify_guide = F,
#              x_angle = 90, x_size = 0,
#              x_hj = 10, x_vj = 1,
#              y_size = 10, y_hj = 1)
#              ggsave("results/figures/cazyme_uniq_rpkm_pectin_assoc_familes.pdf")

# ggplot(cazy_df_annotated %>% filter(sample %in% samples_IN_MY) %>% filter(uniq_rpkm > 0)) +
# # ggplot(cazy_df_annotated %>% filter(sample %in% samples_IN_MY) %>% filter(FUNCTION_AT_DESTINATION_3 != "Unknown") %>% filter(!is.na(FUNCTION_AT_DESTINATION_3))) +
#   geom_tile(aes(x = sample, y = short_names, fill = uniq_rpkm)) +
#   labs(y = "Cazyme family", x = "sample", fill = "uniq_rpkm") +
#   scale_fill_viridis_c(trans = "log10") +
#   new_scale_fill() +
#   geom_tile(aes(x = sample, y = 0, fill = Host)) +
#   scale_fill_manual(values = host_order_color) +
#   make_theme(setFill = FALSE, leg_pos = "right", guide_nrow = 10, modify_guide = F,
#              x_angle = 90, x_size = 0,
#              x_hj = 10, x_vj = 1,
#              y_size = 10, y_hj = 1)
#              ggsave("results/figures/cazyme_uniq_rpkm_glyc_annotation.pdf")

# make this with complex heatmap and add glycan annotation

cazyme_uniq_mat_pec <- cazy_df_annotated %>%
                        select(sample, feature, uniq_rpkm) %>%
                        pivot_wider(names_from = sample, values_from = uniq_rpkm) %>%                        
                        column_to_rownames("feature") %>% as.matrix


column_ba <- HeatmapAnnotation(Host = colnames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color
                              ), border = T
                          )
glycan_anno_names <- rownames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown")
row_la <- rowAnnotation(Glycan = glycan_anno_names,
                        col = list(Glycan = make_col_list(glycan_anno_names)
                        ), border = T, show_legend = c(F)
                    )
column_split <- colnames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host) %>% as.factor %>% as.character
row_split <- rownames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown") %>% as.factor %>% as.character
col_fun = colorRamp2(c(0, 0.1, 1, 10, 1000), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[9]))
heatmap_obj = Heatmap(cazyme_uniq_mat_pec,
        # col = c("#e1ebf4", "#081d58"),
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        column_split = column_split,
        row_split = row_split,
        column_dend_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "RPKM (Unique)"),
        bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        left_annotation = row_la,
        # top_annotation = column_ha,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 7),
          cluster_columns = F,
          cluster_rows = F)
draw(heatmap_obj, merge_legend = TRUE)
pdf("results/figures/cazyme_uniq_rpkm.pdf", width = 20, height = 25)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


# make_without unknown
cazy_df_annotated %>% select(feature, FUNCTION_AT_DESTINATION_3) %>% unique %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(FUNCTION_AT_DESTINATION_3 != "Unknown") %>% pull(FUNCTION_AT_DESTINATION_3) %>% unique
known_families <- cazy_df_annotated %>% select(feature, FUNCTION_AT_DESTINATION_3) %>% unique %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(FUNCTION_AT_DESTINATION_3 != "Unknown" & FUNCTION_AT_DESTINATION_3 != "") %>% pull(feature)
cazyme_uniq_mat_pec_no_na <- cazy_df_annotated %>%
                        filter(feature %in% known_families) %>%
                        select(sample, feature, uniq_rpkm) %>%
                        pivot_wider(names_from = sample, values_from = uniq_rpkm) %>%
                        column_to_rownames("feature") %>% as.matrix

column_ba <- HeatmapAnnotation(Host = colnames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color
                              ), border = T
                          )
glycan_anno_names <- rownames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown")
row_la <- rowAnnotation(Glycan = glycan_anno_names,
                        col = list(Glycan = make_col_list(glycan_anno_names)
                        ), border = T, show_legend = c(F)
                    )
column_split <- colnames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host) %>% as.factor %>% as.character
row_split <- rownames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown") %>% as.factor %>% as.character
col_fun = colorRamp2(c(0, 0.1, 1, 10, 1000), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[9]))
heatmap_obj = Heatmap(cazyme_uniq_mat_pec_no_na,
        # col = c("#e1ebf4", "#081d58"),
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        column_split = column_split,
        row_split = row_split,
        column_dend_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "RPKM (Unique)"),
        bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        left_annotation = row_la,
        # top_annotation = column_ha,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 7),
        # show_row_dend = F,
        # show_column_dend = F,
          cluster_columns = F,
          cluster_rows = F)
draw(heatmap_obj, merge_legend = TRUE)
pdf("results/figures/cazyme_uniq_rpkm_no_unknown.pdf", width = 15, height = 25)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()

# cazy_df_annotated <- cazy_df %>%
#                       left_join(glycan_annotations_cleaned, by = c("feature" = "Subfamily")) %>%
#                       mutate(mod_names_pec = ifelse(grepl("pectin", FUNCTION_AT_DESTINATION_3), 1, 0)) %>%
#                       mutate(short_names = abbreviate(FUNCTION_AT_DESTINATION_3, minlength = 50)) %>%
#                       filter(sample %in% samples_IN_MY)

cazyme_uniq_mat_pec <- cazy_df_annotated %>%
                        select(sample, feature, uniq_rpkm) %>%
                        pivot_wider(names_from = sample, values_from = uniq_rpkm) %>%                        
                        column_to_rownames("feature") %>% as.matrix
# remove sparse rows
min <- 15

cazyme_uniq_mat_pec <- cazyme_uniq_mat_pec[apply(cazyme_uniq_mat_pec, 1, function(x) sum(x > 0, na.rm = T)) > min,]

column_ba <- HeatmapAnnotation(Host = colnames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color
                              ), border = T
                          )
glycan_anno_names <- rownames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown")
row_la <- rowAnnotation(Glycan = glycan_anno_names,
                        col = list(Glycan = make_col_list(glycan_anno_names)
                        ), border = T, show_legend = c(F)
                    )
column_split <- colnames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host) %>% as.factor %>% as.character
row_split <- rownames(cazyme_uniq_mat_pec) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown") %>% as.factor %>% as.character
col_fun = colorRamp2(c(0, 0.1, 1, 10, 1000), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[9]))
heatmap_obj = Heatmap(cazyme_uniq_mat_pec,
        # col = c("#e1ebf4", "#081d58"),
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        column_split = column_split,
        row_split = row_split,
        column_dend_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "RPKM (Unique)"),
        bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        left_annotation = row_la,
        # top_annotation = column_ha,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 7),
          cluster_columns = F,
          cluster_rows = F)
draw(heatmap_obj, merge_legend = TRUE)
pdf("results/figures/cazyme_uniq_rpkm_IN_MY.pdf", width = 20, height = 25)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


# make_without unknown
cazy_df_annotated %>% select(feature, FUNCTION_AT_DESTINATION_3) %>% unique %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(FUNCTION_AT_DESTINATION_3 != "Unknown") %>% pull(FUNCTION_AT_DESTINATION_3) %>% unique
known_families <- cazy_df_annotated %>% select(feature, FUNCTION_AT_DESTINATION_3) %>% unique %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(!is.na(FUNCTION_AT_DESTINATION_3)) %>% filter(FUNCTION_AT_DESTINATION_3 != "Unknown" & FUNCTION_AT_DESTINATION_3 != "") %>% pull(feature)
cazyme_uniq_mat_pec_no_na <- cazy_df_annotated %>%
                        filter(feature %in% known_families) %>%
                        select(sample, feature, uniq_rpkm) %>%
                        pivot_wider(names_from = sample, values_from = uniq_rpkm) %>%
                        column_to_rownames("feature") %>% as.matrix

column_ba <- HeatmapAnnotation(Host = colnames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color
                              ), border = T
                          )
glycan_anno_names <- rownames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown")
row_la <- rowAnnotation(Glycan = glycan_anno_names,
                        col = list(Glycan = make_col_list(glycan_anno_names)
                        ), border = T, show_legend = c(F)
                    )
column_split <- colnames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host) %>% as.factor %>% as.character
row_split <- rownames(cazyme_uniq_mat_pec_no_na) %>% as.data.frame() %>% left_join(cazy_df_annotated %>% select(feature, short_names) %>% filter(!is.na(short_names)) %>% unique, by = c("." = "feature")) %>% pull(short_names) %>% replace_na(., "Unknown") %>% as.factor %>% as.character
col_fun = colorRamp2(c(0, 0.1, 1, 10, 1000), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[9]))
heatmap_obj = Heatmap(cazyme_uniq_mat_pec_no_na,
        # col = c("#e1ebf4", "#081d58"),
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        column_split = column_split,
        row_split = row_split,
        column_dend_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "RPKM (Unique)"),
        bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        left_annotation = row_la,
        # top_annotation = column_ha,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 7),
        # show_row_dend = F,
        # show_column_dend = F,
          cluster_columns = F,
          cluster_rows = F)
draw(heatmap_obj, merge_legend = TRUE)
pdf("results/figures/cazyme_uniq_rpkm_no_unknown.pdf", width = 15, height = 25)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()

```

```{r}

# # devtools::install_github("noriakis/ggkegg", dependencies=F)
# library(ggkegg)
# library(ggfx)
# library(igraph)
# library(tidygraph)
# library(dplyr)


# pathway("ko00071") |>
# pathway("ko00020") |>
# pathway("ko01100") |>
#     process_line() |>
#     # highlight_module(module("M00021")) |>
#     # highlight_module(module("M00338")) |>
#     ggraph(x=x, y=y) +
#         geom_node_point(size=1, aes(color=I(fgcolor),
#             filter=fgcolor!="none" & type!="line")) +
#         geom_edge_link0(width=0.1, aes(color=I(fgcolor),
#             filter=type=="line"& fgcolor!="none")) +
#         # with_outer_glow(
#         #     geom_edge_link0(width=1,
#         #         aes(color=I(fgcolor),
#         #             filter=(M00021 | M00338))),
#         #     colour="red", expand=5
#         # ) +
#         # with_outer_glow(
#         #     geom_node_point(size=1.5,
#         #         aes(color=I(fgcolor),
#         #             filter=(M00021 | M00338))),
#         #     colour="red", expand=5
#         # ) +
#         geom_node_text(size=2,
#             aes(x=x, y=y,
#                 label=graphics_name,
#                 filter=name=="path:ko00270"),
#             repel=TRUE, family="sans", bg.colour="white") +
#         theme_void()

# export_kos <- all_sample_genes_df %>%
#             filter(sample %in% c("C3-2")) %>%
#             filter(grepl("Snodgrassella", species)) %>%
#             filter(!is.na(ko)) %>%
#             filter(ko != "") %>%
#             pull(ko)
#             # select(gene, ko)


# pathway("ko01100") |>
#     process_line() |>
#     ggraph(x=x, y=y) +
#         geom_node_point(size=1, aes(color=I(fgcolor),
#             filter=fgcolor!="none" & type!="line")) +
#         geom_edge_link0(width=0.1, aes(color=I(fgcolor),
#             filter=type=="line"& fgcolor!="none")) +
#         geom_node_text(size=12,
#             aes(x=x, y=y,
#                 label=graphics_name,
#                 filter=name=="path:00020"),
#             repel=TRUE, family="sans", bg.colour="white") +
#         theme_void()

# # K04343 - AMR in melli snod

# write.table(export_kos, "results/figures/08-summarize_functions/KOs.csv", row.names = FALSE, quote = FALSE)
```

## 1d. OGs

```{r}
df_og_per_sample <- read.csv("results/figures/08-summarize_functions/ogs_per_sample.tsv", sep = "\t") %>%
                          mutate(host = Vectorize(get_host_from_sample_name)(sample)) %>%
                          mutate(location = Vectorize(get_location_from_sample_name)(sample)) %>%
                            filter(sample %in% samples)
ggplot() +
  geom_boxplot(data = df_og_per_sample, aes(y = ogs, x = factor(host, host_order), 
                fill = factor(host, host_order),
                ),
                outlier.shape = NA
              ) +
  geom_point(data = df_og_per_sample, aes(y = ogs, x = factor(host, host_order), group = host,
                 shape = location, fill = factor(host)), 
              size = 3,
              position = position_jitterdodge(0.75)) +
  scale_fill_manual(values = host_order_color_dark) +
  labs(x = "Host", y = "Number of OGs", fill = "Host") +
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
  # add stats
  geom_signif(data = df_og_per_sample,
              aes(x = host, y = ogs),
              test = "wilcox.test",
              step_increase = 0.05,
              tip_length = 0,
              comparisons = list(c("Apis florea", "Apis andreniformis"), 
                                #  c("Apis mellifera", "Apis florea"),
                                #  c("Apis mellifera", "Apis andreniformis"),
                                #  c("Apis cerana", "Apis dorsata"),
                                #  c("Apis cerana", "Apis florea"),
                                #  c("Apis dorsata", "Apis florea"),
                                 c("Apis mellifera", "Apis cerana"),
                                 c("Apis dorsata", "Apis andreniformis"),
                                 c("Apis dorsata", "Apis cerana"),
                                 c("Apis mellifera", "Apis dorsata")),
              map_signif_level = TRUE, textsize = 6) +
  make_theme(setFill = FALSE, leg_pos = "bottom", 
             x_angle = 90, x_size = 0,
             x_hj = 10, x_vj = 1,
             y_size = 10, y_hj = 1)
ggsave("results/figures/Number_of_ogs_per_sample.pdf")

df_og_per_sample_split <- df_og_per_sample  %>%
                      mutate(accessory = ogs-core) %>%
                      pivot_longer(c("accessory", "core"), names_to = "type", values_to = "number")

ggplot(df_og_per_sample_split) +
  geom_boxplot(aes(y = number, x = factor(host, host_order), 
                fill = type
                ),
                outlier.shape = NA
              ) +
  # geom_point(aes(y = number, x = factor(host, host_order), group = host,
  #                shape = location, fill = type),
  #             size = 3,
  #             position = position_jitterdodge(0.75)) +
  scale_fill_manual(values = c("accessory" = "grey", "core" = "black")) +
  labs(x = "Host", y = "Number of OGs", fill = "Type") +
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))
  ggsave("results/figures/Number_of_ogs_per_sample_core.pdf")

df_og_cum_data <- read.csv("results/figures/08-summarize_functions/og_cum_curve.tsv", sep = "\t")
ggplot() +
  geom_smooth(data = df_og_cum_data %>% filter(host == "Apis mellifera"), aes(x = size, y = ogs, color = host), method = "loess", se = T) +
  geom_point(data = df_og_cum_data %>% filter(host == "Apis mellifera"), aes(x = size, y = ogs, color = host)) +
  geom_smooth(data = df_og_cum_data %>% filter(host == "Apis cerana"), aes(x = size, y = ogs, color = host), method = "loess", se = T) +
  geom_point(data = df_og_cum_data %>% filter(host == "Apis cerana"), aes(x = size, y = ogs, color = host)) +
  geom_smooth(data = df_og_cum_data %>% filter(host == "Apis dorsata"), aes(x = size, y = ogs, color = host), method = "loess", se = T) +
  geom_point(data = df_og_cum_data %>% filter(host == "Apis dorsata"), aes(x = size, y = ogs, color = host)) +
  geom_smooth(data = df_og_cum_data %>% filter(host == "Apis florea"), aes(x = size, y = ogs, color = host), method = "loess", se = T) +
  geom_point(data = df_og_cum_data %>% filter(host == "Apis florea"), aes(x = size, y = ogs, color = host)) +
  geom_smooth(data = df_og_cum_data %>% filter(host == "Apis andreniformis"), aes(x = size, y = ogs, color = host), method = "loess", se = T) +
  geom_point(data = df_og_cum_data %>% filter(host == "Apis andreniformis"), aes(x = size, y = ogs, color = host)) +
  scale_color_manual(values = host_order_color_dark) +
  labs(x = "# Bees", y = "# OGs", color = "Host") +
  # scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
  make_theme(setFill = F, setCol = F, guide_nrow = 1)
ggsave("results/figures/OG_cum_curve.pdf")
```

```{r}
df_og_sharing <- read.csv("results/figures/08-summarize_functions/ogs_per_host.tsv", sep = "\t") %>%
                  mutate(number = as.numeric(number)) %>%
                  filter(number > 0) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  # group_by(genus, kos_C) %>%
                  # mutate(num_hosts = n_distinct(host)) %>%
                  # filter(num_hosts > 1) %>%
                  # ungroup() %>%
                  filter(host != "shared") %>%
                  group_by(genus, kos_C) %>%
                  filter(!genus %in% c("Klebsiella", "Enterobacter")) %>%
                  mutate(n = n()) %>%
                  filter(kos_C != "Unknown") %>%
                  filter(number > 10)
ggplot(df_og_sharing) +
  geom_bar(aes(y = kos_C, x = number, fill = host), stat = "identity") +
  facet_wrap(~genus, scales = "free") +
  scale_fill_manual(values = host_order_color_dark, na.value = "grey") +
  labs(x = "Host", y = "Number of OGs", fill = "Host") +
  make_theme(setFill = F, y_size = 6,
             y_hj = 1, y_vj = 0.5,
             x_hj = 1, x_vj = 0.5,
             x_angle = 90, x_size = 5)
```

# 2a. KOs

```{r}
# get feature matrix and plot as heatmap (KO, cazy, og?)
  # make sample names the row names and make it a numeric matrix
ko_df <- read.csv("results/figures/08-summarize_functions/ko_matrix.csv") %>%
          filter(sample %in% samples) %>%
          column_to_rownames("sample") %>%
            as.matrix()
# dist_mat <- vegdist(ko_df, method = "bray")
# pcoa_plot(dist_mat, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)

# presence-absence matrix of 1 and 0
ko_df_pa <- ko_df
ko_df_pa[ko_df_pa > 0] <- 1

# dist_mat <- beta.pair(ko_df_pa, "sorensen")
# pcoa_plot(dist_mat$beta.sim, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)
# pcoa_plot(dist_mat$beta.sor, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)

dist_mat <- beta.pair(ko_df_pa, "jaccard")
res <- adonis2(dist_mat$beta.jac ~ Species, 
          rownames(dist_mat$beta.jac %>% as.matrix) %>% as.data.frame() %>%
            mutate(Species = Vectorize(get_host_from_sample_name)(.)))
pdf("results/figures/KO_jaccard_pcoa.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat$beta.jac, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark) +
  annotate("text", x = -0.15,
                      y = 0.2,
                      size = 8,
          label = paste0("omega2 = ", round(adonis_OmegaSq(res)$parOmegaSq[[1]], 3),
                         ", \np = ", adonis_OmegaSq(res)$`Pr(>F`[[1]])
          )
dev.off()
res <- adonis2(dist_mat$beta.jtu ~ Species, 
          rownames(dist_mat$beta.jtu %>% as.matrix) %>% as.data.frame() %>%
            mutate(Species = Vectorize(get_host_from_sample_name)(.)))
pdf("results/figures/KO_jaccard_turnover.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat$beta.jtu, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark) +
  annotate("text", x = -0.1,
                      y = 0.13,
                      size = 8,
          label = paste0("omega2 = ", round(adonis_OmegaSq(res)$parOmegaSq[[1]], 3),
                         ", \np = ", adonis_OmegaSq(res)$`Pr(>F`[[1]])
          )
dev.off()

dist_mat <- beta.pair(ko_df_pa, "sorensen")
res <- adonis2(dist_mat$beta.sor ~ Species, 
          rownames(dist_mat$beta.sor %>% as.matrix) %>% as.data.frame() %>%
            mutate(Species = Vectorize(get_host_from_sample_name)(.)))
res_omega <- adonis_OmegaSq(res)
pdf("results/figures/KO_sorensen_pcoa.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat$beta.sor, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark) +
  annotate("text", x = -0.15,
                      y = 0.2,
                      size = 8,
          label = paste0("omega2 = ", round(adonis_OmegaSq(res)$parOmegaSq[[1]], 3),
                         ", \np = ", adonis_OmegaSq(res)$`Pr(>F`[[1]])
          )
dev.off()
res <- adonis2(dist_mat$beta.sim ~ Species, 
          rownames(dist_mat$beta.sim %>% as.matrix) %>% as.data.frame() %>%
            mutate(Species = Vectorize(get_host_from_sample_name)(.)))
res_omega <- adonis_OmegaSq(res)  
pdf("results/figures/KO_sorensen_turnover_pcoa.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat$beta.sim, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark) +
  annotate("text", x = -0.1,
                      y = 0.13,
                      size = 8,
          label = paste0("omega2 = ", round(adonis_OmegaSq(res)$parOmegaSq[[1]], 3),
                         ", \np = ", adonis_OmegaSq(res)$`Pr(>F`[[1]])
          )
dev.off()


ko_matrix_norm <- decostand(ko_df, method = "normalize")
dist_mat <- vegdist(ko_matrix_norm, method = "robust.aitchison")
ado_res <- adonis2(dist_mat ~ Species, 
          rownames(dist_mat %>% as.matrix) %>% as.data.frame() %>%
            mutate(Species = Vectorize(get_host_from_sample_name)(.)))
pdf("results/figures/KO_robust_aitchison_pcoa.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark) +
  annotate("text", x = -25,
                      y = -10,
                      size = 8,
          label = paste0("omega2 = ", round(adonis_OmegaSq(ado_res)$parOmegaSq[[1]], 3),
                         ", \np = ", adonis_OmegaSq(ado_res)$`Pr(>F`[[1]])
          )
dev.off()



# # load host info from scripts/visualization/09-plot_community_composition.Rmd
# sample_bac_div <- beta.pair(ko_df_pa, "sorensen")$beta.sim
# # heatmap(sample_bac_div %>% as.matrix)
# sample_host_div <- sample_host_div %>% 
#   filter(row.names(.) %in% as.vector(as.matrix(sample_bac_div) %>% row.names)) %>%
#   select(as.vector(as.matrix(sample_bac_div) %>% row.names))
# # heatmap(sample_host_div %>% as.matrix)

# sample_host_div_df <- sample_host_div %>%
#                         as.matrix %>%
#                         as.data.frame %>%
#                         rownames_to_column("Sample") %>%
#                         filter(Sample %in% samples_IN_MY) %>%
#                         pivot_longer(!Sample, names_to = "Compared", values_to = "dist_h")
# sample_bac_div_df <- sample_bac_div %>%
#                         as.matrix %>%
#                         as.data.frame %>%
#                         rownames_to_column("Sample") %>%
#                         filter(Sample %in% samples_IN_MY) %>%
#                           pivot_longer(!Sample, names_to = "Compared", values_to = "dist_b") %>%
#                           mutate(Host_s = Vectorize(get_host_from_sample_name)(Sample)) %>%
#                           mutate(Host_c = Vectorize(get_host_from_sample_name)(Compared)) %>%
#                           mutate(Type_b = ifelse(Host_s == Host_c, "Within", "Between")) %>%
#                           mutate(Hosts_compared = paste0(Host_s, "_", Host_c))
# dissimilarities_df_median <- left_join(sample_bac_div_df, sample_host_div_df) %>%
#                         group_by(Hosts_compared) %>%
#                         summarise(dist_h, median_b = median(dist_b)) %>%
#                         unique
# dissimilarities_df_mean <- left_join(sample_bac_div_df, sample_host_div_df) %>%
#                         group_by(Hosts_compared) %>%
#                         summarise(dist_h, mean_b = median(dist_b)) %>%
#                         unique
# dissimilarities_df <- left_join(sample_bac_div_df, sample_host_div_df) %>%
#                         group_by(Hosts_compared)
#                         # filter(Sample %in% c("M1.2", "C1.2", "D1.2", "F1.2", "A1.2"))

# anno_side_b <- rowAnnotation(Host = sample_bac_div %>% as.matrix %>% rownames %>% as.data.frame() %>%
#                                 mutate(host = Vectorize(get_host_from_sample_name)(.)) %>% pull(host),
#                                col = list(Host = host_order_color_dark),
#                                border = TRUE
# )
# anno_side_h <- rowAnnotation(Host = sample_host_div %>% as.matrix %>% rownames %>% as.data.frame() %>%
#                                 mutate(host = Vectorize(get_host_from_sample_name)(.)) %>% pull(host),
#                                col = list(Host = host_order_color_dark),
#                                border = TRUE
# )
# anno_top_b <- HeatmapAnnotation(Host = sample_bac_div %>% as.matrix %>% colnames %>% as.data.frame() %>%
#                                 mutate(host = Vectorize(get_host_from_sample_name)(.)) %>% pull(host),
#                                col = list(Host = host_order_color_dark),
#                                border = TRUE
# )
# anno_top_h <- HeatmapAnnotation(Host = sample_host_div %>% as.matrix %>% colnames %>% as.data.frame() %>%
#                                 mutate(host = Vectorize(get_host_from_sample_name)(.)) %>% pull(host),
#                                col = list(Host = host_order_color_dark),
#                                border = TRUE
# )
# pdf("results/figures/KO_bac_divergence_sorenson_turnover_dend.pdf")
# Heatmap(sample_bac_div %>% as.matrix,
#         rect_gp = gpar(type = "none"),
#         top_annotation = anno_top_b,
#         right_annotation = anno_side_b
#       )
# dev.off()
# pdf("results/figures/KO_host_divergence_sorenson_turnover_dend.pdf")
# Heatmap(sample_host_div %>% as.matrix,
#         rect_gp = gpar(type = "none"),
#         top_annotation = anno_top_h,
#         right_annotation = anno_side_h
#       )
# dev.off()

# mantel(sample_host_div, sample_bac_div, method="spearman", permutations=1000)

# cor_res <- cor.test(dissimilarities_df$dist_h, dissimilarities_df$dist_b, method = "spearman")

# ggplot() +
#   geom_jitter(data=dissimilarities_df,
#              aes(x = dist_h,
#                  y = dist_b,
#                  group = Hosts_compared,
#                  fill = Host_s,
#                  color = Host_c,
#                 #  shape = Type_b
#              ), width=0.0001, size = 4, shape = 21, stroke = 2
#             ) +
#   # geom_jitter(data=dissimilarities_df_mean,
#   #            aes(x = dist_h,
#   #                y = mean_b
#   #            ), size = 3, width = 0
#   #           ) +
#   geom_jitter(data=dissimilarities_df_median,
#              aes(x = dist_h,
#                  y = median_b
#              ), size = 5, width = 0, color = "black", shape = 23, fill = "red", stroke = 2
#             ) +
#   # geom_smooth(data=dissimilarities_df_median,
#   #            aes(x = dist_h,
#   #                y = median_b
#   #            ), method = "lm", se = T, color = "black"
#   #           ) +
#     labs(x = "Host divergence (based on 16S and 12S)", y = "Microbiome dissimilarity (sorensen - turnover)", color = "Host species") +
#     scale_fill_manual(values = host_order_color_dark) +
#     scale_color_manual(values = host_order_color_dark) +
#     annotate("text", x = 0.02,
#                       y = 0.4,
#                       size = 5,
#           label = paste0("Spearman's rho = ", round(cor_res$estimate, 3),
#                          ", \np = ", cor_res$p.value)
#           ) +
#     make_theme(setCol = F, setFill = F, guide_nrow = 2,
#                leg_size = 12,
#                axis_x_title = 20, axis_y_title = 20
#     )
#   ggsave("results/figures/KO_host_dissimilarities_sorenson_turnover.pdf")

# mantel(as.dist(sample_host_div), as.dist(sample_bac_div), method = "pearson", permutations = 9999)

```

# 2b. CAZy

```{r}
# get feature matrix and plot as heatmap (KO, cazy, og?)
  # make sample names the row names and make it a numeric matrix
# unlist(lapply(colnames(cazy_df), function(x) paste(c(strsplit(x, "")[[1]])[1:2], collapse = ""))) %>% unique
cazy_df <- read.csv("results/figures/08-summarize_functions/cazyme_matrix.csv") %>%
          filter(sample %in% samples) %>%
          select(-starts_with("AA"), -starts_with("GT")) %>%
          column_to_rownames("sample") %>%          
            as.matrix()
# dist_mat <- vegdist(cazy_df, method = "bray")
# pcoa_plot(dist_mat, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)

# presence-absence matrix of 1 and 0
# remove columns with no data
cazy_df <- cazy_df[, colSums(cazy_df) > 0]
heatmap(cazy_df)
cazy_df_pa <- cazy_df
cazy_df_pa[cazy_df_pa > 0] <- 1
cazy_df_pa <- cazy_df_pa[apply(cazy_df_pa, 1, function(x) sum(x) > 0),]

dist_mat <- beta.pair(cazy_df_pa, "jaccard")
pdf("results/figures/CAZy_jaccard_pcoa.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat$beta.jac, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)
dev.off()
adonis2(dist_mat$beta.jac ~ Species, 
          rownames(dist_mat$beta.jac %>% as.matrix) %>% as.data.frame() %>%
            mutate(Species = Vectorize(get_host_from_sample_name)(.)))

dist_mat <- beta.pair(cazy_df_pa, "sorensen")
pdf("results/figures/CAZy_sorensen_pcoa.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat$beta.sor, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)
# pcoa_plot(dist_mat$beta.sim, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)
dev.off()
# pcoa_plot(dist_mat$beta.jne, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)
# pcoa_plot(dist_mat$beta.sim, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark)

cazy_df_norm <- decostand(cazy_df, method = "normalize")
dist_mat <- vegdist(cazy_df_norm, method = "robust.aitchison")
ado_res <- adonis2(dist_mat ~ Species, 
          rownames(dist_mat %>% as.matrix) %>% as.data.frame() %>%
            mutate(Species = Vectorize(get_host_from_sample_name)(.)))
pdf("results/figures/CAZy_robust_aitchison_pcoa.pdf", , width = 10, height = 15)
pcoa_plot(dist_mat, df_meta_complete, "Species", color_add = T, color_list = host_order_color_dark) +
  annotate("text", x = -5,
                      y = 4,
                      size = 8,
          label = paste0("omega2 = ", round(adonis_OmegaSq(ado_res)$parOmegaSq[[1]], 3),
                         ", \np = ", adonis_OmegaSq(ado_res)$`Pr(>F`[[1]])
          )
dev.off()

```

# Finding significant associations using Maaslin2 and lsmeans

```{r}
# set R path to:
# /work/FAC/FBM/DMF/pengel/spirit/aprasad/miniconda3/envs/Apis_project_r_env/bin/R
  # used to use,
  # /work/FAC/FBM/DMF/pengel/spirit/aprasad/miniconda3/envs/masslin2/bin/R
  # /work/FAC/FBM/DMF/pengel/spirit/aprasad/miniconda3/envs/rmd-env/bin/R
# Defaults:
# Normalization : TSS
# Transform : Log
# Analysis method : LM
# Min prev : 0.1
# Min abundance : 0
# Correction : BH
# Apply z-score so continuous metadata are 
#         on the same scale [ Default: TRUE ]

# # mellifera as 'reference' - according to the order of levels in factor
# ko_matrix <- read.table("results/figures/08-summarize_functions/ko_matrix.csv", sep = ",", header = T, row.names = 1)
# input_data <- ko_matrix[!rownames(ko_matrix) %in% samples_to_exclude,]
# meta_table <- data.frame(Host = Vectorize(get_host_from_sample_name)(rownames(input_data)))
# meta_table$Host <- factor(meta_table$Host, levels = c("Apis mellifera", "Apis cerana",
#                                                       "Apis dorsata", "Apis florea", "Apis andreniformis"))
# input_metadata <- meta_table
# fit_data = Maaslin2(
#     input_data = input_data, 
#     input_metadata = input_metadata, 
#     output = "results/figures/08-summarize_functions/maaslin2_results_ko", 
#     analysis_method = "LM", # consider
#     transform = "LOG", # consider
#     normalization = "TSS", # consider
#     fixed_effects = c("Host"),
#     min_abundance = 0,
#     save_models = T,
#     plot_scatter = F,
#     plot_heatmap = F)

# # get pairwise comparisons using lsmeans
# # and write to file

# # get the models
# df_pairwise <- data.frame()
# models_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_ko/fits/models.rds")
# for (i in 1:length(models_data)) {
#   model <- models_data[[i]]
#   conts <- summary(lsmeans(model, pairwise ~ Host, adjust = "FDR"))
#   df_to_add <- data.frame(
#                     feature = names(models_data)[[i]],
#                     contrast = conts$contrasts$contrast,
#                     p.value = conts$contrasts$p.value,
#                     estimate = conts$contrasts$estimate
#                   )
#   df_pairwise <- rbind(df_pairwise, df_to_add)
# }
# write.csv(df_pairwise, "results/figures/08-summarize_functions/maaslin2_results_ko_pairwise.csv", row.names = F, quote = F)


# cazyme_matrix <- read.table("results/figures/08-summarize_functions/cazyme_matrix.csv", sep = ",", header = T, row.names = 1)
# input_data <- cazyme_matrix[!rownames(cazyme_matrix) %in% samples_to_exclude,]
# meta_table <- data.frame(Host = Vectorize(get_host_from_sample_name)(rownames(input_data)))
# meta_table$Host <- factor(meta_table$Host, levels = c("Apis mellifera", "Apis cerana",
#                                                       "Apis dorsata", "Apis florea", "Apis andreniformis"))
# input_metadata <- meta_table
# fit_data = Maaslin2(
#     input_data = input_data, 
#     input_metadata = input_metadata, 
#     output = "results/figures/08-summarize_functions/maaslin2_results_cazyme", 
#     analysis_method = "LM", # consider
#     transform = "LOG", # consider
#     normalization = "TSS", # consider
#     fixed_effects = c("Host"),
#     min_abundance = 0,
#     save_models = T,
#     plot_scatter = F,
#     plot_heatmap = F)

# # get pairwise comparisons using lsmeans
# # and write to file
# df_pairwise <- data.frame()
# models_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_cazyme/fits/models.rds")
# for (i in 1:length(models_data)) {
#   print(paste0("Processing ", i, " of ", length(models_data)))
#   model <- models_data[[i]]
#   conts <- summary(lsmeans(model, pairwise ~ Host, adjust = "FDR"))
#   df_to_add <- data.frame(
#                     feature = names(models_data)[[i]],
#                     contrast = conts$contrasts$contrast,
#                     p.value = conts$contrasts$p.value,
#                     estimate = conts$contrasts$estimate
#                   )
#   df_pairwise <- rbind(df_pairwise, df_to_add)
# }
# write.csv(df_pairwise, "results/figures/08-summarize_functions/maaslin2_results_cazyme_pairwise.csv", row.names = F, quote = F)

og_matrix <- read.table("results/figures/08-summarize_functions/og_matrix.csv", sep = ",", header = T, row.names = 1)
input_data <- og_matrix[!rownames(og_matrix) %in% samples_to_exclude,]
meta_table <- data.frame(Host = Vectorize(get_host_from_sample_name)(rownames(input_data)))
meta_table$Host <- factor(meta_table$Host, levels = c("Apis mellifera", "Apis cerana",
                                                      "Apis dorsata", "Apis florea", "Apis andreniformis"))
input_metadata <- meta_table
fit_data = Maaslin2(
    input_data = input_data, 
    input_metadata = input_metadata, 
    output = "results/figures/08-summarize_functions/maaslin2_results_og", 
    analysis_method = "LM", # consider
    transform = "LOG", # consider
    normalization = "TSS", # consider
    fixed_effects = c("Host"),
    min_abundance = 0,
    save_models = T,
    plot_scatter = F,
    plot_heatmap = F)

# get pairwise comparisons using lsmeans
# and write to file
df_pairwise <- data.frame()
models_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_og/fits/models.rds")
for (i in 1:length(models_data)) {
  model <- models_data[[i]]
  conts <- summary(lsmeans(model, pairwise ~ Host, adjust = "FDR"))
  df_to_add <- data.frame(
                    feature = names(models_data)[[i]],
                    contrast = conts$contrasts$contrast,
                    p.value = conts$contrasts$p.value,
                    estimate = conts$contrasts$estimate
                  )
  df_pairwise <- rbind(df_pairwise, df_to_add)
}
write.csv(df_pairwise, "results/figures/08-summarize_functions/maaslin2_results_og_pairwise.csv", row.names = F, quote = F)


# Just some sanity checks
# "g__Apibacter--OG0000704" # is listed as significant in a cerana - florea comparison but florea has none of this.
# Try to validate (to some extent, the model used)
residuals_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_ko/fits/residuals.rds")
models_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_ko/fits/models.rds")
fitted_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_ko/fits/fitted.rds")

models[[1]]
plot(residuals_data[24,] ~ fitted_data[24,], col = alpha("blue", 0.4))
plot(residuals_data[14,] ~ fitted_data[14,], col = alpha("blue", 0.4))
plot(residuals_data[50,] ~ fitted_data[50,], col = alpha("blue", 0.4))
hist(residuals_data[4,])
slopes <- c()
for (i in 1:dim(fitted_data)[1]) {
  slopes <- c(slopes, lm(fitted_data[i,] ~ residuals_data[i,])$coef[[2]])
}
# so long as the slope is close to zero, we can assume that the model is good
# because the residuals are not dependent on the fitted values which
# would be the case if the distribution of the residuals was not normal(ish)
# and the model was not good....
hist(slopes,100)

residuals_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_cazyme/fits/residuals.rds")
models_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_cazyme/fits/models.rds")
fitted_data <- readRDS("results/figures/08-summarize_functions/maaslin2_results_cazyme/fits/fitted.rds")

models[[1]]
plot(residuals_data[24,] ~ fitted_data[24,], col = alpha("blue", 0.4))
plot(residuals_data[14,] ~ fitted_data[14,], col = alpha("blue", 0.4))
plot(residuals_data[50,] ~ fitted_data[50,], col = alpha("blue", 0.4))
hist(residuals_data[4,])
hist(residuals_data[30,])
slopes <- c()
for (i in 1:dim(fitted_data)[1]) {
  slopes <- c(slopes, lm(fitted_data[i,] ~ residuals_data[i,])$coef[[2]])
}
# in sime cases the slope is a little higher or lower than zero
# and this means that the model is not perfect but the p-value
# is still significant so we can still use the results to inform
# our understanding of the data
hist(slopes,100)

```

# Using random forest to find important features as done in microbiomeSeq

```{r}
# conda activate masslin2
# using code from a function in microbiomeSeq
# https://github.com/umerijaz/microbiomeSeq/blob/master/R/kruskal_abundance.R
ko_matrix <- read.table("results/figures/08-summarize_functions/ko_matrix.csv", sep = ",", header = T, row.names = 1)
pvalue.threshold <- 0.05
abund_table <- ko_matrix[!rownames(ko_matrix) %in% samples_to_exclude,]
# filter and normalize (log-relative) data
abund_table <- abund_table[rowSums(abund_table) > 0,]
abund_table <- abund_table / rowSums(abund_table)
abund_table <- log(abund_table + 1)
meta_table <- data.frame(Groups = Vectorize(get_host_from_sample_name)(rownames(abund_table)))
meta_table$Groups <- as.factor(meta_table$Groups)

kruskal.wallis.table <- data.frame()
data <- as.data.frame(abund_table)
for (i in 1:dim(data)[2]){
    ks.test <- kruskal.test(data[,i], g=meta_table$Groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id=names(data)[i],p.value=ks.test$p.value))
}
kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing=FALSE), ]
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
rownames(kruskal.wallis.table) <- kruskal.wallis.table$id

#==================significant feature selection =================================#
last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
selected <- 1:last.significant.element
sig_res <-kruskal.wallis.table$id[selected]

#==random forest classifier ==#
subset.data<-data.frame(data[,as.character(kruskal.wallis.table[rownames(kruskal.wallis.table),"id"])])
kruskal.wallis.table$id <- colnames(subset.data) #enforce that ids and colnames of subset data remain the same for easy indexing later on
subset.data <- subset.data[,sig_res]
rf_res <- randomforest_res(subset.data, meta_table$Groups)
df_accuracy <- rf_res$importance

df <- NULL
for(i in df_accuracy$Sample){
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," p.adj = ",sprintf("%.10g",kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"]),sep=""),dim(data)[1]))
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    df <- na.omit(df)
}

out <- list("SignfeaturesTable"=kruskal.wallis.table, "plotdata"=df, "importance"=df_accuracy)
write.csv(out$SignfeaturesTable, "results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_table.csv")
write.csv(out$importance, "results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance.csv")


cazyme_DRAM_matrix <- read.table("results/figures/08-summarize_functions/cazyme-DRAM_matrix.csv", sep = ",", header = T, row.names = 1)
pvalue.threshold <- 0.05
abund_table <- cazyme_DRAM_matrix[!rownames(cazyme_DRAM_matrix) %in% samples_to_exclude,]
# filter and normalize (log-relative) data
abund_table <- abund_table[rowSums(abund_table) > 0,]
abund_table <- abund_table / rowSums(abund_table)
abund_table <- log(abund_table + 1)
meta_table <- data.frame(Groups = Vectorize(get_host_from_sample_name)(rownames(abund_table)))
meta_table$Groups <- as.factor(meta_table$Groups)

kruskal.wallis.table <- data.frame()
data <- as.data.frame(abund_table)
for (i in 1:dim(data)[2]){
    ks.test <- kruskal.test(data[,i], g=meta_table$Groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id=names(data)[i],p.value=ks.test$p.value))
}
kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing=FALSE), ]
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
rownames(kruskal.wallis.table) <- kruskal.wallis.table$id

#==================significant feature selection =================================#
last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
selected <- 1:last.significant.element
sig_res <-kruskal.wallis.table$id[selected]

#==random forest classifier ==#
subset.data<-data.frame(data[,as.character(kruskal.wallis.table[rownames(kruskal.wallis.table),"id"])])
kruskal.wallis.table$id <- colnames(subset.data) #enforce that ids and colnames of subset data remain the same for easy indexing later on
subset.data <- subset.data[,sig_res]
rf_res <- randomforest_res(subset.data, meta_table$Groups)
df_accuracy <- rf_res$importance

df <- NULL
for(i in df_accuracy$Sample){
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," p.adj = ",sprintf("%.10g",kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"]),sep=""),dim(data)[1]))
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    df <- na.omit(df)
}

out <- list("SignfeaturesTable"=kruskal.wallis.table, "plotdata"=df, "importance"=df_accuracy)
write.csv(out$SignfeaturesTable, "results/figures/08-summarize_functions/cazyme_DRAM_compared_kruskal_wallis_table.csv")
write.csv(out$importance, "results/figures/08-summarize_functions/cazyme_DRAM_compared_kruskal_wallis_importance.csv")

cazyme_matrix <- read.table("results/figures/08-summarize_functions/cazyme_matrix.csv", sep = ",", header = T, row.names = 1) %>%
                    select(-starts_with("AA"), -starts_with("GT"))
pvalue.threshold <- 0.05
abund_table <- cazyme_matrix[!rownames(cazyme_matrix) %in% samples_to_exclude,]
# filter and normalize (log-relative) data
abund_table <- abund_table[rowSums(abund_table) > 0,]
abund_table <- abund_table / rowSums(abund_table)
abund_table <- log(abund_table + 1)
meta_table <- data.frame(Groups = Vectorize(get_host_from_sample_name)(rownames(abund_table)))
meta_table$Groups <- as.factor(meta_table$Groups)

kruskal.wallis.table <- data.frame()
data <- as.data.frame(abund_table)
for (i in 1:dim(data)[2]){
    ks.test <- kruskal.test(data[,i], g=meta_table$Groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id=names(data)[i],p.value=ks.test$p.value))
}
kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing=FALSE), ]
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
rownames(kruskal.wallis.table) <- kruskal.wallis.table$id

#==================significant feature selection =================================#
last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
selected <- 1:last.significant.element
sig_res <-kruskal.wallis.table$id[selected]

#==random forest classifier ==#
subset.data<-data.frame(data[,as.character(kruskal.wallis.table[rownames(kruskal.wallis.table),"id"])])
kruskal.wallis.table$id <- colnames(subset.data) #enforce that ids and colnames of subset data remain the same for easy indexing later on
subset.data <- subset.data[,sig_res]
rf_res <- randomforest_res(subset.data, meta_table$Groups)
df_accuracy <- rf_res$importance

df <- NULL
for(i in df_accuracy$Sample){
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," p.adj = ",sprintf("%.10g",kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"]),sep=""),dim(data)[1]))
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    df <- na.omit(df)
}

out <- list("SignfeaturesTable"=kruskal.wallis.table, "plotdata"=df, "importance"=df_accuracy)
write.csv(out$SignfeaturesTable, "results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_table.csv")
write.csv(out$importance, "results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_importance.csv")

og_matrix <- read.table("results/figures/08-summarize_functions/og_matrix.csv", sep = ",", header = T, row.names = 1)
pvalue.threshold <- 0.05
abund_table <- og_matrix[!rownames(og_matrix) %in% samples_to_exclude,]
# filter and normalize (log-relative) data
abund_table <- abund_table[rowSums(abund_table) > 0,]
abund_table <- abund_table / rowSums(abund_table)
abund_table <- log(abund_table + 1)
meta_table <- data.frame(Groups = Vectorize(get_host_from_sample_name)(rownames(abund_table)))
meta_table$Groups <- as.factor(meta_table$Groups)

kruskal.wallis.table <- data.frame()
data <- as.data.frame(abund_table)
for (i in 1:dim(data)[2]){
    ks.test <- kruskal.test(data[,i], g=meta_table$Groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id=names(data)[i],p.value=ks.test$p.value))
}
kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing=FALSE), ]
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
rownames(kruskal.wallis.table) <- kruskal.wallis.table$id

#==================significant feature selection =================================#
last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
selected <- 1:last.significant.element
sig_res <-kruskal.wallis.table$id[selected]

#==random forest classifier ==#
subset.data<-data.frame(data[,as.character(kruskal.wallis.table[rownames(kruskal.wallis.table),"id"])])
kruskal.wallis.table$id <- colnames(subset.data) #enforce that ids and colnames of subset data remain the same for easy indexing later on
subset.data <- subset.data[,sig_res]
rf_res <- randomforest_res(subset.data, meta_table$Groups)
df_accuracy <- rf_res$importance

df <- NULL
for(i in df_accuracy$Sample){
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," p.adj = ",sprintf("%.10g",kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"]),sep=""),dim(data)[1]))
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    df <- na.omit(df)
}

out <- list("SignfeaturesTable"=kruskal.wallis.table, "plotdata"=df, "importance"=df_accuracy)
write.csv(out$SignfeaturesTable, "results/figures/08-summarize_functions/og_compared_kruskal_wallis_table.csv")
write.csv(out$importance, "results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance.csv")

```

```{r}
# plot masslin2 results
"results/figures/08-summarize_functions/maaslin2_results/significant_results_annotated.tsv"

```
```{r}
plot_signif <- function(df=NULL,top.taxa=20,...){

  #==plot the significant features information and random classifier results
  p<-NULL
  df <- na.omit(df)
  #pick the top significant taxa, this is to avoid overwhelming clutter.
  df <- df[which(as.numeric(df$Rank)%in%(1:top.taxa)),]
  max.rank <- max(as.numeric(df$Rank))
  df$bar_height <- max.rank-as.numeric(df$Rank)+1
  if(max.rank >=20){
  df$bar_height <- 1 + max.rank-round(as.numeric(df$Rank)/(top.taxa/10))
  }
  df$rank_label <-NULL
  for (i in 1:dim(df)[1]){
    rank <- as.numeric(df$Rank[i])
    df$rank_label[i] <- paste(paste(rep("\u25ac", df$bar_height[i]),collapse=""),rank) #"\u25aa" "\u25ae" "\u25ac"
  }

  if(!is.null(df)){
    p<-ggplot(df,aes(Groups,Value,colour=Groups))
    p<-p+geom_boxplot(outlier.size=NA)+geom_jitter(position = position_jitter(height = 0, width=0))+theme_bw()
    p<-p+ facet_grid( ~rank_label+Taxa, scales="free_x")
    p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text = element_text(size = 10, colour = "black", angle = 90,vjust=0)) #vjust=0 aligns the labels to the bottom in this case
    p<-p+theme(strip.background =element_rect(fill="white"))+theme(plot.margin = unit(c(1, 1, 0, 1), "lines"))#set background colour for facet lbels
  }

  return(p)
}
library(ggplot2)
plot_signif(out$plotdata, top.taxa = 10)
plot_MDA <- function(df_accuracy, top.taxa=20){
  mda_plot <-NULL
  if(!is.null(df_accuracy)){
    df_accuracy <- df_accuracy[which(df_accuracy$rank%in%c(1:top.taxa)),]
    mda_plot <- ggplot(data = df_accuracy,aes(x=Sample,y=Value)) + theme_bw()
    mda_plot <- mda_plot+geom_bar(stat = "identity",fill="darkblue",width = 0.5)
    mda_plot <- mda_plot + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    mda_plot <- mda_plot + xlab("Taxa description") + ylab("Mean Decrease in Accuracy")
  }
  return(mda_plot)
}
plot_MDA(df_accuracy, top.taxa = 10)
plot_corrections <- function(corrections_table, pvalue.Cutoff){

  if(!is.null(corrections_table)){
    plot(corrections_table$p.value, corrections_table$E.value,main='Multitesting corrections',
         xlab='Nominal p-value',ylab='Multitesting-corrected statistics',log='xy',col='blue',panel.first=grid(col='#BBBBBB',lty='solid'))
    lines(corrections_table$p.value,corrections_table$FWER,pch=20,col='darkgreen', type='p')
    lines(corrections_table$p.value,corrections_table$q.value,pch='+',col='darkred', type='p')
    abline(h=pvalue.Cutoff, col='red', lwd=2)
    legend('topleft', legend=c('E-value', 'p-value', 'q-value'), col=c('blue', 'darkgreen','darkred'), lwd=2,bg='white',bty='o')
  }

}
plot_corrections(kruskal.wallis.table, pvalue.threshold)
```

# Annontate restults

Use the script `make_gene_info_summaries.py` to get gene information summaries that are written next to the results files.

# Plot random forest results

```{r}
# name = "beta-Alanine metabolism [PATH:ko00410]"
remove_fluff_ko_C <- function(name){
  if (grepl("\\[", name)) {
    return(strsplit(name, "\\ \\[")[[1]][1])
  } else {
    return(name)
  }
}

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv")
max_rank <- max(df_kos_ann$rank)
df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                filter(rank < 500) %>%
                  group_by(ko_B, ko_A) %>%
                  summarise(mean_value = mean(Value))
ko_b_order <- df_kos_ann %>% arrange(mean_value) %>% pull(ko_B) 
ggplot(df_kos_ann, aes(y = factor(ko_B, ko_b_order), x = mean_value, fill = ko_A))+
  # geom_boxplot(aes(y = ko_B, x = Value)) +
  # geom_point(aes(y = ko_B, x = Value)) +
  labs(x = "Average mean decrease in accuracy", y = "KO Category B (Top 500)") +
  geom_bar(stat = "identity") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=T, guide_nrow = 10)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosB_Top500.pdf", width = 10, height = 10, dpi = 300)


df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv")
max_rank <- max(df_kos_ann$rank)
df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                filter(!is.na(ko_B)) %>%
                  group_by(ko_B, ko_A) %>%
                  summarise(mean_value = mean(Value))
ko_b_order <- df_kos_ann %>% arrange(mean_value) %>% pull(ko_B) 
ggplot(df_kos_ann, aes(y = factor(ko_B, ko_b_order), x = mean_value, fill = ko_A))+
  # geom_boxplot(aes(y = ko_B, x = Value)) +
  # geom_point(aes(y = ko_B, x = Value)) +
  labs(x = "Average mean decrease in accuracy", y = "KO Category B") +
  geom_bar(stat = "identity") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=T, guide_nrow = 10)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosB_All.pdf", width = 10, height = 10, dpi = 300)

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                filter(rank < 500) %>%
                  group_by(ko_B)
ggplot(df_kos_ann)+
  geom_boxplot(aes(y = ko_B, x = Value, fill = ko_A)) +
  geom_point(aes(y = ko_B, x = Value)) +
  labs(x = "Mean decrease in accuracy", y = "KO Category B (Top 500)") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=T, guide_nrow = 10)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosB_boxplot_Top500.pdf", width = 10, height = 10, dpi = 300)

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                  group_by(ko_B) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 1)
ggplot(df_kos_ann)+
  geom_boxplot(aes(y = ko_B, x = Value, fill = ko_A)) +
  geom_point(aes(y = ko_B, x = Value)) +
  labs(x = "Mean decrease in accuracy", y = "KO Category B (Top 500)") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=T, guide_nrow = 10)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosB_boxplot_Ngrtr1.pdf", width = 15, height = 15, dpi = 300)

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                filter(rank < 500) %>%
                  group_by(ko_C) %>%
                  mutate(ko_C = Vectorize(remove_fluff_ko_C)(ko_C)) %>%
                  summarize(mean_value = mean(Value)) %>%
                  unique()
ko_c_order <- df_kos_ann %>% arrange(mean_value) %>% pull(ko_C) 
ggplot(df_kos_ann, aes(y = factor(ko_C, ko_c_order), x = mean_value))+
  # geom_boxplot(aes(y = ko_C, x = Value)) +
  # geom_point(aes(y = ko_C, x = Value)) +
  labs(x = "Average mean decrease in accuracy", y = "KO Category C (Top 500)") +
  geom_bar(stat = "identity") +
  make_theme(y_size = 10, y_hj = 1)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosC_Top500.pdf", width = 15, height = 15, dpi = 300)

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                # filter(rank < 500) %>%
                  group_by(ko_C) %>%
                  mutate(ko_C = Vectorize(remove_fluff_ko_C)(ko_C)) %>%
                  summarize(mean_value = mean(Value)) %>%
                  unique()
ko_c_order <- df_kos_ann %>% arrange(mean_value) %>% pull(ko_C) 
ggplot(df_kos_ann, aes(y = factor(ko_C, ko_c_order), x = mean_value))+
  # geom_boxplot(aes(y = ko_C, x = Value)) +
  # geom_point(aes(y = ko_C, x = Value)) +
  labs(x = "Average mean decrease in accuracy", y = "KO Category C") +
  geom_bar(stat = "identity") +
  make_theme(y_size = 10, y_hj = 1)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosC_All.pdf", width = 10, height = 22, dpi = 300)

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                filter(rank < 500) %>%
                  group_by(ko_C) %>%
                  mutate(ko_C = Vectorize(remove_fluff_ko_C)(ko_C)) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 1)
ko_c_order <- df_kos_ann %>% group_by(ko_C) %>% mutate(mean_value = mean(Value)) %>% arrange(mean_value) %>% pull(ko_C)
ggplot(df_kos_ann)+
  geom_boxplot(aes(y = ko_C, x = Value, fill = ko_A)) +
  geom_point(aes(y = ko_C, x = Value)) +
  labs(x = "Mean decrease in accuracy", y = "KO Category C (Top 500 and group with N > 2)") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=T, guide_nrow = 5)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosC_boxplot_Top500Ngtr1.pdf", width = 10, height = 10, dpi = 300)

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                  group_by(ko_C) %>%
                  mutate(ko_C = Vectorize(remove_fluff_ko_C)(ko_C)) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 1)
# ko_c_order <- df_kos_ann %>% group_by(ko_C) %>% mutate(mean_value = mean(Value)) %>% arrange(mean_value) %>% pull(ko_C)
ggplot(df_kos_ann)+
  geom_boxplot(aes(y = ko_C, x = Value, fill = ko_A)) +
  geom_point(aes(y = ko_C, x = Value)) +
  labs(x = "Mean decrease in accuracy", y = "KO Category C (Top 500 and group with N > 1)") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=T, guide_nrow = 10)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosC_boxplot_Ngtr1.pdf", width = 10, height = 22, dpi = 300)

df_kos_ann <- read.csv("results/figures/08-summarize_functions/KOs_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_rank) %>%
                  group_by(ko_C) %>%
                  mutate(ko_C = Vectorize(remove_fluff_ko_C)(ko_C)) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 5)
ko_c_order <- df_kos_ann %>% group_by(ko_C) %>% mutate(mean_value = mean(Value)) %>% arrange(mean_value) %>% pull(ko_C)
ggplot(df_kos_ann)+
  geom_boxplot(aes(y = ko_C, x = Value, fill = ko_A)) +
  geom_point(aes(y = ko_C, x = Value)) +
  labs(x = "Mean decrease in accuracy", y = "KO Category C (Top 500 and group with N > 5)") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=T, guide_nrow = 10)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_kosC_boxplot_Ngtr5.pdf", width = 10, height = 20, dpi = 300)

df_cazy_ann <- read.csv("results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_importance_annotated.csv")
max_rank <- max(df_cazy_ann$rank)
df_cazy_ann <- read.csv("results/figures/08-summarize_functions/cazyme_compared_kruskal_wallis_importance_annotated.csv") %>%
                mutate(cazy = Sample) %>%
                filter(Glycan_annotation != "") %>%
                arrange(rank) %>%
                filter(rank < max_rank)
cazy_colors <- c("Unknown" = "#f2f2f2",
                 "PG" = "#a6cee3",
                 "Amylose_and_Amylopectin" = "#fbb4ae",
                 "Amylose_and_Amylopectin,Gum,Pectin" = "#bdbdbdbd",
                 "Amylose_and_Amylopectin,Gum,Hemicellulose,Mannans_and_Heteromannans,Pectin" = "#bdbdbdbd",
                 "Amylose_and_Amylopectin,Pectin" = "#bdbdbdbd",
                 "Amylose_and_Amylopectin,Gum,Mannans_and_Heteromannans,Pectin" = "#bdbdbdbd",
                 "Gum" = "#fdbf6f",
                 "Gum,Hemicellulose,Mannans_and_Heteromannans,Pectin" = "#bdbdbdbd",
                 "Gum,Glycoprotein,Hemicellulose,Pectin" = "#bdbdbdbd",
                 "Gum,Glycoprotein,Pectin" = "#bdbdbdbd",
                 "Gum,Pectin" = "#bdbdbdbd",
                 "Hemicellulose,Mannans_and_Heteromannans" = "#bdbdbdbd",
                 "Hemicellulose,Pectin" = "#bdbdbdbd",
                 "Hemicellulose" = "#ff7f00",
                 "Pectin" = "#e31a1c",
                 "Cellulose" = "#1f78b4",
                 "Cellulose,Gum,Glycoprotein,Hemicellulose,Mannans_and_Heteromannans,Pectin" = "#bdbdbdbd",
                 "Cellulose,Hemicellulose,Mannans_and_Heteromannans" = "#bdbdbdbd",
                 "Glycoprotein,Gum,Hemicellulose,Mannans_and_Heteromannans,Pectin" = "#bdbdbdbd"
)
cazy_order <- df_cazy_ann %>% arrange(Value) %>% pull(cazy)
ggplot(df_cazy_ann, aes(y = factor(cazy, cazy_order), x = Value, fill = Function_at_destination_3))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cazy_colors) +
  labs(x = "Mean decrease in accuracy", y = "CAZyme") +
  make_theme(y_size = 10, y_hj = 1, leg_pos = "right", setFill=F, guide_nrow = 20, leg_size = 12)
  ggsave("results/figures/08-summarize_functions/kruskal_wallis_importance_cazy_All.pdf", width = 15, height = 10, dpi = 300)


get_genus_from_og <- function(og_species) {
  og_species = str_split(og_species, " ")[[1]][1]
  if (grepl("g__", og_species)) {
    return(strsplit(og_species, "g__")[[1]][2])
  } else {
    return(og_species)
  }
}

df_og_ann <- read.csv("results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance_annotated.csv")
max_og_rank <- max(df_og_ann$rank)
df_og_ann <- read.csv("results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_og_rank) %>%
                mutate(genus = Vectorize(get_genus_from_og)(og_species)) %>%
                mutate(genus = ifelse(genus %in% names(genusColorsClean), genus, NA)) %>%
                filter(rank < max_og_rank) %>%
                  group_by(kos_C, genus) %>%
                  # mutate(num_in_group = n()) %>%
                  # filter(num_in_group > 1) %>%
                  filter(rank < 500) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  filter(kos_C != "") %>%
                  summarize(genus, mean_value = mean(Value)) %>%
                  unique()
kos_c_order <- df_og_ann %>% group_by(kos_C) %>% summarize(mean_value_sum = sum(mean_value)) %>% arrange(mean_value_sum) %>% pull(kos_C) %>% unique
ggplot(df_og_ann, aes(y = factor(kos_C, kos_c_order), x = mean_value, fill = genus))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genusColorsClean) +
  labs(x = "Mean decrease in accuracy", y = "OG and KO Category C (Top 500)") +
  make_theme(setFill = F, y_size = 10, y_hj = 1,
              leg_pos = "right", guide_nrow = 22)
df_og_ann <- read.csv("results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_og_rank) %>%
                mutate(genus = Vectorize(get_genus_from_og)(og_species)) %>%
                mutate(genus = ifelse(genus %in% names(genusColorsClean), genus, NA)) %>%
                  group_by(kos_C, genus) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 1) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  summarize(genus, mean_value = mean(Value)) %>%
                  filter(kos_C != "") %>%
                  unique()
kos_c_order <- df_og_ann %>% group_by(kos_C) %>% summarize(mean_value_sum = sum(mean_value)) %>% arrange(mean_value_sum) %>% pull(kos_C) %>% unique
ggplot(df_og_ann, aes(y = factor(kos_C, kos_c_order), x = mean_value, fill = genus))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genusColorsClean) +
  labs(x = "Mean decrease in accuracy", y = "OG and KO Category C") +
  make_theme(setFill = F, y_size = 10, y_hj = 1,
              leg_pos = "right", guide_nrow = 22)

df_og_ann <- read.csv("results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(coreness == "core") %>%
                filter(rank < max_og_rank) %>%
                mutate(genus = Vectorize(get_genus_from_og)(og_species)) %>%
                mutate(genus = ifelse(genus %in% names(genusColorsClean), genus, NA)) %>%
                  group_by(kos_C, genus) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 1) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  summarize(genus, mean_value = mean(Value)) %>%
                  filter(kos_C != "") %>%
                  unique()
kos_c_order <- df_og_ann %>% group_by(kos_C) %>% summarize(mean_value_sum = sum(mean_value)) %>% arrange(mean_value_sum) %>% pull(kos_C) %>% unique
ggplot(df_og_ann, aes(y = factor(kos_C, kos_c_order), x = mean_value, fill = genus))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genusColorsClean) +
  labs(x = "Mean decrease in accuracy", y = "OG and KO Category C (core)") +
  make_theme(setFill = F, y_size = 10, y_hj = 1,
              leg_pos = "right", guide_nrow = 22)

df_og_ann <- read.csv("results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(coreness != "core") %>%
                filter(rank < max_og_rank) %>%
                mutate(genus = Vectorize(get_genus_from_og)(og_species)) %>%
                mutate(genus = ifelse(genus %in% names(genusColorsClean), genus, NA)) %>%
                  group_by(kos_C, genus) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 1) %>%
                  # mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  summarize(genus, mean_value = mean(Value)) %>%
                  filter(kos_C != "") %>%
                  unique()
kos_c_order <- df_og_ann %>% group_by(kos_C) %>% summarize(mean_value_sum = sum(mean_value)) %>% arrange(mean_value_sum) %>% pull(kos_C) %>% unique
ggplot(df_og_ann, aes(y = factor(kos_C, kos_c_order), x = mean_value, fill = genus))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genusColorsClean) +
  labs(x = "Mean decrease in accuracy", y = "OG and KO Category C (NOT core)") +
  make_theme(setFill = F, y_size = 10, y_hj = 1,
              leg_pos = "right", guide_nrow = 22)

df_og_ann <- read.csv("results/figures/08-summarize_functions/og_compared_kruskal_wallis_importance_annotated.csv") %>%
                arrange(rank) %>%
                filter(rank < max_og_rank) %>%
                mutate(genus = Vectorize(get_genus_from_og)(og_species)) %>%
                mutate(genus = ifelse(genus %in% names(genusColorsClean), genus, NA)) %>%
                  group_by(kos_C, genus) %>%
                  mutate(num_in_group = n()) %>%
                  filter(num_in_group > 1) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  summarize(genus, mean_value = mean(Value)) %>%
                  filter(kos_C != "") %>%
                  filter(genus == "Bifidobacterium") %>%
                  unique()
kos_c_order <- df_og_ann %>% group_by(kos_C) %>% summarize(mean_value_sum = sum(mean_value)) %>% arrange(mean_value_sum) %>% pull(kos_C) %>% unique
ggplot(df_og_ann, aes(y = factor(kos_C, kos_c_order), x = mean_value, fill = genus))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genusColorsClean) +
  labs(x = "Mean decrease in accuracy", y = "OG and KO Category C") +
  make_theme(setFill = F, y_size = 10, y_hj = 1,
              leg_pos = "right", guide_nrow = 22)

```

# Plot Maaslin2 pairwise results

```{r}
system("mkdir results/figures/08_plot_functional_comparisons")
system("mkdir results/figures/08_plot_functional_comparisons/masslin2_KO")
df_ko_maaslin <- read.csv("results/figures/08-summarize_functions/maaslin2_results_ko_pairwise_annotated.tsv", sep = "\t")
colnames(df_ko_maaslin)
all_pairs <- df_ko_maaslin$contrast %>% unique()
for (pair in all_pairs) {
  print(pair)
  # pair <- all_pairs[1]
  df_ko_subset <- df_ko_maaslin %>% filter(contrast == pair) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C))
  # make the volcano plot
  ggplot(df_ko_subset, aes(x = estimate, y = -log10(p.value), color = ifelse(p.value < 0.05, "significant", "not significant"))) +
    geom_point() +
    geom_text_repel(aes(label = kos_C), box.padding = 0.05, point.padding = 0.05, segment.color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("significant" = "red", "not significant" = "grey")) +
    labs(x = "Effect size", y = "-log10(p-value)", color = "Significane") +
    ggtitle(pair) +
    make_theme(setCol = F)
    ggsave(paste0("results/figures/08_plot_functional_comparisons/masslin2_KO/maaslin2_ko_", pair, ".pdf"))

  # summarize which categories are significantly different for each pair
  # as a bar plot with y axis as kos_C and x axis as the number of significant differences
  df_ko_subset_more  <- df_ko_maaslin %>% filter(contrast == pair) %>%
                  mutate(pair1 = str_split(contrast, " - ") %>% sapply(function(x) x[1])) %>%
                  mutate(pair2 = str_split(contrast, " - ") %>% sapply(function(x) x[2])) %>%
                  mutate(host = ifelse(estimate > 0, paste(pair1), paste(pair2))) %>%
                  filter(p.value < 0.05)
  ko_b_order <- df_ko_subset_more %>% group_by(kos_B) %>% summarise(n = n()) %>% unique %>% arrange(n) %>% pull(kos_B) %>% as.character
  ggplot(df_ko_subset_more) +
    geom_point(aes(x = estimate, y = factor(kos_B, ko_b_order), color = host)) +
    make_theme(y_size = 10, y_hj = 1, setCol = F) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = host_order_color_dark) +
    ggtitle(pair) +
    labs(x = "Effect size", y = "KO Category")
    ggsave(paste0("results/figures/08_plot_functional_comparisons/masslin2_KO/maaslin2_koB_", pair, "_host_compared_All.pdf"), width = 15, height = 15, dpi = 300)
  df_ko_subset_more  <- df_ko_maaslin %>% filter(contrast == pair) %>%
                  mutate(pair1 = str_split(contrast, " - ") %>% sapply(function(x) x[1])) %>%
                  mutate(pair2 = str_split(contrast, " - ") %>% sapply(function(x) x[2])) %>%
                  mutate(host = ifelse(estimate > 0, paste(pair1), paste(pair2))) %>%
                  filter(p.value < 0.05) %>%
                  group_by(kos_B) %>%
                  mutate(host, n = n()) %>%
                  filter(n >= 10)
  ko_b_order <- df_ko_subset_more %>% group_by(kos_B) %>% summarise(n = n()) %>% unique %>% arrange(n) %>% pull(kos_B) %>% as.character
  ggplot(df_ko_subset_more) +
    geom_point(aes(x = estimate, y = factor(kos_B, ko_b_order), color = host)) +
    make_theme(y_size = 10, y_hj = 1, setCol = F) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = host_order_color_dark) +
    ggtitle(pair) +
    geom_text_repel(data = df_ko_subset_more %>% filter(estimate > 1), aes(label = feature, x = estimate, y = factor(kos_B, ko_b_order), color = host ), hjust = 1.5, vjust = 0.5) +
    geom_text_repel(data = df_ko_subset_more %>% filter(estimate < -1), aes(label = feature, x = estimate, y = factor(kos_B, ko_b_order), color = host ), hjust = 1.5, vjust = 0.5) +
    labs(x = "Effect size", y = "KO Category")
    ggsave(paste0("results/figures/08_plot_functional_comparisons/masslin2_KO/maaslin2_koB_", pair, "_host_compared_Ngeq10.pdf"), width = 15, height = 15, dpi = 300)
  df_ko_subset_more  <- df_ko_maaslin %>% filter(contrast == pair) %>%
                  mutate(pair1 = str_split(contrast, " - ") %>% sapply(function(x) x[1])) %>%
                  mutate(pair2 = str_split(contrast, " - ") %>% sapply(function(x) x[2])) %>%
                  mutate(host = ifelse(estimate > 0, paste(pair1), paste(pair2))) %>%
                  filter(p.value < 0.05) %>%
                  group_by(kos_B) %>%
                  mutate(host, n = n()) %>%
                  filter(n < 10)
  ko_b_order <- df_ko_subset_more %>% group_by(kos_B) %>% summarise(n = n()) %>% unique %>% arrange(n) %>% pull(kos_B) %>% as.character
  ggplot(df_ko_subset_more) +
    geom_point(aes(x = estimate, y = factor(kos_B, ko_b_order), color = host)) +
    make_theme(y_size = 10, y_hj = 1, setCol = F) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = host_order_color_dark) +
    ggtitle(pair) +
    labs(x = "Effect size", y = "KO Category") +
    geom_text_repel(aes(label = feature, x = estimate, y = factor(kos_B, ko_b_order), color = host ), hjust = 1.5, vjust = 0.5)
    ggsave(paste0("results/figures/08_plot_functional_comparisons/masslin2_KO/maaslin2_koB_", pair, "_host_compared_Nless10.pdf"), width = 15, height = 15, dpi = 300)
  
  df_ko_subset_more  <- df_ko_maaslin %>% filter(contrast == pair) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  mutate(pair1 = str_split(contrast, " - ") %>% sapply(function(x) x[1])) %>%
                  mutate(pair2 = str_split(contrast, " - ") %>% sapply(function(x) x[2])) %>%
                  mutate(host = ifelse(estimate > 0, paste(pair1), paste(pair2))) %>%
                  filter(p.value < 0.05) %>%
                  group_by(kos_C)
  ko_c_order <- df_ko_subset_more %>% group_by(kos_C) %>% summarise(n = n()) %>% unique %>% arrange(n) %>% pull(kos_C) %>% as.character
  ggplot(df_ko_subset_more) +
    geom_point(aes(x = estimate, y = factor(kos_C, ko_c_order), color = host)) +
    make_theme(y_size = 10, y_hj = 1, setCol = F) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = host_order_color_dark) +
    labs(x = "Effect size", y = "KO Category") +
    ggtitle(pair)
  ggsave(paste0("results/figures/08_plot_functional_comparisons/masslin2_KO/maaslin2_koC_", pair, "_host_compared_All.pdf"), width = 15, height = 15, dpi = 300)
  df_ko_subset_more  <- df_ko_maaslin %>% filter(contrast == pair) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  mutate(pair1 = str_split(contrast, " - ") %>% sapply(function(x) x[1])) %>%
                  mutate(pair2 = str_split(contrast, " - ") %>% sapply(function(x) x[2])) %>%
                  mutate(host = ifelse(estimate > 0, paste(pair1), paste(pair2))) %>%
                  filter(p.value < 0.05) %>%
                  group_by(kos_C) %>%
                  mutate(host, n = n()) %>%
                  filter(n >= 10)
  ko_c_order <- df_ko_subset_more %>% group_by(kos_C) %>% summarise(n = n()) %>% unique %>% arrange(n) %>% pull(kos_C) %>% as.character
  ggplot(df_ko_subset_more) +
    geom_point(aes(x = estimate, y = factor(kos_C, ko_c_order), color = host)) +
    make_theme(y_size = 10, y_hj = 1, setCol = F) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = host_order_color_dark) +
    labs(x = "Effect size", y = "KO Category") +
    geom_text_repel(data = df_ko_subset_more %>% filter(estimate > 2), aes(label = feature, x = estimate, y = factor(kos_C, ko_c_order), color = host ), hjust = 1.5, vjust = 0.5) +
    geom_text_repel(data = df_ko_subset_more %>% filter(estimate < -2), aes(label = feature, x = estimate, y = factor(kos_C, ko_c_order), color = host ), hjust = 1.5, vjust = 0.5) +
    ggtitle(pair)
  ggsave(paste0("results/figures/08_plot_functional_comparisons/masslin2_KO/maaslin2_koC_", pair, "_host_compared_Ngeq10.pdf"), width = 15, height = 15, dpi = 300)
  df_ko_subset_more  <- df_ko_maaslin %>% filter(contrast == pair) %>%
                  mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                  mutate(pair1 = str_split(contrast, " - ") %>% sapply(function(x) x[1])) %>%
                  mutate(pair2 = str_split(contrast, " - ") %>% sapply(function(x) x[2])) %>%
                  mutate(host = ifelse(estimate > 0, paste(pair1), paste(pair2))) %>%
                  filter(p.value < 0.05) %>%
                  group_by(kos_C) %>%
                  mutate(host, n = n()) %>%
                  filter(n < 10)
  ko_c_order <- df_ko_subset_more %>% group_by(kos_C) %>% summarise(n = n()) %>% unique %>% arrange(n) %>% pull(kos_C) %>% as.character
  ggplot(df_ko_subset_more) +
    geom_point(aes(x = estimate, y = factor(kos_C, ko_c_order), color = host)) +
    make_theme(y_size = 10, y_hj = 1, setCol = F) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = host_order_color_dark) +
    ggtitle(pair) +
    geom_text_repel(data = df_ko_subset_more %>% filter(estimate > 1), aes(label = feature, x = estimate, y = factor(kos_C, ko_c_order), color = host ), hjust = 1.5, vjust = 0.5) +
    geom_text_repel(data = df_ko_subset_more %>% filter(estimate < -1.5), aes(label = feature, x = estimate, y = factor(kos_C, ko_c_order), color = host ), hjust = 1.5, vjust = 0.5) +
    labs(x = "Effect size", y = "KO Category")
  ggsave(paste0("results/figures/08_plot_functional_comparisons/masslin2_KO/maaslin2_koC_", pair, "_host_compared_Nless10.pdf"), width = 15, height = 15, dpi = 300)
}

df_cazy_maaslin <- read.csv("results/figures/08-summarize_functions/maaslin2_results_cazyme_pairwise_annotated.csv")
colnames(df_cazy_maaslin)
all_pairs <- df_cazy_maaslin$contrast %>% unique()
pair <- all_pairs[1]
for (pair in all_pairs) {
  df_cazy_subset <- df_cazy_maaslin %>% filter(contrast == pair)
  # make the volcano plot
  ggplot(df_cazy_subset, aes(x = estimate, y = -log10(p.value), color = ifelse(p.value < 0.05, "significant", "not significant"))) +
    geom_point() +
    geom_text_repel(aes(label = paste0(feature, " - ", Function_at_destination_3)), box.padding = 0.5, point.padding = 0.5, segment.color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("significant" = "red", "not significant" = "grey")) +
    labs(x = "Effect size", y = "-log10(p-value)", color = "Significane") +
    ggtitle(pair) +
    make_theme(setCol = F)
  ggsave(paste0("results/figures/08_plot_functional_comparisons/maaslin2_cazyme_", pair, ".pdf"))
}


df_og_maaslin <- read.csv("results/figures/08-summarize_functions/maaslin2_results_og_pairwise_annotated.csv")
colnames(df_og_maaslin)
all_pairs_og <- df_og_maaslin$contrast %>% unique()
pair <- all_pairs_og[1]
df_og_subset <- df_og_maaslin %>% filter(contrast == pair) %>%
                mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                mutate(genus = Vectorize(get_genus_from_og)(og_species))
# make the volcano plot
for (genus_iter in df_og_subset$genus %>% unique()) {
  df_og_subset_genus <- df_og_subset %>% filter(genus == genus_iter)
  ggplot(df_og_subset_genus, aes(x = estimate, y = -log10(p.value), color = ifelse(p.value < 0.05, "significant", "not significant"))) +
    geom_point() +
    geom_text_repel(aes(label = kos_C), box.padding = 0.05, point.padding = 0.05, segment.color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("significant" = "red", "not significant" = "grey")) +
    labs(x = "Effect size", y = "-log10(p-value)", color = "Significane") +
    ggtitle(paste(pair, genus_iter)) +
    make_theme(setCol = F)
  ggsave(paste0("results/figures/08_plot_functional_comparisons/maaslin2_og_", pair, "_", genus_iter, ".pdf"))
}
# summarize for each pair in each genus which kos_C are significantly different
# as a bar plot with y axis as kos_C and x axis as the number of significant differences
# for each genus with the color as the honeybee species in which it is found more 
# (+ve means 1st in pair and -ve means 2nd in pair)
df_og_subset_more <- df_og_maaslin %>% filter(contrast == pair) %>%
                mutate(kos_C = Vectorize(remove_fluff_ko_C)(kos_C)) %>%
                mutate(genus = Vectorize(get_genus_from_og)(og_species)) %>%
                filter(p.value < 0.05) %>%
                mutate(pair1 = str_split(contrast, " - ") %>% sapply(function(x) x[1])) %>%
                mutate(pair2 = str_split(contrast, " - ") %>% sapply(function(x) x[2])) %>%
                mutate(host = ifelse(estimate > 0, paste(pair1), paste(pair2))) %>%
                group_by(kos_C, genus) %>%
                summarize(host, n = n()) %>%
                filter(n > 10) %>%
                filter(kos_C != "") %>%
                arrange(n) %>%
                filter(!is.na(kos_C))
ko_c_order <- df_og_subset_more %>% pull(kos_C) %>% unique
ggplot(df_og_subset_more, aes(x = n, y = factor(kos_C, ko_c_order), fill = host)) +
  geom_bar(stat = "identity") +
  facet_wrap(~genus) +
  scale_fill_manual(values = host_order_color) +
  ggtitle(pair) +
  scale_x_continuous(trans = "log") +
  make_theme(y_size = 10, x_size = 10, y_hj = 1, setFill = F)

```

```{r}
```
# Module completeness using MicrobeAnnotator results

```{r}
df_by_sample <- read.csv("results/figures/08-summarize_functions/microbeannotator_out/by_sample_module_completeness.tab", check.names = F, sep = "\t") %>%
                  mutate(name = abbreviate(name, 55))
format_colnames <- function(list_cols) {
  sapply(list_cols, function(x) if (x %in% c("name", "module", "pathway group")) {return(x)} else {return(strsplit(x, "_kos.csv")[[1]][1])})
}
colnames(df_by_sample) <- df_by_sample %>% colnames %>% format_colnames
by_sample_mat <- df_by_sample %>% 
  select(!c(module, `pathway group`), ) %>%
  column_to_rownames("name") %>% as.matrix
# remove sparse rows
by_sample_mat <- by_sample_mat[apply(by_sample_mat, 1, function(x) sum(x > 0, na.rm = T)) > 10,]
anno_top_host <- HeatmapAnnotation(Host = by_sample_mat %>% colnames %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                              col = list(Host = host_order_color_dark
                              ), border = T
                          )
anno_left_module <- rowAnnotation(Group = by_sample_mat %>% rownames %>% as.data.frame() %>% left_join(df_by_sample %>% select(name, `pathway group`), by = c("." = "name")) %>% pull(`pathway group`),
                              col = list(Group = make_col_list(by_sample_mat %>% rownames %>% as.data.frame() %>% left_join(df_by_sample %>% select(name, `pathway group`), by = c("." = "name")) %>% pull(`pathway group`), "Set3")
                              ), border = T,
                              show_legend = c(F)
                          )
row_split <- by_sample_mat %>% rownames %>% as.data.frame() %>% left_join(df_by_sample %>% select(name, `pathway group`), by = c("." = "name")) %>% pull(`pathway group`) %>% as.factor
column_split <- by_sample_mat %>% colnames %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host) %>% as.factor
col_fun = colorRamp2(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), rev(colorRampPalette(c(brewer.pal(9, "Blues")[9], "#FFFFFF"))(11)))
hp <- Heatmap(by_sample_mat,
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        column_split = column_split,
        row_split = row_split,
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "Completeness (%)"),
        top_annotation = anno_top_host,
        left_annotation = anno_left_module,
        # bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 5),
        cluster_columns = T,
        cluster_rows = T
        # column_dend_side = "bottom",
)
pdf("results/figures/11-gene_function/heatmap_by_sample_module_completeness.pdf", width = 15, height = 20)
draw(hp, heatmap_legend_side = "top", merge_legend = TRUE)
dev.off()
hp <- Heatmap(by_sample_mat,
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        # column_split = column_split,
        row_split = row_split,
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "Completeness (%)"),
        top_annotation = anno_top_host,
        left_annotation = anno_left_module,
        # bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 5),
        cluster_columns = T,
        cluster_rows = T
        # column_dend_side = "bottom",
)
pdf("results/figures/11-gene_function/heatmap_by_sample_module_completeness_clust.pdf", width = 15, height = 20)
draw(hp, heatmap_legend_side = "top", merge_legend = TRUE)
dev.off()
```

```{r}
# genus <- "g__Gilliamella"
for (genus in c(genera, "g__Saezia_snod")) {
  if (file.exists(paste0("results/figures/08-summarize_functions/microbeannotator_out/by_genus/", genus, "/", genus,"__module_completeness.tab"))) {
    if (1==1) {
    # if (!dir.exists(paste0("results/figures/11-gene_function/by_genus/", genus))) {
      dir.create(paste0("results/figures/11-gene_function/by_genus/", genus))
      df_by_genus <- read.csv(paste0("results/figures/08-summarize_functions/microbeannotator_out/by_genus/", genus, "/", genus,"__module_completeness.tab"), check.names = F, sep = "\t") %>%
                        mutate(name = abbreviate(name, 55))
      format_colnames <- function(list_cols) {
        sapply(list_cols, function(x) if (x %in% c("name", "module", "pathway group")) {return(x)} else {return(strsplit(x, "_kos.csv")[[1]][1])})
      }
      colnames(df_by_genus) <- df_by_genus %>% colnames %>% format_colnames
      by_genus_mat <- df_by_genus %>% 
        select(!c(module, `pathway group`), ) %>%
        column_to_rownames("name") %>% as.matrix
      # for each species / genus get a set of samples that have enough coverage of the species / genus to proceed
      # subset those columns here - for now only removing "samples_to_exclude"
      by_genus_mat <- by_genus_mat[, colnames(by_genus_mat) %in% samples]
      # remove sparse rows
      min <- 10
      if (dim(by_genus_mat)[[2]] < 10) {
        min <- 1
      }
      by_genus_mat <- by_genus_mat[apply(by_genus_mat, 1, function(x) sum(x > 0, na.rm = T)) > min,]
      anno_top_host <- HeatmapAnnotation(Host = by_genus_mat %>% colnames %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host),
                                    col = list(Host = host_order_color_dark
                                    ), border = T
                                )
      anno_left_module <- rowAnnotation(Group = by_genus_mat %>% rownames %>% as.data.frame() %>% left_join(df_by_genus %>% select(name, `pathway group`), by = c("." = "name")) %>% pull(`pathway group`),
                                    col = list(Group = make_col_list(by_genus_mat %>% rownames %>% as.data.frame() %>% left_join(df_by_genus %>% select(name, `pathway group`), by = c("." = "name")) %>% pull(`pathway group`), "Set3")
                                    ), border = T,
                                    show_legend = c(F)
                                )
      row_split <- by_genus_mat %>% rownames %>% as.data.frame() %>% left_join(df_by_genus %>% select(name, `pathway group`), by = c("." = "name")) %>% pull(`pathway group`) %>% as.factor
      column_split <- by_genus_mat %>% colnames %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(.)) %>% pull(Host) %>% as.factor
      col_fun = colorRamp2(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), rev(colorRampPalette(c(brewer.pal(9, "Blues")[9], "#FFFFFF"))(11)))
      hp <- Heatmap(by_genus_mat,
              col = col_fun,
              clustering_method_columns = 'ward.D2',
              column_split = column_split,
              row_split = row_split,
              column_title_gp = grid::gpar(fontsize = 0),
              row_title_gp = grid::gpar(fontsize = 12),
              row_title_rot = 0,
              # border = T,
              heatmap_legend_param = list(title = "Completeness (%)"),
              top_annotation = anno_top_host,
              left_annotation = anno_left_module,
              # bottom_annotation = column_ba, 
              # right_annotation = row_ra,
              column_names_gp = grid::gpar(fontsize = 0),
              row_names_gp = grid::gpar(fontsize = 5),
              cluster_columns = T,
              cluster_rows = T
              # column_dend_side = "bottom",
      )
      draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0(genus, " module completeness"))
      pdf(paste0("results/figures/11-gene_function/by_genus/", genus, "/", genus, "_module_completeness.pdf"), width = 15, height = 20)
      draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0(genus, " module completeness"))
      dev.off()
      hp <- Heatmap(by_genus_mat,
              col = col_fun,
              clustering_method_columns = 'ward.D2',
              # column_split = column_split,
              row_split = row_split,
              column_title_gp = grid::gpar(fontsize = 0),
              row_title_gp = grid::gpar(fontsize = 12),
              row_title_rot = 0,
              # border = T,
              heatmap_legend_param = list(title = "Completeness (%)"),
              top_annotation = anno_top_host,
              left_annotation = anno_left_module,
              # bottom_annotation = column_ba, 
              # right_annotation = row_ra,
              column_names_gp = grid::gpar(fontsize = 0),j
              row_names_gp = grid::gpar(fontsize = 5),
              cluster_columns = T,
              cluster_rows = T
              # column_dend_side = "bottom",
      )
      draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0(genus, " module completeness"))
      pdf(paste0("results/figures/11-gene_function/by_genus/", genus, "/", genus, "_module_completeness_clust.pdf"), width = 15, height = 20)
      draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0(genus, " module completeness"))
      dev.off()
    } else {
      print(paste0("directory already exists for ", genus))
    }
  } else {
    print(paste0("no file for ", genus, " at ", paste0("results/figures/08-summarize_functions/microbeannotator_out/by_genus/", genus, "/", genus,"__module_completeness.tab")))
  }
}
```

## For each species

```{r}
get_species_of_mag_split <- function(name) {
  splits = strsplit(name, "_")[[1]]
  spec_name = paste0(splits[[1]], " ", splits[[2]])
  return(spec_name)
}
dir.create("results/figures/11-gene_function/by_mag")

df_by_mag <- read.csv(paste0("results/figures/08-summarize_functions/microbeannotator_out/by_mag_bombi_module_completeness.tab"), check.names = F, sep = "\t") %>%
                  mutate(name = abbreviate(name, 55))
format_colnames <- function(list_cols) {
  sapply(list_cols, function(x) if (x %in% c("name", "module", "pathway group")) {return(x)} else {return(strsplit(x, "_kos.csv")[[1]][1])})
}
colnames(df_by_mag) <- df_by_mag %>% colnames %>% format_colnames
by_mag_mat <- df_by_mag %>% 
        select(!c(module, `pathway group`), ) %>%
        column_to_rownames("name") %>% as.matrix

# remove sparse rows
min <- 10
if (dim(by_mag_mat)[[2]] < 10) {
    min <- 1
}
by_mag_mat <- by_mag_mat[apply(by_mag_mat, 1, function(x) sum(x > 0, na.rm = T)) > min,]
specs_list <- by_mag_mat %>% colnames %>% as.data.frame() %>% mutate(Spec = Vectorize(get_species_of_mag_split)(.)) %>% pull(Spec)
anno_top_spec <- HeatmapAnnotation(Spec = specs_list,
                              col = list(Host = make_col_list(specs_list, "Set3")
                              ), border = T
                          )
group_list <- by_mag_mat %>% rownames %>% as.data.frame() %>% left_join(df_by_mag %>% select(name, `pathway group`), by = c("." = "name")) %>% pull(`pathway group`)
anno_left_module <- rowAnnotation(Group = group_list,
                              col = list(Group = make_col_list(group_list, "Set3")
                              ), border = T,
                              show_legend = c(F)
                          )
row_split <- group_list %>% as.factor
column_split <- specs_list %>% as.factor
col_fun = colorRamp2(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), rev(colorRampPalette(c(brewer.pal(9, "Blues")[9], "#FFFFFF"))(11)))
hp <- Heatmap(by_mag_mat,
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        column_split = column_split,
        row_split = row_split,
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "Completeness (%)"),
        top_annotation = anno_top_spec,
        left_annotation = anno_left_module,
        # bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 5),
        cluster_columns = T,
        cluster_rows = T
        # column_dend_side = "bottom",
)
draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0("Bombilactobacillus_mellifer module completeness"))
pdf(paste0("results/figures/11-gene_function/by_mag/Bombilactobacillus_mellifer_module_completeness.pdf"), width = 15, height = 20)
draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0("Bombilactobacillus_mellifer module completeness"))
dev.off()


hp <- Heatmap(by_mag_mat,
        col = col_fun,
        clustering_method_columns = 'ward.D2',
        # column_split = column_split,
        row_split = row_split,
        column_title_gp = grid::gpar(fontsize = 0),
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_rot = 0,
        # border = T,
        heatmap_legend_param = list(title = "Completeness (%)"),
        top_annotation = anno_top_spec,
        left_annotation = anno_left_module,
        # bottom_annotation = column_ba, 
        # right_annotation = row_ra,
        column_names_gp = grid::gpar(fontsize = 0),
        row_names_gp = grid::gpar(fontsize = 5),
        cluster_columns = T,
        cluster_rows = T
        # column_dend_side = "bottom",
)
draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0("Bombilactobacillus_mellifer module completeness"))
pdf(paste0("results/figures/11-gene_function/by_mag/Bombilactobacillus_mellifer_module_completeness_clust.pdf"), width = 15, height = 20)
draw(hp, heatmap_legend_side = "bottom", merge_legend = TRUE, column_title = paste0("Bombilactobacillus_mellifer module completeness"))
dev.off()
```


# Comparing omega2 between different groups of KOs

```{r}
ko_set_files <- list.files("results/figures/08-summarize_functions/exploreKOs/bycategory/B", full.names = T)
#  if the table inside the files has <= 10 columns, drop them
ko_set_files <- ko_set_files[!sapply(ko_set_files, function(x) {
  read.table(x, sep = ",", header = T) %>% dim %>% .[2] <= 10
})]
# print the category name after "ko_matrix" and before ".csv" and the omega2 when comparing by host species
for (ko_file in ko_set_files) {
  ko_table <- t(read.table(ko_file, sep = ",", header = T, row.names = 1))
  dist_kos <- vegdist(t(ko_table), method = "bray")
  df_meta_info <- dist_kos %>% as.data.frame() %>% mutate(Host = Vectorize(get_host_from_sample_name)(rownames(.))) %>% mutate(Host = as.factor(Host))
  adonis_ko <- adonis2(dist_kos ~ Host, data = df_meta_info, permutations = 999)
  res <- adonis_OmegaSq(adonis_ko)
  print(paste0("Category: ", gsub("ko_matrix_", "", gsub(".csv", "", basename(ko_file))), " Omega2: ", res$omega2))
}
```



# heatmap with TMI
```{r}
genera_clean_sub_all <- c(
    "Lactobacillus",
    "Bifidobacterium",
    "Bombilactobacillus",
    "CALYQJ01",
    "Snodgrassella",
    "Gilliamella",
    "CALYQQ01",
    "Apibacter",
    "Bartonella_A",
    "Frischella",
    "Commensalibacter",
    "Apilactobacillus",
    # "Bombella",
    "Dysgonomonas",
    # "Enterobacter",
    "Pectinatus",
    # "Spiroplasma",
    # "Zymobacter",
    "Entomomonas",
    "Saezia",
    # "Parolsenella",
    "WRHT01",
    # "Cronobacter",
    # "Melissococcus",
    # "Rosenbergiella",
    "Pluralibacter"
    # "KACC-21507"
  )
# dist_matrix <- vegdist(rpkm_matrix, method = "robust.aitchison")
# matrix <- as.matrix(dist_matrix)
# dist <- as.dist(matrix)
# res_pcoa <- pcoa(dist)
# ev1 <- res_pcoa$vectors[,1]
# ev2 <- res_pcoa$vectors[,2]
# df_pcoa_new <- data.frame(cbind(ev1,ev2)) %>% rownames_to_column("Sample")
# rownames(df_pcoa_new) <- NULL
# df_pcoa_new <- left_join(df_pcoa_new, all_mag_metadata, by = c("Sample" = "ID"))
# perc_axis <- round(((res_pcoa$values$Relative_eig[c(1,2)])*100), digits=1)
# axis_x_title <- paste0("PCo1 (",perc_axis[1],"%)")
# axis_y_title <- paste0("PCo2 (",perc_axis[2],"%)")
# ggplot(df_pcoa_new) +
#     geom_mark_ellipse(aes(x = ev1, y = ev2, group = Genus),
#                         color = NA, fill = "#000000", geom = "polygon",
#                         label.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
#                         con.type = "straight",
#                         label.fontsize = 5, con.cap = unit(1, "mm"),
#                         alpha = 0.05, level = 0.95, size = 0.5) +
#     geom_point(aes(x = ev1,
#                    y = ev2,
#                    shape = factor(Host, host_order),
#                    size = factor(Host, host_order),
#                    fill = factor(Genus, genera_clean_sub_all),
#                    color = factor(Genus, genera_clean_sub_all))
#                 ) +
#     geom_text_repel(aes(x = ev1, y = ev2, label = ifelse(Representative == 1, MAG_species_name_final, "")),
#                       box.padding = 0.1, color = "#000000", size = 3) +
#     labs(x=axis_x_title, y = axis_y_title, color = "Genus", shape = "Host") +
#     facet_wrap(~factor(Host, host_order)) +
#     scale_color_manual(values=genusColorsClean) +
#     scale_fill_manual(values=genusColorsClean) +
#     scale_shape_manual(values=c(
#       "Apis mellifera" = 16,
#       "Apis cerana" = 16,
#       "Apis dorsata" = 16,
#       "Apis florea" = 16,
#       "Apis andreniformis" = 16
#     )) +
#     scale_size_manual(values=c(
#       "Apis mellifera" = 3,
#       "Apis cerana" = 3,
#       "Apis dorsata" = 3.5,
#       "Apis florea" = 3,
#       "Apis andreniformis" = 2
#     )) +
#     guides(fill = "none", size = "none") +
#     make_theme(setFill = F, setCol = F, 
#     leg_pos = "bottom", leg_size = 8, guide_nrow = 30,
#     x_size = 10, y_size = 10, modify_guide = F)
#     ggsave("results/figures/11-figures/tax_func_plots/all_kos_vitamin_metabolism.pdf", width = 15, height = 15)
df_rpkms_by_func_all_pa <- df_rpkms_by_func_all
df_rpkms_by_func_all_pa[df_rpkms_by_func_all_pa > 0] <- 1
# complex heatmap sep by species / host
genus = "Apibacter"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_C %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Apibacter.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


# complex heatmap sep by species / host
genus = "Snodgrassella"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_C %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Snodgrassella.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


# complex heatmap sep by species / host
genus = "Bombilactobacillus"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_B %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Bombilactobacillus.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()

# complex heatmap sep by species / host
genus = "Bifidobacterium"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_C %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Bifidobacterium.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()

# complex heatmap sep by species / host
genus = "Gilliamella"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_C %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Gilliamella.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


genus = "Lactobacillus"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_C %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Lactobacillus.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


genus = "Frischella"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_C %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Lactobacillus.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


genus = "Bartonella"
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all_pa %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_C %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            # col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = 8),
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Lactobacillus.pdf", width = 20, height = 40)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()


# complex heatmap sep by species / host
# genus = ""
mags_to_select <- all_mag_metadata %>%
                      filter(Quality != "low") %>%
                      # filter(Quality == "high") %>%
                      # filter(Completeness >= 95) %>%
                      # filter(Contamination < 4) %>%
                      # filter(Representative == "1") %>%
                      filter(Genus %in% c("Snodgrassella", "Saezia")) %>%
                      # filter(Genus == genus) %>%
                      pull(ID)
df_rpkms_by_func <- df_rpkms_by_func_all %>%
                        filter(mag %in% mags_to_select)
rpkm_matrix <- df_rpkms_by_func %>% column_to_rownames("mag") %>% as.matrix()
# remove genes found in < 5 samples
rpkm_matrix <- rpkm_matrix[, colSums(rpkm_matrix > 0) > 5]
rpkm_matrix <- rpkm_matrix[, !colnames(rpkm_matrix) %in% kos_to_exclude]
# # make a presence absence matrix
# rpkm_matrix_pa <- rpkm_matrix
# rpkm_matrix_pa[rpkm_matrix > 0] <- 1
species_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(MAG_species_name_final)
host_names = rownames(rpkm_matrix) %>% as.data.frame() %>%
                            left_join(all_mag_metadata %>% ungroup(), by = c(. = "ID")) %>%
                              pull(Host)
anno_host_spec_col = HeatmapAnnotation(Host = host_names,
                                       Species = species_names,
                                        col = list(Host = host_order_color,
                                                     Species = make_col_list(species_names)
                                                  ),
                                          border = T
                                  )
ko_names <- colnames(rpkm_matrix)
category_A <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(A)
category_B <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(B)
category_C <- colnames(rpkm_matrix) %>% as.data.frame() %>%
                left_join(ko_category_info, by = c(. = "ko")) %>%
                pull(C)
col_split = species_names %>% as.factor
row_split = category_B %>% as.factor
heatmap_obj = Heatmap(t(rpkm_matrix),
            # col = col_input,
            # col = colorRamp2(c(0, 1), c("#ffffff", "#08306b")),
            col = colorRamp2(c(0, 10, 100, 1000, 10000), c("#ffffff", "#c6dbef", "#6baed6", "#4292c6", "#08306b")),
            # top_annotation = top_annotation_obj,
            # bottom_annotation = bottom_annotation_obj,
            top_annotation = anno_host_spec_col,
            column_names_gp = grid::gpar(fontsize = 0),
            row_names_gp = grid::gpar(fontsize = 5),
            heatmap_legend_param = list(title = "Genes abundance (RPKM)"),
            # cell_fun = my_cell_fun,
            cluster_columns = T,
            column_split = col_split,
            row_split = row_split,
            column_title_rot = 90,
            row_title_rot = 0,
            border = T,
            # row_dend_reorder = T,
            cluster_rows = T
            )
pdf("results/figures/11-figures/tax_func_heatmaps/ko_heatmap_Snod_Saez.pdf", width = 20, height = 30)
draw(heatmap_obj, merge_legend = TRUE)
dev.off()
```