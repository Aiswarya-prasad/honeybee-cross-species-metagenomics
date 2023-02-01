##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)
source("scripts/visualization/00a-plot_DNA_concentrations.R")
source("scripts/visualization/00b-plot_qPCR_values.R")

##############
# files to be read
##############
# parse this later
# read.csv(paste0("02_HostMapping/", sample, "_coverage.tsv"), sep = "\t")

df_reads <- data.frame()
for (sample in samples) {
  if (file.exists(paste0("fastqc/raw/", sample, "_R1_fastqc.zip"))) {
    number_raw_R1 <- read.csv(unz(paste0("fastqc/raw/", sample, "_R1_fastqc.zip"), paste0(sample, "_R1_fastqc/fastqc_data.txt")), sep = "\t") %>%
                  filter(.[[1]] == "Total Sequences") %>%
                    pull() %>%
                      as.integer()
  } else {
    number_raw_R1 <- NA
  }
  if (file.exists(paste0("fastqc/raw/", sample, "_R2_fastqc.zip"))) {
    number_raw_R2 <- read.csv(unz(paste0("fastqc/raw/", sample, "_R2_fastqc.zip"), paste0(sample, "_R2_fastqc/fastqc_data.txt")), sep = "\t") %>%
                  filter(.[[1]] == "Total Sequences") %>%
                    pull() %>%
                      as.integer()
  } else {
    number_raw_R2 <- NA
  }
  raw_reads <- number_raw_R1 + number_raw_R2
  if (file.exists(paste0("fastqc/trim/", sample, "_R1_trim_fastqc.zip"))) {
    number_trimmed_R1 <- read.csv(unz(paste0("fastqc/trim/", sample, "_R1_trim_fastqc.zip"), paste0(sample, "_R1_trim_fastqc/fastqc_data.txt")), sep = "\t") %>%
                  filter(.[[1]] == "Total Sequences") %>%
                    pull() %>%
                      as.integer()
  } else {
    number_trimmed_R1 <- NA
  }
  if (file.exists(paste0("fastqc/trim/", sample, "_R2_trim_fastqc.zip"))) {
    number_trimmed_R2 <- read.csv(unz(paste0("fastqc/trim/", sample, "_R2_trim_fastqc.zip"), paste0(sample, "_R2_trim_fastqc/fastqc_data.txt")), sep = "\t") %>%
                  filter(.[[1]] == "Total Sequences") %>%
                    pull() %>%
                      as.integer()
  } else {
    number_trimmed_R2 <- NA
  }
  trimmed <- number_trimmed_R1 + number_trimmed_R2
  if (file.exists(paste0("09_MagDatabaseProfiling/MAGsDatabaseMapping/", sample, "_flagstat.tsv"))) {
    mapped_full_db <- read.csv(paste0("09_MagDatabaseProfiling/MAGsDatabaseMapping/", sample, "_flagstat.tsv"), sep = "\t") %>%
                        filter(.[[3]] == "with itself and mate mapped") %>%
                          pull(1) %>%
                            as.integer()
  } else {
    mapped_full_db <- NA
  }
  if (file.exists(paste0("09_MagDatabaseProfiling/SNVProfiling/Mapping/", sample, "_flagstat.tsv"))) {
    mapped_rep_db <- read.csv(paste0("09_MagDatabaseProfiling/SNVProfiling/Mapping/", sample, "_flagstat.tsv"), sep = "\t") %>%
                        filter(.[[3]] == "with itself and mate mapped") %>%
                          pull(1) %>%
                            as.integer()
  } else {
    mapped_rep_db <- NA
  }
  if (file.exists(paste0("09_MagDatabaseProfiling/MAGsDatabaseMapping/", sample, "_host_mapping_flagstat.tsv"))) {
    mapped_host_db <- read.csv(paste0("09_MagDatabaseProfiling/MAGsDatabaseMapping/", sample, "_host_mapping_flagstat.tsv"), sep = "\t") %>%
                        filter(.[[3]] == "with itself and mate mapped") %>%
                          pull(1) %>%
                            as.integer()
  } else {
    mapped_host_db <- NA
  }
  host_filtered <- trimmed - mapped_host_db
  unmapped_all <- trimmed - mapped_host_db - mapped_full_db
  if (file.exists(paste0("09_MagDatabaseProfiling/MAGsDatabaseMapping/", sample, "_host-filtered_flagstat.tsv"))) {
    mapped_full_db_host_filtered <- read.csv(paste0("09_MagDatabaseProfiling/MAGsDatabaseMapping/", sample, "_host-filtered_flagstat.tsv"), sep = "\t") %>%
                        filter(.[[3]] == "with itself and mate mapped") %>%
                          pull(1) %>%
                            as.integer()
  } else {
    mapped_full_db_host_filtered <- NA
  }
  values_to_bind <- c(sample, as.integer(c(raw_reads, trimmed, mapped_full_db, mapped_rep_db, host_filtered, mapped_full_db_host_filtered, mapped_host_db, unmapped_all)))
  df_reads <- rbind(df_reads, values_to_bind)
}
df_colnames <- c("Sample", "Raw", "Trimmed", "MAGs_DB", "Mapped_to_rep_MAGs", "Host_filtered_reads", "MAGs_DB_after_host_filtering", "Host_mapped", "Unmapped")
colnames(df_reads) <- df_colnames
df_reads <- df_reads %>%
              mutate(across(!c("Sample"), as.integer))
df_reads %>% filter(Sample %in% "A6.1")
df_reads_plot <- df_reads %>%
                          mutate(sum = Host_mapped + MAGs_DB) %>%
                          mutate(perc_mags = MAGs_DB/Trimmed*100) %>%
                          mutate(perc_mags_no_host = MAGs_DB/(Trimmed-Host_mapped)*100) %>%
                          mutate(perc_unmapped_no_host = 100 - perc_mags_no_host) %>%
                          mutate(perc_host = Host_mapped/Trimmed*100) %>%
                          mutate(perc_unmapped_all = 100-perc_mags-perc_host) %>%
                          mutate(perc_mapped = MAGs_DB/Trimmed*100) %>%
                          mutate(perc_unmapped = 100-perc_mapped) %>%
                          mutate(perc_mapped_hf = MAGs_DB_after_host_filtering/host_filtered*100) %>%
                          mutate(perc_unmapped_hf = 100-perc_mapped_hf) %>%
                          mutate(microbe_to_host_ratio = MAGs_DB/Host_mapped) %>%
                            pivot_longer(!Sample, names_to = "Type", values_to = "Number") %>%
                              mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                              mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) 
##############
# analyse data and plot
##############

Total_reads_MY <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Raw", "Trimmed")), aes(x=factor(Sample, levels = samples_MY), y=Number, fill = Type)) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~Host, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Total_reads_raw_trim_MY_samples.pdf")
Total_reads_IN <- ggplot(filter(df_reads_plot, Sample %in% samples_IN & Type %in% c("Raw", "Trimmed")), aes(x=factor(Sample, levels = samples_IN), y=Number, fill = Type)) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~Host, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Total_reads_raw_trim_IN_samples.pdf")
Total_reads_KE <- ggplot(filter(df_reads_plot, Sample %in% samples_KE & Type %in% c("Raw", "Trimmed")), aes(x=factor(Sample, levels = samples_KE), y=Number, fill = Type)) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~Host, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Total_reads_raw_trim_KE_samples.pdf")

Trimmed_reads_MY <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Trimmed")), aes(x=factor(Sample, levels = samples_MY), y=Number, fill = factor(Host, host_order))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_fill_manual(values=host_order_color)
                                  ggsave("Figures/01-Trimmed_reads_MY_samples.pdf")
Trimmed_reads_IN <- ggplot(filter(df_reads_plot, Sample %in% samples_IN & Type %in% c("Trimmed")), aes(x=factor(Sample, levels = samples_IN), y=Number, fill = factor(Host, host_order))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                scale_fill_manual(values=host_order_color)
                                  ggsave("Figures/01-Trimmed_reads_IN_samples.pdf")
Trimmed_reads_KE <- ggplot(filter(df_reads_plot, Sample %in% samples_KE & Type %in% c("Trimmed")), aes(x=factor(Sample, levels = samples_KE), y=Number, fill = factor(Host, host_order))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_fill_manual(values=host_order_color)
                                  ggsave("Figures/01-Trimmed_reads_KE_samples.pdf")

Mapped_to_MAGs_MY <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Trimmed", "MAGs_DB")), aes(x=factor(Sample, levels = samples_MY), y=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_MAGs_MY.pdf")
Mapped_to_MAGs_IN <- ggplot(filter(df_reads_plot, Sample %in% samples_IN & Type %in% c("Trimmed", "MAGs_DB")), aes(x=factor(Sample, levels = samples_IN), y=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_MAGs_IN.pdf")
Mapped_to_MAGs_KE <- ggplot(filter(df_reads_plot, Sample %in% samples_KE & Type %in% c("Trimmed", "MAGs_DB")), aes(x=factor(Sample, levels = samples_KE), y=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_MAGs_KE.pdf")

Mapped_to_host_MY <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Trimmed", "Host_mapped")), aes(x=factor(Sample, levels = samples_MY), y=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_host_MY.pdf")
Mapped_to_host_IN <- ggplot(filter(df_reads_plot, Sample %in% samples_IN & Type %in% c("Trimmed", "Host_mapped")), aes(x=factor(Sample, levels = samples_IN), y=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_host_IN.pdf")
Mapped_to_host_KE <- ggplot(filter(df_reads_plot, Sample %in% samples_KE & Type %in% c("Trimmed", "Host_mapped")), aes(x=factor(Sample, levels = samples_KE), y=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_host_KE.pdf")

Mapping_numbers_MY <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_MY),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_MY.pdf")
Mapping_numbers_IN <- ggplot(filter(df_reads_plot, Sample %in% samples_IN & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_IN),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_IN.pdf")
Mapping_numbers_KE <- ggplot(filter(df_reads_plot, Sample %in% samples_KE & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_KE),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_KE.pdf")

Mapping_numbers_Am <- ggplot(filter(df_reads_plot, Sample %in% samples_am & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_am),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ Location, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_Am.pdf")
Mapping_numbers_Ac <- ggplot(filter(df_reads_plot, Sample %in% samples_ac & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_ac),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ Location, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_Ac.pdf")
Mapping_numbers_Ad <- ggplot(filter(df_reads_plot, Sample %in% samples_ad & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_ad),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ Location, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_Ad.pdf")
Mapping_numbers_Af <- ggplot(filter(df_reads_plot, Sample %in% samples_af & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_af),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ Location, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_Af.pdf")
Mapping_numbers_Aa  <- ggplot(filter(df_reads_plot, Sample %in% samples_aa & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(x=factor(Sample, levels = samples_aa),
                                  y=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Type") +
                            facet_wrap(~ Location, scale = "free_x") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_Aa.pdf")
                                                             
Mapping_proportions_MY <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_MY),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_MY) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Host) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_MY.pdf")
Mapping_proportions_IN <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_IN),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_IN) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Host) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_IN.pdf")
Mapping_proportions_KE <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_KE & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_KE),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_KE) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Host) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Host, host_order), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_KE.pdf")

Mapping_proportions_Am <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_am & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_am),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_am) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Location) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Location), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_Am.pdf")
Mapping_proportions_Ac <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_ac & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_ac),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_ac) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Location) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Location), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_Ac.pdf")
Mapping_proportions_Ad <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_ad & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_ad),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_ad) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Location) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Location), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_Ad.pdf")
Mapping_proportions_Af <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_af & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_af),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_af) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Location) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Location), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_Af.pdf")
Mapping_proportions_Aa  <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_aa & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(x=factor(Sample, levels = samples_aa),
                                      y=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(x = "Sample",
                                  y = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_aa) %>%
                                                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                  mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                  group_by(Location) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(xmin=SampleID - 0.5,
                                              xmax=SampleID + 0.5,
                                              ymax = 108,
                                              ymin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                    facet_wrap(~ factor(Location), scale = "free_x") +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Proportion_Mapped_to_MAGs_and_host_Aa.pdf")

Total_reads <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Raw", "Trimmed")),
                                aes(y=factor(Sample, levels = samples_IN_MY),
                                    x=Number, fill = Type),
                                stat="identity", position = "dodge") +
                          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                          scale_fill_manual(values = brewer.pal(9, "Pastel1")[1:2]) +
                          ylab("Sample") +
                          xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         setFill = F, setCol = F,
                                         guide_nrow = 1,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Total_reads_raw_trim.pdf")
Trimmed_reads_conc <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Raw", "Trimmed")),
                                aes(y=factor(Sample, levels = samples_IN_MY),
                                    x=Number, fill = Type),
                                stat="identity", position = "dodge") +
                          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                          scale_fill_manual(values = brewer.pal(9, "Pastel1")[1:2]) +
                          ylab("Sample") +
                          xlab("Number of reads (paired end)") +
                          new_scale_fill() +
                          labs(fill = "DNA Concentration") +
                          scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                          geom_rect(data = filter(df_concentrations, ID %in% samples_IN_MY) %>%
                                            arrange(factor(ID, levels = samples_IN_MY)) %>%
                                            mutate(SampleID = row_number()), 
                                    aes(ymin=SampleID - 0.5,
                                        ymax=SampleID + 0.5,
                                        xmin = -5000000,
                                        xmax = -100000,
                                        fill = Concentration)
                                  ) +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         setFill = F, setCol = F,
                                         guide_nrow = 1,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Trimmed_reads_and_concentration.pdf")
Trimmed_reads_qPCR <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Raw", "Trimmed")),
                                aes(y=factor(Sample, levels = samples_IN_MY),
                                    x=Number, fill = Type),
                                stat="identity", position = "dodge") +
                          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                          scale_fill_manual(values = brewer.pal(9, "Pastel1")[1:2]) +
                          ylab("Sample") +
                          xlab("Number of reads (paired end)") +
                          new_scale_fill() +
                          labs(fill = "16S region Ct value") +
                          scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                          geom_rect(data = left_join(data.frame("Sample" = samples_IN_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                              arrange(factor(Sample, levels = samples_IN_MY)) %>%
                                               mutate(SampleID = row_number()), 
                                    aes(ymin=SampleID - 0.5,
                                        ymax=SampleID + 0.5,
                                        xmin = -5000000,
                                        xmax = -100000,
                                        fill = Ct_mean)
                                  ) +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         setFill = F, setCol = F,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Trimmed_reads_and_CT_values.pdf")
Trimmed_reads <- ggplot(filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Trimmed")),
                            aes(y=factor(Sample, levels = samples_IN_MY),
                                x=Number,
                                fill = factor(Host, host_order))) +
                        geom_bar(stat="identity", position = "dodge") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         y_size=7,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_fill_manual(values=host_order_color)
                                  ggsave("Figures/01-Trimmed_reads.pdf")
Mapped_to_MAGs <- ggplot(filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Trimmed", "MAGs_DB")), 
                          aes(y=factor(Sample, levels = samples_IN_MY), x=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         y_size=7,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_MAGs.pdf")
Mapped_to_host <- ggplot(filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Trimmed", "Host_mapped")), aes(y=factor(Sample, levels = samples_IN_MY), x=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         y_size=7,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_host.pdf")
Mapping_numbers <- ggplot(filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(y=factor(Sample, levels = samples_IN_MY),
                                  x=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB")))) +
                        geom_bar(stat="identity", position = "stack") +
                            labs(y = "Sample",
                                 x = "Number of reads (paired end)",
                                 fill = "Type") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         y_size=7,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             ))
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host.pdf")
Mapping_numbers_qPCR <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(y=factor(Sample, levels = samples_IN_MY),
                                  x=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                  stat="identity", position = "stack") +
                            labs(y = "Sample",
                                 x = "Number of reads (paired end)",
                                 fill = "Type") +
                        scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             )) +
                              new_scale_fill() +
                              labs(fill = "16S region Ct value") +
                              scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                              geom_rect(data = left_join(data.frame("Sample" = samples_IN_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                                  arrange(factor(Sample, levels = samples_IN_MY)) %>%
                                                  mutate(SampleID = row_number()), 
                                        aes(ymin=SampleID - 0.5,
                                            ymax=SampleID + 0.5,
                                            xmin = -5000000,
                                            xmax = -100000,
                                            fill = Ct_mean)
                                      ) +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         guide_nrow = 1, setFill = F,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_CT_values.pdf")
Mapping_numbers_Concentration <- ggplot() +
                                  geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                                          aes(y=factor(Sample, levels = samples_IN_MY),
                                              x=Number,
                                              fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                          stat="identity", position = "stack") +
                                      labs(y = "Sample",
                                          x = "Number of reads (paired end)",
                                          fill = "Type") +
                                      scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                      scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                                "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                                "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                                )) +
                                        new_scale_fill() +
                                        labs(fill = "DNA Concentration") +
                                        scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                                        geom_rect(data = filter(df_concentrations, ID %in% samples_IN_MY) %>%
                                                          arrange(factor(ID, levels = samples_IN_MY)) %>%
                                                          mutate(SampleID = row_number()), 
                                                  aes(ymin=SampleID - 0.5,
                                                      ymax=SampleID + 0.5,
                                                      xmin = -5000000,
                                                      xmax = -100000,
                                                      fill = Concentration)
                                                ) +
                                            make_theme(theme_name=theme_few(), leg_pos="right",
                                                      setFill = F, setCol = F,
                                                      guide_nrow = 1,
                                                      y_size=7, modify_guide = F,
                                                      x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                            ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration.pdf")
Mapping_proportions_reads  <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(y=factor(Sample, levels = samples_IN_MY),
                                      x=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(y = "Sample",
                                  x = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_IN_MY) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(ymin=SampleID - 0.5,
                                              ymax=SampleID + 0.5,
                                              xmax = 108,
                                              xmin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                y_size=7,
                                                x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                      ggsave("Figures/01-Mapped_to_MAGs_and_host.pdf")
Mapping_proportions_qPCR  <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(y=factor(Sample, levels = samples_IN_MY),
                                      x=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(y = "Sample",
                                  x = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_IN_MY) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(ymin=SampleID - 0.5,
                                              ymax=SampleID + 0.5,
                                              xmax = 108,
                                              xmin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                y_size=7,
                                                x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                      ggsave("Figures/01-Mapped_to_MAGs_and_host_CT_values.pdf")
Mapping_proportions_concentration  <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(y=factor(Sample, levels = samples_IN_MY),
                                      x=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(y = "Sample",
                                  x = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_IN_MY) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(ymin=SampleID - 0.5,
                                              ymax=SampleID + 0.5,
                                              xmax = 108,
                                              xmin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                y_size=7,
                                                x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                      ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration.pdf")



Total_reads <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Raw", "Trimmed")),
                                aes(y=factor(Sample, levels = samples_MY),
                                    x=Number, fill = Type),
                                stat="identity", position = "dodge") +
                          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                          scale_fill_manual(values = brewer.pal(9, "Pastel1")[1:2]) +
                          ylab("Sample") +
                          xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         setFill = F, setCol = F,
                                         guide_nrow = 1,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Total_reads_raw_trim_MY.pdf")
Trimmed_reads_conc <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Raw", "Trimmed")),
                                aes(y=factor(Sample, levels = samples_MY),
                                    x=Number, fill = Type),
                                stat="identity", position = "dodge") +
                          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                          scale_fill_manual(values = brewer.pal(9, "Pastel1")[1:2]) +
                          ylab("Sample") +
                          xlab("Number of reads (paired end)") +
                          new_scale_fill() +
                          labs(fill = "DNA Concentration") +
                          scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", trans = "log10") +
                          geom_rect(data = filter(df_concentrations, ID %in% samples_MY) %>%
                                            arrange(factor(ID, levels = samples_MY)) %>%
                                            mutate(SampleID = row_number()), 
                                    aes(ymin=SampleID - 0.5,
                                        ymax=SampleID + 0.5,
                                        xmin = -5000000,
                                        xmax = -100000,
                                        fill = Concentration)
                                  ) +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         setFill = F, setCol = F,
                                         guide_nrow = 1,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Trimmed_reads_and_concentration_MY.pdf")
Trimmed_reads_qPCR <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Raw", "Trimmed")),
                                aes(y=factor(Sample, levels = samples_MY),
                                    x=Number, fill = Type),
                                stat="identity", position = "dodge") +
                          scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                          scale_fill_manual(values = brewer.pal(9, "Pastel1")[1:2]) +
                          ylab("Sample") +
                          xlab("Number of reads (paired end)") +
                          new_scale_fill() +
                          labs(fill = "16S region Ct value") +
                          scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Greys"))) +
                          geom_rect(data = left_join(data.frame("Sample" = samples_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                              arrange(factor(Sample, levels = samples_MY)) %>%
                                               mutate(SampleID = row_number()), 
                                    aes(ymin=SampleID - 0.5,
                                        ymax=SampleID + 0.5,
                                        xmin = -5000000,
                                        xmax = -100000,
                                        fill = Ct_mean)
                                  ) +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         setFill = F, setCol = F,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Trimmed_reads_and_CT_values_MY.pdf")
Trimmed_reads <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Trimmed")),
                            aes(y=factor(Sample, levels = samples_MY),
                                x=Number,
                                fill = factor(Host, host_order))) +
                        geom_bar(stat="identity", position = "dodge") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         y_size=7,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_fill_manual(values=host_order_color)
                                  ggsave("Figures/01-Trimmed_reads_MY.pdf")
Mapped_to_MAGs <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Trimmed", "MAGs_DB")), 
                          aes(y=factor(Sample, levels = samples_MY), x=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         y_size=7,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_MAGs_MY.pdf")
Mapped_to_host <- ggplot(filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Trimmed", "Host_mapped")), aes(y=factor(Sample, levels = samples_MY), x=Number, fill = factor(Type))) +
                        geom_bar(stat="identity", position = "dodge") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         y_size=7,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Mapped_to_host_MY.pdf")
Mapping_numbers_qPCR <- ggplot() +
                        geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                              aes(y=factor(Sample, levels = samples_MY),
                                  x=Number,
                                  fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                  stat="identity", position = "stack") +
                            labs(y = "Sample",
                                 x = "Number of reads (paired end)",
                                 fill = "Type") +
                        scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                             )) +
                              new_scale_fill() +
                              labs(fill = "16S region Ct value") +
                              scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                              geom_rect(data = left_join(data.frame("Sample" = samples_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                                  arrange(factor(Sample, levels = samples_MY)) %>%
                                                  mutate(SampleID = row_number()), 
                                        aes(ymin=SampleID - 0.5,
                                            ymax=SampleID + 0.5,
                                            xmin = -5000000,
                                            xmax = -100000,
                                            fill = Ct_mean)
                                      ) +
                              make_theme(theme_name=theme_few(), leg_pos="right",
                                         guide_nrow = 1, setFill = F,
                                         y_size=7, modify_guide = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                  ggsave("Figures/01-Mapped_to_MAGs_and_host_CT_values_MY.pdf")
Mapping_numbers_Concentration <- ggplot() +
                                  geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                                          aes(y=factor(Sample, levels = samples_MY),
                                              x=Number,
                                              fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                          stat="identity", position = "stack") +
                                      labs(y = "Sample",
                                          x = "Number of reads (paired end)",
                                          fill = "Type") +
                                      scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                      scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                                "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                                "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                                )) +
                                        new_scale_fill() +
                                        labs(fill = "DNA Concentration") +
                                        scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                                        geom_rect(data = filter(df_concentrations, ID %in% samples_MY) %>%
                                                          arrange(factor(ID, levels = samples_MY)) %>%
                                                          mutate(SampleID = row_number()), 
                                                  aes(ymin=SampleID - 0.5,
                                                      ymax=SampleID + 0.5,
                                                      xmin = -5000000,
                                                      xmax = -100000,
                                                      fill = Concentration)
                                                ) +
                                            make_theme(theme_name=theme_few(), leg_pos="right",
                                                      setFill = F, setCol = F,
                                                      guide_nrow = 1,
                                                      y_size=7, modify_guide = F,
                                                      x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                            ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration_MY.pdf")
Mapping_numbers_Concentration_qPCR <- ggplot() +
                                  geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                                          aes(y=factor(Sample, levels = samples_MY),
                                              x=Number,
                                              fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                          stat="identity", position = "stack") +
                                      labs(y = "Sample",
                                          x = "Number of reads (paired end)",
                                          fill = "Type") +
                                      scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                      scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                                "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                                "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                                )) +
                                        new_scale_fill() +
                                        labs(fill = "DNA Concentration") +
                                        scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                                        geom_rect(data = filter(df_concentrations, ID %in% samples_MY) %>%
                                                          arrange(factor(ID, levels = samples_MY)) %>%
                                                          mutate(SampleID = row_number()), 
                                                  aes(ymin=SampleID - 0.5,
                                                      ymax=SampleID + 0.5,
                                                      xmin = -5000000,
                                                      xmax = -100000,
                                                      fill = Concentration)
                                                ) +
                                        new_scale_fill() +
                                        labs(fill = "16S region Ct value") +
                                        scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                                        geom_rect(data = left_join(data.frame("Sample" = samples_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                                            arrange(factor(Sample, levels = samples_MY)) %>%
                                                            mutate(SampleID = row_number()), 
                                                  aes(ymin=SampleID - 0.5,
                                                      ymax=SampleID + 0.5,
                                                      xmin = -10000000,
                                                      xmax = -5001000,
                                                      fill = Ct_mean)
                                                ) +
                                            make_theme(theme_name=theme_few(), leg_pos="right",
                                                      setFill = F, setCol = F,
                                                      guide_nrow = 1,
                                                      y_size=7, modify_guide = F,
                                                      x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                            ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration_qPCR_MY.pdf")
Mapping_proportions_reads  <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(y=factor(Sample, levels = samples_MY),
                                      x=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(y = "Sample",
                                  x = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_MY) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(ymin=SampleID - 0.5,
                                              ymax=SampleID + 0.5,
                                              xmax = 108,
                                              xmin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Greys"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                y_size=7,
                                                x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                      ggsave("Figures/01-Mapped_to_MAGs_and_host_MY.pdf")
Mapping_proportions_qPCR  <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(y=factor(Sample, levels = samples_MY),
                                      x=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(y = "Sample",
                                  x = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_MY) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(ymin=SampleID - 0.5,
                                              ymax=SampleID + 0.5,
                                              xmax = 108,
                                              xmin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Blues"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                y_size=7,
                                                x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                      ggsave("Figures/01-Mapped_to_MAGs_and_host_CT_values_MY.pdf")
Mapping_proportions_concentration  <- ggplot() +
                            geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                  aes(y=factor(Sample, levels = samples_MY),
                                      x=Number,
                                      fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                      stat="identity", position = "stack") +
                            labs(y = "Sample",
                                  x = "Number of reads (paired end)",
                                  fill = "Type") +
                              scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                        "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                        "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                        )) +
                                new_scale_fill() +
                                geom_rect(data = filter(df_reads, Sample %in% samples_MY) %>%
                                                    mutate(SampleID = row_number()), 
                                          aes(ymin=SampleID - 0.5,
                                              ymax=SampleID + 0.5,
                                              xmax = 108,
                                              xmin = 102,
                                              fill = Trimmed)
                                        ) +
                                    labs(fill = "#Trimmed reads") +
                                    scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", labels=unit_format(unit = "M", scale = 1e-6)) +
                                      make_theme(setFill = F, setCol = F, theme_name=theme_few(), leg_pos="right",
                                                guide_nrow = 1, modify_guide = F,
                                                y_size=7,
                                                x_angle=0 ,x_vj=0, x_hj=0, x_size=12)
                                      ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration_MY.pdf")

Mapping_numbers_Concentration_qPCR <- ggplot() +
                                  geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                                          aes(x=factor(Sample, levels = samples_MY),
                                              y=Number,
                                              fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                          stat="identity", position = "stack") +
                                      labs(x = "Sample",
                                          y = "Number of reads (paired end)",
                                          fill = "Type") +
                                      scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                      scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                                "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                                "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                                )) +
                                        new_scale_fill() +
                                        labs(fill = "DNA Concentration") +
                                        scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                                        geom_rect(data = filter(df_concentrations, ID %in% samples_MY) %>%
                                                          arrange(factor(ID, levels = samples_MY)) %>%
                                                          mutate(Host = Vectorize(get_host_from_sample_name)(ID)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()), 
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -5000000,
                                                      ymax = -100000,
                                                      fill = Concentration)
                                                ) +
                                        new_scale_fill() +
                                        labs(fill = "16S region Ct value") +
                                        scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                                        geom_rect(data = left_join(data.frame("Sample" = samples_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                                            arrange(factor(Sample, levels = samples_MY)) %>%
                                                            mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()),
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -10000000,
                                                      ymax = -5001000,
                                                      fill = Ct_mean)
                                                ) +
                                            facet_wrap(~ factor(Host, levels = host_order), scale = "free_x") +
                                            make_theme(theme_name=theme_few(), leg_pos="right",
                                                      setFill = F, setCol = F,
                                                      guide_nrow = 1,
                                                      modify_guide = F,
                                                      x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                      y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                            ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration_qPCR_by_host_MY.pdf")

Mapping_proportion_Concentration_qPCR <- ggplot() +
                                  geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("perc_host", "perc_mags", "perc_unmapped_all")),
                                          aes(x=factor(Sample, levels = samples_MY),
                                              y=Number,
                                              fill = factor(Type, levels = c("perc_unmapped_all", "perc_host", "perc_mags"))),
                                          stat="identity", position = "stack") +
                                      labs(x = "Sample",
                                          y = "Number of reads (paired end)",
                                          fill = "Type") +
                                      scale_y_continuous(labels=unit_format(unit = "%", scale = 1e-0)) +
                                      scale_fill_manual(values=c("perc_unmapped_all" = brewer.pal(9, "Pastel1")[1],
                                                                "perc_mags" = brewer.pal(9, "Pastel1")[3],
                                                                "perc_host" = brewer.pal(9, "Pastel1")[2]
                                                                )) +
                                        new_scale_fill() +
                                        labs(fill = "DNA Concentration") +
                                        scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                                        geom_rect(data = filter(df_concentrations, ID %in% samples_MY) %>%
                                                          arrange(factor(ID, levels = samples_MY)) %>%
                                                          mutate(Host = Vectorize(get_host_from_sample_name)(ID)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()), 
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = 0,
                                                      ymax = -5,
                                                      fill = Concentration)
                                                ) +
                                        new_scale_fill() +
                                        labs(fill = "16S region Ct value") +
                                        scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                                        geom_rect(data = left_join(data.frame("Sample" = samples_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                                            arrange(factor(Sample, levels = samples_MY)) %>%
                                                            mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()),
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -5,
                                                      ymax = -10,
                                                      fill = Ct_mean)
                                                ) +
                                            facet_wrap(~ factor(Host, levels = host_order), scale = "free_x") +
                                            make_theme(theme_name=theme_few(), leg_pos="right",
                                                      setFill = F, setCol = F,
                                                      guide_nrow = 1,
                                                      modify_guide = F,
                                                      x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                                      y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                            ggsave("Figures/01-Proportion_mapped_to_MAGs_and_host_concentration_qPCR_by_host_MY.pdf")


ggplot(df_reads %>% filter((Sample %in% samples)) %>% 
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample))
        ) +
    geom_point(aes(x = factor(Sample, samples),
                   y = MAGs_DB,
                   shape = Location,
                   color = factor(Host, host_order))) +
    labs(x = "Sample", y = "Mapped to microbial MAGs DB", color = "Host_mapped", shape = "Location") +
    geom_hline(yintercept = 1e+06) +
    scale_color_manual(values=host_order_color_dark) +
    scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6), limits = c(0, 1e+08)) +
    make_theme(setCol = F, guide_nrow = 5, x_size = 0)
    ggsave("Figures/01-Total_reads_MAG_mapped_reads_per_samples.pdf")

ggplot(df_reads %>% filter((Sample %in% samples)) %>% 
        mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(Sample))
        ) +
    geom_point(aes(x = factor(Sample, samples),
                   y = Trimmed,
                   shape = Location,
                   color = factor(Host, host_order))) +
    labs(x = "Sample", y = "Total trimmed reads", color = "Host_mapped", shape = "Location") +
    geom_hline(yintercept = 1e+06) +
    scale_color_manual(values=host_order_color_dark) +
    scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6), limits = c(0, 1e+08)) +
    make_theme(setCol = F, guide_nrow = 5, x_size = 0)
      ggsave("Figures/01-Total_reads_trimmed_reads_per_samples.pdf")

Mapping_numbers_Concentration_qPCR_all <- ggplot() +
                                  geom_bar(data = filter(df_reads_plot, Sample %in% samples & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                                          aes(x=factor(Sample, levels = samples),
                                              y=Number,
                                              fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                          stat="identity", position = "stack") +
                                      labs(x = "Sample",
                                          y = "Number of reads (paired end)",
                                          fill = "Type") +
                                      scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                      scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                                "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                                "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                                )) +
                                        new_scale_fill() +
                                        labs(fill = "DNA Concentration") +
                                        scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                                        geom_rect(data = filter(df_concentrations, ID %in% samples) %>%
                                                          arrange(factor(ID, levels = samples)) %>%
                                                          mutate(Host = Vectorize(get_host_from_sample_name)(ID)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()), 
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -10000000,
                                                      ymax = -5001000,
                                                      fill = Concentration)
                                                ) +
                                        new_scale_fill() +
                                        labs(fill = "16S region Ct value") +
                                        scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                                        geom_rect(data = left_join(data.frame("Sample" = samples), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                                            arrange(factor(Sample, levels = samples)) %>%
                                                            mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()),
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -15000000,
                                                      ymax = -10001000,
                                                      fill = Ct_mean)
                                                ) +
                                        new_scale_fill() +
                                        labs(fill = "Location") +
                                        scale_fill_brewer(palette="Set1") +
                                        geom_rect(data = data.frame("Sample" = samples) %>%
                                                            arrange(factor(Sample, levels = samples)) %>%
                                                            mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                            mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()),
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -5000000,
                                                      ymax = -100000,
                                                      fill = Location)
                                                ) +
                                            facet_wrap(~ factor(Host, levels = host_order), scale = "free_x") +
                                            geom_hline(yintercept = 1e+06, linetype = "solid") +
                                            annotate("text", y=2e+06, x=8, label="1 M", size = 2) +
                                            geom_hline(yintercept = 5e+06, linetype = "dashed") +
                                            annotate("text", y=6e+06, x=8, label="5 M", size = 2) +
                                            geom_hline(yintercept = 1e+07, linetype = "dotted") +
                                            annotate("text", y=1.1e+07, x=8, label="10 M", size = 2) +
                                            make_theme(theme_name=theme_few(), leg_pos="right",
                                                      setFill = F, setCol = F,
                                                      guide_nrow = 1,
                                                      modify_guide = F,
                                                      x_angle=45 ,x_vj=1.2, x_hj=1, x_size=0,
                                                      y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                            ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration_qPCR_by_host.pdf")

Mapping_numbers_Concentration_qPCR_MY_IN <- ggplot() +
                                  geom_bar(data = filter(df_reads_plot, Sample %in% samples_IN_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
                                          aes(x=factor(Sample, levels = samples_IN_MY),
                                              y=Number,
                                              fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))),
                                          stat="identity", position = "stack") +
                                      labs(x = "Sample",
                                          y = "Number of reads (paired end)",
                                          fill = "Type") +
                                      scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                      scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                                                                "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                                                                "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                                                                )) +
                                        new_scale_fill() +
                                        labs(fill = "DNA Concentration") +
                                        scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "Reds"), guide = "colourbar", trans = "log10") +
                                        geom_rect(data = filter(df_concentrations, ID %in% samples_IN_MY) %>%
                                                          arrange(factor(ID, levels = samples_IN_MY)) %>%
                                                          mutate(Host = Vectorize(get_host_from_sample_name)(ID)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()), 
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -10000000,
                                                      ymax = -5001000,
                                                      fill = Concentration)
                                                ) +
                                        new_scale_fill() +
                                        labs(fill = "16S region Ct value") +
                                        scale_fill_gradientn(na.value = "transparent", colors = rev(brewer.pal(4, "Blues"))) +
                                        geom_rect(data = left_join(data.frame("Sample" = samples_IN_MY), qpcr_df_info, by = c("Sample" = "Sample.Name")) %>%
                                                            arrange(factor(Sample, levels = samples_IN_MY)) %>%
                                                            mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()),
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -15000000,
                                                      ymax = -10001000,
                                                      fill = Ct_mean)
                                                ) +
                                        new_scale_fill() +
                                        labs(fill = "Location") +
                                        scale_fill_brewer(palette="Set1") +
                                        geom_rect(data = data.frame("Sample" = samples_IN_MY) %>%
                                                            arrange(factor(Sample, levels = samples_IN_MY)) %>%
                                                            mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                                                            mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                                                            group_by(Host) %>%
                                                              mutate(SampleID = row_number()),
                                                  aes(xmin=SampleID - 0.5,
                                                      xmax=SampleID + 0.5,
                                                      ymin = -5000000,
                                                      ymax = -100000,
                                                      fill = Location)
                                                ) +
                                            facet_wrap(~ factor(Host, levels = host_order), scale = "free_x") +
                                            geom_hline(yintercept = 1e+06, linetype = "solid") +
                                            annotate("text", y=2e+06, x=8, label="1 M", size = 2) +
                                            geom_hline(yintercept = 5e+06, linetype = "dashed") +
                                            annotate("text", y=6e+06, x=8, label="5 M", size = 2) +
                                            geom_hline(yintercept = 1e+07, linetype = "dotted") +
                                            annotate("text", y=1.1e+07, x=8, label="10 M", size = 2) +
                                            make_theme(theme_name=theme_few(), leg_pos="right",
                                                      setFill = F, setCol = F,
                                                      guide_nrow = 1,
                                                      modify_guide = F,
                                                      x_angle=45 ,x_vj=1.2, x_hj=1, x_size=0,
                                                      y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                            ggsave("Figures/01-Mapped_to_MAGs_and_host_concentration_qPCR_IN_MY_by_host.pdf")
# ratio of host to microbe DNA
ggplot() +
  geom_bar(data = df_reads_plot %>%
                    filter(Type %in% c("microbe_to_host_ratio") & Sample %in% samples_MY) %>%
                    left_join(df_reads_plot %>% 
                              filter(Type %in% c("Trimmed") & Sample %in% samples_MY) %>%
                                summarise(Sample, Trimmed = Number)
                              ),
           aes(x = Sample,
               y = Number,
               fill = Trimmed),
               stat="identity") +
  scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "RdYlGn"), guide = "legend", labels=unit_format(unit = "M", scale = 1e-6), trans = "log10") +
  labs(x = "Sample", y = "Microbe to Host DNA Ratio", fill = "#Trimmed") +
  new_scale_fill() +
  geom_rect(data = filter(df_reads, Sample %in% samples_MY) %>%
                    mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                    mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                    group_by(Host) %>%
                      mutate(SampleID = row_number()) %>%
                        left_join(df_reads_plot %>% 
                                                              filter(Type %in% c("microbe_to_host_ratio") & Sample %in% samples_MY) %>%
                                                                group_by(Host) %>%
                                                                mutate(ratio_max = max(Number))
                                                              ), 
            aes(xmin=SampleID - 0.5,
                xmax=SampleID + 0.5,
                ymin = 0,
                ymax = -10,
                # ymin = -ratio_max/200,
                # ymax = -ratio_max/25,
                fill = MAGs_DB)
          ) +
  labs(fill = "#Mapped to MAGs DB") +
  scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "RdYlGn"), guide = "legend", labels=unit_format(unit = "M", scale = 1e-6)) +
  facet_wrap(~factor(Host, host_order), scale = "free_x") +
  make_theme(theme_name=theme_few(), leg_pos="right",
             setFill = F, setCol = F,
             guide_nrow = 1,
             modify_guide = F,
             x_angle=45 ,x_vj=1.2, x_hj=1, x_size=0,
             y_angle=0 ,y_vj=0, y_hj=0, y_size=12)

ggplot() +
  geom_bar(data = df_reads_plot %>%
                    filter(Type %in% c("microbe_to_host_ratio") & Sample %in% samples_MY & Host != "Apis mellifera") %>%
                    left_join(df_reads_plot %>% 
                              filter(Type %in% c("Trimmed") & Sample %in% samples_MY) %>%
                                summarise(Sample, Trimmed = Number)
                              ),
           aes(x = Sample,
               y = Number,
               fill = Trimmed),
               stat="identity") +
  scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "RdYlGn"), guide = "legend", labels=unit_format(unit = "M", scale = 1e-6), trans = "log10") +
  labs(x = "Sample", y = "Microbe to Host DNA Ratio", fill = "#Trimmed") +
  new_scale_fill() +
  geom_rect(data = filter(df_reads, Sample %in% samples_MY) %>%
                    mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                    mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                    group_by(Host) %>%
                    filter(Host != "Apis mellifera") %>%
                      mutate(SampleID = row_number()) %>%
                        left_join(df_reads_plot %>% 
                                                              filter(Type %in% c("microbe_to_host_ratio") & Sample %in% samples_MY) %>%
                                                                group_by(Host) %>%
                                                                mutate(ratio_max = max(Number))
                                                              ), 
            aes(xmin=SampleID - 0.5,
                xmax=SampleID + 0.5,
                ymin = 0,
                ymax = -10,
                # ymin = -ratio_max/200,
                # ymax = -ratio_max/25,
                fill = MAGs_DB)
          ) +
  labs(fill = "#Mapped to MAGs DB") +
  scale_fill_gradientn(na.value = "transparent", colors = brewer.pal(4, "RdYlGn"), guide = "legend", labels=unit_format(unit = "M", scale = 1e-6)) +
  facet_wrap(~factor(Host, host_order), scale = "free_x") +
  make_theme(theme_name=theme_few(), leg_pos="right",
             setFill = F, setCol = F,
             guide_nrow = 1,
             modify_guide = F,
             x_angle=45 ,x_vj=1.2, x_hj=1, x_size=0,
             y_angle=0 ,y_vj=0, y_hj=0, y_size=12)

df_reads_updated <-  df_reads %>%
  select(Sample, MAGs_DB, Host_mapped, Unmapped, Trimmed) %>%
    mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
    mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
      filter(Location %in% c("Malaysia")) %>%
        mutate(MAGs_DB_range = ifelse(MAGs_DB >= 30e+06, "High", NA)) %>%
        mutate(MAGs_DB_range = ifelse(MAGs_DB >= 15e+06 & MAGs_DB < 30e+06, "Sufficient", MAGs_DB_range)) %>%
        mutate(MAGs_DB_range = ifelse(MAGs_DB >= 0 & MAGs_DB < 15e+06, "Low", MAGs_DB_range)) %>%
          mutate(Total_depth = ifelse(Trimmed >= 50e+06, "Target_reached", "Insufficient"))

ggplot() +
  geom_bar(data = filter(df_reads_plot, Sample %in% samples_MY & Type %in% c("Host_mapped", "MAGs_DB", "Unmapped")),
           aes(x=factor(Sample, levels = samples_MY),
               y=Number,
               fill = factor(Type, levels = c("Unmapped", "Host_mapped", "MAGs_DB"))
              ),
           stat="identity", position = "stack") +
  labs(x = "Sample",
       y = "Number of reads (paired end)",
       fill = "Type") +
  scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
  scale_fill_manual(values=c("Unmapped" = brewer.pal(9, "Pastel1")[1],
                             "MAGs_DB" = brewer.pal(9, "Pastel1")[3],
                             "Host_mapped" = brewer.pal(9, "Pastel1")[2]
                            ),
                    labels=c("Unmapped",
                             "Mapped to MAGs",
                             "Mapped to Host"
                            )
                   ) +     
  new_scale_fill() +
  geom_rect(data = filter(df_reads_updated, Sample %in% samples_MY) %>%
                    mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                    mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                    group_by(Host) %>%
                      mutate(SampleID = row_number()), 
            aes(xmin=SampleID - 0.5,
                xmax=SampleID + 0.5,
                ymin = -100000,
                ymax = -5000000,
                fill = factor(MAGs_DB_range, levels = c("Low", "Sufficient", "High")))
          ) +
  labs(fill = "#Mapped to MAGs DB") +
  scale_fill_manual(values=c("Low" = brewer.pal(9, "Paired")[6],
                             "Sufficient" = brewer.pal(9, "Paired")[3],
                             "High" = brewer.pal(9, "Paired")[4]
                            ),
                    labels=c("Low (< 15M)",
                             "Sufficient (15M - 30M)",
                             "High (> 30M)"
                            ),
                   ) +
  new_scale_fill() +
  geom_rect(data = filter(df_reads_updated, Sample %in% samples_MY) %>%
                    mutate(Host = Vectorize(get_host_from_sample_name)(Sample)) %>%
                    mutate(Location = Vectorize(get_location_from_sample_name)(Sample)) %>%
                    group_by(Host) %>%
                      mutate(SampleID = row_number()), 
            aes(xmin=SampleID - 0.5,
                xmax=SampleID + 0.5,
                ymin = -5000000,
                ymax = -10000000,
                fill = factor(Total_depth, levels = c("Target_reached", "Insufficient")))
          ) +
  labs(fill = "#Mapped to MAGs DB") +
  scale_fill_manual(values=c("Target_reached" = brewer.pal(9, "Paired")[4],
                             "Insufficient" = brewer.pal(9, "Paired")[6]
                            ),
                    labels=c("Target reached (>50M)",
                             "Insufficient"
                            ),
                   ) +
  facet_wrap(~factor(Host, host_order), scale = "free_x") +
  make_theme(theme_name=theme_few(), leg_pos="right",
             setFill = F, setCol = F,
             guide_nrow = 1,
             modify_guide = F,
             x_angle=45 ,x_vj=1.2, x_hj=1, x_size=0,
             y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
             ggsave("Figures/01-Sequecning_Depth_target_and_mapped_reads.pdf")
