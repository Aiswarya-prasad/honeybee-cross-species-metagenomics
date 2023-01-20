source('scripts/visualization/utilities.R', chdir = TRUE)

############
# read info
############

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
  raw_reads <- number_raw_R1 + number_raw_R2
  trimmed <- number_trimmed_R1 + number_trimmed_R2
  if (file.exists(paste0("02_HostMapping/", sample, "_flagstat.tsv"))) {
    mapped_host_db <- read.csv(paste0("02_HostMapping/", sample, "_flagstat.tsv"), sep = "\t") %>%
                      filter(.[[3]] == "with itself and mate mapped") %>%
                        pull(1) %>%
                          as.integer()
  } else {
   mapped_host_db <- NA 
  }
  unmapped_host_db <- trimmed - mapped_host_db
  if (file.exists(paste0("04_MicrobiomeMappingDirect/", sample, "_flagstat.tsv"))) {
    mapped_mic_db <- read.csv(paste0("04_MicrobiomeMappingDirect/", sample, "_flagstat.tsv"), sep = "\t") %>%
                      filter(.[[3]] == "with itself and mate mapped") %>%
                        pull(1) %>%
                          as.integer()
  } else {
   mapped_mic_db <- NA 
  }
  unmapped_mic_db <- trimmed - mapped_mic_db
  if (file.exists(paste0("03_MicrobiomeMapping/", sample, "_flagstat.tsv"))) {
    mapped_mic_db_host_filtered <- read.csv(paste0("03_MicrobiomeMapping/", sample, "_flagstat.tsv"), sep = "\t") %>%
                      filter(.[[3]] == "with itself and mate mapped") %>%
                        pull(1) %>%
                          as.integer()
  } else {
   mapped_mic_db_host_filtered <- NA 
  }
  unmapped_sequential <- unmapped_host_db - mapped_mic_db_host_filtered
  values_to_bind <- c(sample, as.integer(c(raw_reads, trimmed, mapped_host_db, unmapped_host_db, mapped_mic_db, unmapped_mic_db, mapped_mic_db_host_filtered, unmapped_sequential)))
  df_reads <- rbind(df_reads, values_to_bind)
}
# parse this later
# read.csv(paste0("02_HostMapping/", sample, "_coverage.tsv"), sep = "\t")
df_colnames <- c("Sample", "Raw", "Trimmed", "Mapped_host", "Unmapped_host", "Mapped_microbiome", "Unmapped_microbiome", "Mapped_microbiome_filtered", "Unmapped_filtered")
colnames(df_reads) <- df_colnames
df_reads <- df_reads %>%
              mutate(across(!c("Sample"), as.integer))
df_meta <- read.csv("config/Metadata_211018_Medgenome_india_samples.csv", sep = ',')
colnames(df_meta)[which(colnames(df_meta) == "ID")] <- "Sample"
df_meta$SpeciesID <- recode(df_meta$SpeciesID, "Am" = "Apis mellifera", "Ac" = "Apis cerana", "Af" = "Apis florea", "Ad" = "Apis dorsata")
df_meta <- df_meta %>%
            filter(Sample %in% samples) %>%
              arrange(match(Sample, samples))
df_meta %>% group_by(SpeciesID) %>% tally()
df_meta %>% group_by(SpeciesID, Country) %>% tally()
df_meta %>% group_by(SpeciesID, Country, Colony) %>% tally()
df_plot_reads <- pivot_longer(df_reads, !Sample, values_to = "Number", names_to = "Type") %>%
                  merge(select(df_meta, Sample, SpeciesID), by="Sample")
df_plot_reads$Number <- as.integer(df_plot_reads$Number)
df_meta_complete

Total_reads <- ggplot(filter(df_plot_reads, Type %in% c("Raw", "Trimmed")), aes(x=factor(Sample, levels = samples), y=Number, fill = Type)) +
                        geom_bar(stat="identity", position = "dodge") +
                          # ggtitle("Total reads per sample") +
                            xlab("Sample") +
                            ylab("Number of reads (paired end)") +
                            geom_vline(xintercept = 5.5, linetype="solid") +
                            geom_vline(xintercept = 20.5, linetype="solid") +
                            geom_vline(xintercept = 35.5, linetype="solid") +
                            geom_hline(yintercept = 30e+6, linetype="dotted") +
                            geom_hline(yintercept = 50e+6, linetype="dotted") +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6))
                                  ggsave("Figures/01-Total_reads_trimming.pdf")

Total_reads_species <- ggplot(filter(df_plot_reads, Type %in% c("Trimmed")), aes(x=factor(Sample, levels = samples), y=Number, fill = SpeciesID)) +
                        geom_bar(stat="identity", position = "dodge") +
                          # ggtitle("Total reads per sample") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                            geom_vline(xintercept = 5.5, linetype="solid") +
                            geom_vline(xintercept = 20.5, linetype="solid") +
                            geom_vline(xintercept = 35.5, linetype="solid") +
                            geom_hline(yintercept = 30e+6, linetype="dotted") +
                            make_theme(theme_name=theme_few(), leg_pos="bottom",
                                       guide_nrow = 1,
                                       x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                       y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_fill_manual(values=host_order_color)
                                  ggsave("Figures/01-Total_reads_species_after_trimming.pdf")

df_plot_reads_prop <- select(df_reads, Sample, Trimmed, Mapped_host, Mapped_microbiome_filtered, Unmapped_filtered) %>%
  summarise(Sample, "None" = Unmapped_filtered/Trimmed*100,
            "Host" = Mapped_host/Trimmed*100,
            "Microbiome" = Mapped_microbiome_filtered/Trimmed*100) %>%
              pivot_longer(!Sample, names_to = "Type", values_to = "Number")
Mapped_Unmapped_reads_prop <- ggplot(df_plot_reads_prop, aes(x=factor(Sample, levels = samples),
                                                            y=Number,
                                                            fill=factor(Type, levels = c("None", "Host", "Microbiome")))) +
                        geom_bar(stat="identity", position = "stack") +
                          # ggtitle("Total reads per sample") +
                            labs(y = "Sample",
                                 x = "Percentage of reads (paired end)",
                                 fill = "Mapped to") +
                            geom_vline(xintercept = 5.5, linetype="solid") +
                            geom_vline(xintercept = 20.5, linetype="solid") +
                            geom_vline(xintercept = 35.5, linetype="solid") +
                            make_theme(theme_name=theme_few(), leg_pos="bottom",
                                       guide_nrow = 1,
                                       x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                       y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Mapped_Unmapped_reads_prop.pdf")

df_plot_reads_num <- select(df_reads, Sample, Trimmed, Mapped_host, Mapped_microbiome_filtered, Unmapped_filtered) %>%
  summarise(Sample, "None" = Unmapped_filtered,
            "Host" = Mapped_host,
            "Microbiome" = Mapped_microbiome_filtered) %>%
              pivot_longer(!Sample, names_to = "Type", values_to = "Number")
Mapped_Unmapped_reads_num <- ggplot(df_plot_reads_num, aes(x=factor(Sample, levels = samples),
                                                            y=Number,
                                                            fill=factor(Type, levels = c("None", "Host", "Microbiome")))) +
                        geom_bar(stat="identity", position = "stack") +
                          # ggtitle("Total reads per sample") +
                            labs(x = "Sample",
                                 y = "Percentage of reads (paired end)",
                                 fill = "Mapped to") +
                            geom_vline(xintercept = 5.5, linetype="solid") +
                            geom_vline(xintercept = 20.5, linetype="solid") +
                            geom_vline(xintercept = 35.5, linetype="solid") +
                            geom_hline(yintercept = 10e+6, linetype="dotted") +
                            geom_hline(yintercept = 1e+6, linetype="dotted") +
                            scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                            make_theme(theme_name=theme_few(), leg_pos="bottom",
                                       guide_nrow = 1,
                                       x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                       y_angle=0 ,y_vj=0, y_hj=0, y_size=12)
                                      ggsave("Figures/01-Mapped_Unmapped_reads_num.pdf")

df_plot_reads_non_specific <- df_reads %>%
  mutate(Unmapped_microbiome_filtered = Mapped_host + Unmapped_filtered) %>%
    select(Sample, Mapped_microbiome, Unmapped_microbiome, Mapped_microbiome_filtered, Unmapped_microbiome_filtered) %>%
              pivot_longer(!Sample, names_to = "Type", values_to = "Number") %>%
                mutate(Host_filtered = ifelse(Type %in% c("Mapped_microbiome", "Unmapped_microbiome"), "Direct mapping", "Host filtered")) %>%
                  mutate(Type = ifelse(Type %in% c("Mapped_microbiome", "Mapped_microbiome_filtered"), "Mapped", "Unmapped"))


Non_specific_reads <- ggplot(df_plot_reads_non_specific %>% filter(Type == "Mapped"),
                             aes(x=factor(Sample, samples),
                                 y=Number,
                                 fill=Host_filtered)) +
                        geom_bar(stat="identity", position = "dodge") +
                            labs(x = "Sample",
                                 y = "Number of reads (paired end)",
                                 fill = "Mapping strategy") +
                            geom_vline(xintercept = 5.5, linetype="solid") +
                            geom_vline(xintercept = 20.5, linetype="solid") +
                            geom_vline(xintercept = 35.5, linetype="solid") +
                            geom_hline(yintercept = 10e+6, linetype="dotted") +
                            geom_hline(yintercept = 1e+6, linetype="dotted") +
                              scale_y_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                              make_theme(theme_name=theme_few(), leg_pos="bottom",
                                         guide_nrow = 1, setFill = F,
                                         x_angle=45 ,x_vj=1.2, x_hj=1, x_size=7,
                                         y_angle=0 ,y_vj=0, y_hj=0, y_size=12) +
                                scale_fill_manual(values=brewer.pal(4, "Paired")[c(3,4)])
                                      ggsave("Figures/01-Mapping_non_specific.pdf")