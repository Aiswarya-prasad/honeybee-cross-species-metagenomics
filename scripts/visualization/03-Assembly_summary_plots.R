##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

##############
# files to be read
##############

df_assembly <- read.csv("05_Assembly/MapToAssembly/Assembly_mapping_summary.csv", sep = ',')
df_assembly <- merge(df_assembly, df_meta, by="Sample")
df_assembly$SpeciesID <- recode(df_assembly$SpeciesID, "Am" = "Apis mellifera", "Ac" = "Apis cerana", "Af" = "Apis florea", "Ad" = "Apis dorsata")

##############
# analyse data and plot
##############

assembly_sizes <- ggplot(df_assembly, aes(y = factor(Sample, levels=samples),
                                          x = Assembly.size,
                                          fill = factor(SpeciesID, host_order)
                                        )
                                      ) +
                    geom_bar(stat = "identity") +
                    geom_vline(xintercept = 2000000) +
                      labs(x = "Total Assembly Size", y = "Sample", fill = "Host Species") +
                        scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                          make_theme(setFill = F,
                                     leg_size = 6,
                                     x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                                     y_angle=0 ,y_vj=0, y_hj=1, y_size=7,
                                     guide_nrow = 1) +
                            scale_fill_manual(values=host_order_color)
                          # ggsave("Figures/03-AssemblySizes.pdf")

df_assembly_plot <- df_assembly %>%
  select(Sample, Number.of.reads, Number.mapped) %>%
    mutate(Number.unmapped = Number.of.reads - Number.mapped) %>%
    rename(Unmapped = Number.unmapped, Mapped = Number.mapped) %>%
     pivot_longer(cols = 3:4, names_to ="Type", values_to = "Number")

number_mapped_assembly <- ggplot(df_assembly, aes(y = factor(Sample, levels=samples), x = Number.mapped, fill = factor(SpeciesID, host_order))) +
                            geom_bar(stat = "identity") +
                            labs(x = "Number of reads mapped to assembly", y = "Sample", fill = "Type") +
                            scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                              make_theme(setFill = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                                         y_angle=0 ,y_vj=0, y_hj=1, y_size=7,
                                         guide_nrow = 1) +
                                scale_fill_manual(values=host_order_color)
                              # ggsave("Figures/03-NumberAssembled.pdf")

df_assembly_plot <- df_assembly %>%
  select(Sample, Number.of.reads, Number.mapped) %>%
    mutate(Number.unmapped = Number.of.reads - Number.mapped) %>%
    mutate(percent_mapped = Number.mapped/Number.of.reads*100) %>%
    mutate(percent_unmapped = Number.unmapped/Number.of.reads*100) %>%
    rename(Unmapped = percent_unmapped, Mapped = percent_mapped, Unmapped_number = Number.unmapped, Mapped_number = Number.mapped) %>%
    select(Sample, Mapped, Unmapped) %>%
     pivot_longer(cols = 2:3, names_to ="Type", values_to = "Percent")

df_assembly_plot_number <- df_assembly %>%
  select(Sample, Number.of.reads, Number.mapped) %>%
    mutate(Number.unmapped = Number.of.reads - Number.mapped) %>%
    mutate(percent_mapped = Number.mapped/Number.of.reads*100) %>%
    mutate(percent_unmapped = Number.unmapped/Number.of.reads*100) %>%
    rename(Unmapped = Number.unmapped, Mapped = Number.mapped) %>%
    select(Sample, Mapped, Unmapped) %>%
     pivot_longer(cols = 2:3, names_to ="Type", values_to = "Number")

number_mapped_assembly <- ggplot(df_assembly_plot_number, aes(y = factor(Sample, levels=samples), x = Number, fill = factor(Type, levels = c("Unmapped", "Mapped")))) +
                            geom_bar(stat = "identity") +
                            labs(x = "Number of reads (paired end)", y = "Sample", fill = "Type") +
                            scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                              make_theme(x_angle = 60, x_hj = 1, x_vj = 1, max_colors = length(unique(df_assembly_plot$Type)), x_size = 7, guide_nrow = 1)
                              # ggsave("Figures/03-NumberMappedtoAssembly.pdf")

Number_contig_plot <- ggplot(df_assembly, aes(y = factor(Sample, levels=samples), x = Number.of.filtered.scaffolds, fill = factor(SpeciesID, levels = rev(host_order)))) +
                            geom_bar(stat = "identity") +
                            labs(y = "Sample", x = "Number of scaffolds", fill = "Host") +
                            scale_x_continuous(labels=unit_format(unit = "K", scale = 1e-4)) +
                            scale_fill_manual(values=host_order_color) +
                              make_theme(setFill = F,
                                         leg_size = 6,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                                         y_angle=0 ,y_vj=0, y_hj=1, y_size=7,
                                         guide_nrow = 1)
                              ggsave("Figures/03-NumberOfFilteredcontigs.pdf")

percent_mapped_assembly <- ggplot(df_assembly_plot, aes(y = factor(Sample, levels=samples), x = Percent, fill = factor(Type, levels = c("Unmapped", "Mapped")))) +
                            geom_bar(stat = "identity") +
                            labs(y = "Percentage of reads", x = "Sample", fill = "Type") +
                              make_theme(x_angle = 60, x_hj = 1, x_vj = 1, max_colors = length(unique(df_assembly_plot$Type)), x_size = 7, guide_nrow = 1)
                              # ggsave("Figures/03-ProportionMappedtoAssembly.pdf")

Total_data_species <- ggplot(mutate(filter(df_plot_reads, Type == "Raw" ), amount_data = Number*2*150), aes(y=factor(Sample, levels = samples),
                                      x=amount_data,
                                      fill = SpeciesID)) +
                        geom_bar(stat="identity") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(setFill = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                                         y_angle=0 ,y_vj=0, y_hj=1, y_size=7,
                                         guide_nrow = 1, leg_pos = "none") +
                                scale_x_continuous(labels=unit_format(unit = "G", scale = 1e-9)) +
                                scale_fill_manual(values=host_order_color)

Total_reads_species_temp <- ggplot(filter(df_plot_reads, Type == "Raw"), aes(y=factor(Sample, levels = samples),
                                      x=Number,
                                      fill = SpeciesID)) +
                        geom_bar(stat="identity") +
                          # ggtitle("Total reads per sample") +
                            ylab("Sample") +
                            xlab("Number of reads (paired end)") +
                              make_theme(setFill = F,
                                         x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                                         y_angle=0 ,y_vj=0, y_hj=1, y_size=7,
                                         guide_nrow = 1) +
                                scale_x_continuous(labels=unit_format(unit = "M", scale = 1e-6)) +
                                scale_fill_manual(values=host_order_color)
g <- grid.arrange(
      arrangeGrob(
        assembly_sizes + make_theme(setFill = F, setCol = F, leg_size = 6, y_size = 5, y_hj = 1, y_vj = 1, x_size = 8, guide_nrow = 1, leg_pos = "none"),
        number_mapped_assembly + make_theme(setFill = F, setCol = F, leg_size = 6, y_size = 5, y_hj = 1, y_vj = 1, x_size = 8, guide_nrow = 1, leg_pos = "none"),
        nrow = 1
      ),
      get_only_legend(assembly_sizes),
      arrangeGrob(
        Total_data_species + make_theme(setFill = F, setCol = F, leg_size = 6, y_size = 5, y_hj = 1, y_vj = 1, x_size = 8, guide_nrow = 2, leg_pos = "none"),
        Total_reads_species_temp + make_theme(setFill = F, setCol = F, leg_size = 6, y_size = 5, y_hj = 1, y_vj = 1, x_size = 8, guide_nrow = 1, leg_pos = "none"),
        nrow = 1
      ),
      nrow = 3, heights = c(8,1,8)
 )
 g
 ggsave("Figures/03-Assembly_summary.pdf", g)
Number_contig_plot