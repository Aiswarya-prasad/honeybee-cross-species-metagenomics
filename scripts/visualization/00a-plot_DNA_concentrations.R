##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

##############
# files to be read
##############

##############
# analyse data and plot
##############

df_concentrations <- df_meta_complete %>%
  select(ID, Concentration) %>%
    left_join(df_meta, by = c("ID" = "Sample")) %>%
      mutate(Host = Vectorize(get_host_from_sample_name)(ID)) %>%
        mutate(Location = Vectorize(get_location_from_sample_name)(ID))

ggplot(df_concentrations %>% filter(ID %in% c(samples_IN, samples_MY)), 
        aes(y = factor(ID, rev(c(samples_IN, samples_MY))),
            x = Concentration,
            # fill = factor(Host, host_order),
            shape = Location,
            color = factor(Host, host_order)
          )
        ) +
  geom_point() +
  facet_wrap(~Host, scales = "free") +
  geom_hline(yintercept = 6.5, linetype = "solid") +
   geom_hline(yintercept = 14.5, linetype = "solid") +
   geom_hline(yintercept = 22.5, linetype = "solid") +
   geom_hline(yintercept = 30.5, linetype = "solid") +
    labs(color = "Host species", y = "Sample name", x = "Concentration (ng/uL) - Aliquot of 10 uL") +
    make_theme(setFill = F, setCol = F, guide_nrow = 3) +
    geom_vline(xintercept = 0, linetype = "solid", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = "solid", color = "#1a9850", alpha = 1) +
    geom_vline(xintercept = 5, linetype = "solid", color = "#fee08b", alpha = 1) +
    geom_vline(xintercept = 10, linetype = "solid", color = "#d73027", alpha = 1) +
      scale_fill_manual(values=host_order_color) +
      scale_color_manual(values=host_order_color_dark)
      ggsave("Figures/00a-DNA_Concentrations_compared.pdf")

# write.csv(df_concentrations %>% filter(ID %in% samples_MY) %>% select(ID, Concentration, Host), "Figures/Concentrations_for_plotting_edited_manually.csv", row.names = F)
df_concentrations_edited <- read.csv("Figures/Concentrations_for_plotting_edited_manually.csv")
ggplot(df_concentrations_edited %>% filter(ID %in% samples_MY) %>% mutate(Host_name = Vectorize(get_host_from_sample_name)(ID)), 
        aes(y = factor(ID, rev(samples_MY)),
            x = Concentration,
            # fill = factor(Host, host_order),
            color = factor(Host_name, host_order)
          )
        ) +
  geom_point() +
  facet_wrap(~Host, scales = "free") +
  geom_hline(yintercept = 6.5, linetype = "solid") +
   geom_hline(yintercept = 14.5, linetype = "solid") +
   geom_hline(yintercept = 22.5, linetype = "solid") +
   geom_hline(yintercept = 30.5, linetype = "solid") +
    labs(color = "Host species", y = "Sample name", x = "Concentration (ng/uL) - Aliquot of 10 uL") +
    make_theme(setFill = F, setCol = F, guide_nrow = 3) +
    geom_vline(xintercept = 0, linetype = "solid", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = "solid", color = "#d73027", alpha = 1) +
    geom_vline(xintercept = 5, linetype = "solid", color = "#fee08b", alpha = 1) +
    geom_vline(xintercept = 10, linetype = "solid", color = "#1a9850", alpha = 1) +
      scale_fill_manual(values=host_order_color) +
      scale_color_manual(values=host_order_color_dark)
      ggsave("Figures/00a-DNA_Concentrations_with_subset.pdf")

ggplot(df_concentrations %>% filter(ID %in% samples_MY),
        aes(y = factor(ID, rev(samples_MY)),
            x = Concentration,
            # fill = factor(Host, host_order),
            color = factor(Host, host_order)
          )
        ) +
  geom_point() +
  facet_wrap(~Host, scales = "free") +
  geom_hline(yintercept = 6.5, linetype = "solid") +
   geom_hline(yintercept = 14.5, linetype = "solid") +
   geom_hline(yintercept = 22.5, linetype = "solid") +
   geom_hline(yintercept = 30.5, linetype = "solid") +
    labs(color = "Host species", y = "Sample name", x = "Concentration (ng/uL) - Aliquot of 10 uL") +
    make_theme(setFill = F, setCol = F, guide_nrow = 3) +
    geom_vline(xintercept = 0, linetype = "solid", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = "solid", color = "#d73027", alpha = 1) +
    geom_vline(xintercept = 5, linetype = "solid", color = "#fee08b", alpha = 1) +
    geom_vline(xintercept = 10, linetype = "solid", color = "#1a9850", alpha = 1) +
      scale_fill_manual(values=host_order_color) +
      scale_color_manual(values=host_order_color_dark)
      ggsave("Figures/00a-DNA_Concentrations.pdf")