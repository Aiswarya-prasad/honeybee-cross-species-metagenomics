##############
# functions used
##############

source('scripts/visualization/utilities.R', chdir = TRUE)

##############
# files to be read
##############

qpcr_results <- c(
  "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/qPCRs/plate1/plate1-2022-11-17_115607_Results.xls",
  "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/qPCRs/plate2/plate2-2022-11-17_142003_Results.xls",
  "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/qPCRs/plate3/plate3-2022-11-17_161529_Results.xls",
  "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/qPCRs/plate4/plate4-2022-11-17_181034_Results.xls"
  # "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/qPCRs/plate5/plate5-2023-01-xx_xxxxxx_Results.xls"
  # "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/qPCRs/plate6/plate6-2023-01-xx_xxxxxx_Results.xls"
  # "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20211018_aprasad_ApisCrossSpeciesAnalysis/Analysis/qPCRs/plate7/plate7-2023-01-xx_xxxxxx_Results.xls"
)

##############
# analyse data and plot
##############

qpcr_df <- data.frame()
plate_num <- 1
for(x in qpcr_results){
  # read just the rows and columns with the values and leave out metadata cells
  # Limits the cells from rows=(36 to 132) and columns=(2 to 11)
  t <- read_excel(x, sheet = 1, range = cell_limits(c(36,2), c(132,11)))
  t <- t %>%
        mutate(`Sample Name` = ifelse(`Well Position` %in% c("H1", "H2", "H3"), paste0("Control_", plate_num), `Sample Name`)) %>%
        mutate(`Sample Name` = ifelse(`Well Position` %in% c("H4", "H5", "H6"), paste0("Control_", plate_num), `Sample Name`)) %>%
        mutate(`Sample Name` = ifelse(`Well Position` %in% c("H7", "H8", "H9"), paste0("Control_", plate_num), `Sample Name`)) %>%
        mutate(`Sample Name` = ifelse(`Well Position` %in% c("H10", "H11", "H12"), paste0("Control_", plate_num), `Sample Name`))
  # get only columns "Well Position", "Sample Name", "Target Name" and "CT" from the xls report
  t <- t[,c(1,3,4,8,9,10)]
  names(t)[1:4] <- c("Well", "Sample.Name","Target.Name","Ct")
  t$Ct <- as.numeric(t$Ct) # NAs will be introduced at the place of 'undetermined' values
  t$Target.Name <- factor(t$Target.Name)
  t$Sample.Name <- factor(t$Sample.Name)
  # t$sd <- factor(t$sd)
  t <- as.data.frame(t)
  # add in the plate number
  t <- cbind(t, "Plate" = as.factor(rep(plate_num, dim(t)[1])))
  # bind the previous plate values with the next before moving to the next file
  qpcr_df <- rbind(qpcr_df, t)
  plate_num = plate_num + 1
}
qpcr_df_info <- qpcr_df %>%
          mutate(Target.Name = ifelse(Target.Name == "UV_0356/0772", "Bacteria", "Host")) %>%
                filter(!is.na(Sample.Name)) %>%
                filter(!is.na(`Ct Mean`)) %>%
                  mutate(Host = ifelse(grepl("Control", Sample.Name), Sample.Name, NA)) %>%
                  mutate(Host = Vectorize(get_host_from_sample_name)(Sample.Name)) %>%
                  mutate(Host = as.factor(Host)) %>%
                  filter(Target.Name == "Bacteria") %>%
                  group_by(Sample.Name) %>%
                    summarise(`Ct Mean`, Host, `Ct SD`, Plate) %>%
                      group_by(Host) %>%
                      unique() %>%
                        summarise(Sample.Name, `Ct Mean`, `Ct SD`, Plate) %>%
                          rename(Ct_mean = `Ct Mean`) %>%
                          rename(Ct_SD = `Ct SD`)
ggplot(qpcr_df_info %>% filter(!is.na(Host)) %>% unique(),
       aes(x = factor(Host, levels = c(host_order, "Control_1", "Control_2", "Control_3")),
           y = Ct_mean,
           fill = factor(Host, levels = c(host_order, "Control_1", "Control_2", "Control_3")),
         )
     ) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size = 2, aes(shape = Plate)) +
  labs(y = "Ct value with UV primers (mean)", x = "Host species") +
    make_theme(leg_pos = "none", x_angle = 30, setFill = F, x_vj = 1, x_hj = 1, ) +
      scale_fill_manual(values=host_order_color)
      ggsave(Figures/00-qPCR_CT_vales.pdf)

pairwise.wilcox.test(qpcr_df_info$Ct_mean, qpcr_df_info$Host, p.adjust = "fdr")
summary(glm(data = qpcr_df_info, Ct_mean ~ Host, family = "gaussian"))
