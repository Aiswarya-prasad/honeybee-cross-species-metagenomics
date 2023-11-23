##############
# load libraries
##############

working_dir <- "/scratch/aprasad/20230313_apis_species_comparison"
setwd(working_dir)

library(ggsignif)
library(plotly)
library(ggpubr)
library(htmlwidgets)
library(ggplot2)
library(readxl)
library(knitr)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(dplyr)
library(gridExtra)
library(ggVennDiagram)
library(vegan)
library(ape)
library(ggforce)
library(VennDiagram)
library(circlize)
library(phyloseq)
library(ggdendro)
library(dendextend)
library(ComplexHeatmap)
library(ggnewscale)

##############
# define functions used
##############

# function to edit the look of ggplot objects to make plots uniform
# add to any ggplot object with desired parameters modified
make_theme <- function(theme_name=theme_classic(), max_colors=0, palettefill="Pastel1", palettecolor="Dark2", modify_guide = T,
                        setFill=TRUE, setCol=TRUE,
                        guide_nrow=2, guide_nrow_byrow=TRUE, leg_pos="top", leg_size=12,
                        axis_x_title = 12, axis_y_title = 12,
                        x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                        y_angle=0 ,y_vj=0, y_hj=0, y_size=12){
  n_11 = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  n_12 = c("Paired", "Set3")
  n_8 = c("Accent", "Dark2", "Pastel2", "Set2")
  if (palettefill %in% n_12) {
    n_f = 12
  } else {
    if (palettefill %in% n_11) {
      n_f = 11
    } else {
      if (palettefill %in% n_8) {
        n_f  = 8
      } else {
        n_f = 9
      }
    }
  }
  if (palettecolor %in% n_12) {
    n_c = 12
  } else {
    if (palettecolor %in% n_11) {
      n_c = 11
    } else {
      if (palettecolor %in% n_8) {
        n_c  = 8
      } else {
        n_c = 9
      }
    }
  }
  getFill = colorRampPalette(brewer.pal(n_f, palettefill))
  getColor = colorRampPalette(brewer.pal(n_c, palettecolor))
  theme_params <- theme(axis.text.x = element_text(angle = x_angle,
    vjust = x_vj, hjust=x_hj,
    size = x_size),
    axis.text.y = element_text(angle = y_angle,
      vjust = y_vj, hjust=y_hj,
      size = y_size),
      axis.title.x = element_text(size=axis_x_title),
      axis.title.y = element_text(size=axis_y_title),
      legend.position=leg_pos,
      legend.text = element_text(size=leg_size)
    )
  if (modify_guide == T) {
    guide_params <- guides(fill = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  ),
                          col = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  )
                    )
  my_theme <- list(
                theme_name,
                theme_params,
                guide_params
              )
  } else {
    my_theme <- list(
                  theme_name,
                  theme_params
                )
  }
  if(setFill) {
    if (n_f < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_fill_manual(values = getFill(max_colors), na.value="grey")
                  )

    } else {
      my_theme <- list(
                    my_theme,
                    scale_fill_brewer(palette=palettefill, na.value="grey")
                  )
    }
  }
  if(setCol) {
    if (n_c < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_color_manual(values = getColor(max_colors), na.value="grey")
                  )

    } else {
      my_theme <- list(
                    my_theme,
                    scale_color_brewer(palette=palettecolor, na.value="grey")
                  )
    }
  }
  return(my_theme)
}

get_phylotype <- function(SDP){
  # only works if sdps are written in the right format
  # phylotype_x (eg. firm5_1)
  phy = strsplit(SDP, "_")[[1]][1]
  return(phy)
}

# function to be used (with Vectorize)
# to mutate in dataframes
get_host_from_colony <- function(colony_name){
  # only works if sdps are written in the right format
  # phylotype_x (eg. Am_xx)
  host_name = strsplit(colony_name, "_")[[1]][1]
  if (host_name == "Am"){
    return("Apis mellifera")
  }
  if (host_name == "Ac"){
    return("Apis cerana")
  }
  if (host_name == "Ad"){
    return("Apis dorsata")
  }
  if (host_name == "Af"){
    return("Apis florea")
  }
  if (host_name == "Aa"){
    return("Apis andreniformis")
  }
  return(NA)
}

# function to be used (with Vectorize)
# to mutate in dataframes
get_sample_name <- function(magname){
  if (startsWith(magname, "MAG_")){
    paste0(head(strsplit(strsplit(magname, "MAG_")[[1]][2], "_")[[1]], -1), collapse="_")
  } else {
    strsplit(magname, "_MAG")[[1]][1]
  }
}


# function to be used (with Vectorize)
# to mutate in dataframes
get_host_name <- function(magname){
    sample_name = strsplit(magname, "_")[[1]][1]
    if (sample_name %in% samples_to_exclude) {
    return(NA)
    }
    if (grepl("Dr|Gr", sample_name)) {
      return("Apis mellifera")
    }
    if (grepl("Am", sample_name)) {
      return("Apis mellifera")
    }
    if (grepl("Ac", sample_name)) {
      return("Apis cerana")
    }
    if (grepl("M1-|M2-|M3-|M4-|M5-|M6-|M7-|M8-|M9-|DrY|GrY|AmAi|AmIu", sample_name)) {
      return("Apis mellifera")
    }
    if (grepl("C1-|C2-|C3-|C4-|C5-|C6-|C7-|C8-|C9-|AcCh|AcKn", sample_name)) {
      return("Apis cerana")
    }
    if (grepl("D1-|D2-|D3-|D4-|D5-|D6-|D7-|D8-|D9-", sample_name)) {
      return("Apis dorsata")
    }
    if (grepl("F1-|F2-|F3-|F4-|F5-|F6-|F7-|F8-|F9-", sample_name)) {
      return("Apis florea")
    }
    if (grepl("A1-|A2-|A3-|A4-|A5-|A6-|A7-|A8-|A9-", sample_name)) {
      return("Apis andreniformis")
    }
    return(NA)
}
# function to be used (with Vectorize)
# to mutate in dataframes
get_host_from_sample_name <- function(sample_name){
  if (sample_name %in% samples_to_exclude) {
    return(NA)
  }
  if (grepl("Am", sample_name)) {
    return("Apis mellifera")
  }
  if (grepl("Ac", sample_name)) {
    return("Apis cerana")
  }
  if (grepl("M1-|M2-|M3-|M4-|M5-|M6-|M7-|M8-|M9-|DrY|GrY|AmAi|AmIu", sample_name)) {
    return("Apis mellifera")
  }
  if (grepl("C1-|C2-|C3-|C4-|C5-|C6-|C7-|C8-|C9-|AcCh|AcKn", sample_name)) {
    return("Apis cerana")
  }
  if (grepl("D1-|D2-|D3-|D4-|D5-|D6-|D7-|D8-|D9-", sample_name)) {
    return("Apis dorsata")
  }
  if (grepl("F1-|F2-|F3-|F4-|F5-|F6-|F7-|F8-|F9-", sample_name)) {
    return("Apis florea")
  }
  if (grepl("A1-|A2-|A3-|A4-|A5-|A6-|A7-|A8-|A9-", sample_name)) {
    return("Apis andreniformis")
  }
  return(sample_name)
}

# function to be used (with Vectorize)
# convert sample name in the old format to new format
convert_sample_name <- function(sample_name){
  if (grepl("AcCh|AcKn|DrY|GrY|AmAi|AmIu", sample_name)) {
    return(sample_name)
  }
  return(sub('\\.', '-', sample_name))
}

# function to be used (with Vectorize)
# to mutate in dataframes
get_location_from_sample_name <- function(sample_name){
  if (grepl("Am|Ac", sample_name)) {
    return("Japan")
  }
  if (grepl("DrY|GrY", sample_name)) {
    return("Switzerland")
  }
  if (grepl("M1-", sample_name)) {
    return("India")
  }
  if (grepl("M2-|M3-|M4-|M5-|M6-|M7-|M8-|M9-", sample_name)) {
    return("Malaysia")
  }
  if (grepl("C1-|C2-|C3-", sample_name)) {
    return("India")
  }
  if (grepl("C4-|C5-|C6-|C7-|C8-|C9-", sample_name)) {
    return("Malaysia")
  }
  if (grepl("D1-|D2-|D3-", sample_name)) {
    return("India")
  }
  if (grepl("D4-|D5-|D6-|D7-|D8-|D9-", sample_name)) {
    return("Malaysia")
  }
  if (grepl("F1-|F2-|F3-", sample_name)) {
    return("India")
  }
  if (grepl("F4-|F5-|F6-|F7-|F8-|F9-", sample_name)) {
    return("Malaysia")
  }
  if (grepl("A1-|A2-|A3-|A4-|A5-|A6-|A7-|A8-|A9-", sample_name)) {
    return("Malaysia")
  }
  return(NA)
}

# function to be used (with Vectorize)
# to mutate in dataframes
get_origin_name <- function(magname){
    sample_name = strsplit(magname, "_")[[1]][1]
    return(get_location_from_sample_name(sample_name))
}

# function to be used (with Vectorize)
# to mutate in dataframes
get_colony_from_sample_name <- function(sample_name){
  if (grepl("AmAi", sample_name)) {
    return("AmAi")
  }
  if (grepl("AmIu", sample_name)) {
    return("AmIu")
  }
  if (grepl("AcKn", sample_name)) {
    return("AcKn")
  }
  if (grepl("AcCh", sample_name)) {
    return("AcCh")
  }
  if (grepl("DrY1", sample_name)) {
    return("DrY1")
  }
  if (grepl("DrY2", sample_name)) {
    return("DrY2")
  }
  if (grepl("GrY1", sample_name)) {
    return("GrY1")
  }
  if (grepl("GrY2", sample_name)) {
    return("GrY2")
  }
  return(strsplit(sample_name, "\\.")[[1]][[1]])
}

# function to be used (with Vectorize)
# to mutate in dataframes
get_colonyid_from_sample_name <- function(sample_name){
  if (grepl("AmAi", sample_name)) {
    return("1")
  }
  if (grepl("AmIu", sample_name)) {
    return("2")
  }
  if (grepl("AcKn", sample_name)) {
    return("1")
  }
  if (grepl("AcCh", sample_name)) {
    return("2")
  }
  if (grepl("DrY1", sample_name)) {
    return("1")
  }
  if (grepl("DrY2", sample_name)) {
    return("3")
  }
  if (grepl("GrY1", sample_name)) {
    return("2")
  }
  if (grepl("GrY2", sample_name)) {
    return("4")
  }
  split_sample = strsplit(sample_name, "\\.")[[1]][[1]]
  return(strsplit(split_sample, "")[[1]][[2]])
}
get_only_legend <- function(plot) {
  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot))
  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  # extract legend
  legend <- plot_table$grobs[[legend_plot]]
  # return legend
  return(legend)
}
# function to be used (with Vectorize)
# to mutate in dataframes
format_genome_name <- function(genome){
  MAG = strsplit(genome, ".fa")[[1]]
  return(MAG)
}
# function to be used (with Vectorize)
# to mutate in dataframes
get_sdp <- function(cluster){
  if (is.na(clusters_sdp[cluster])) {
    return(cluster)
  } else{
    return(clusters_sdp[cluster])
  }
}
# extend color definitions for a vector given a named vector
# and add to it names not already defined make it grey
extend_colors <- function(names_vec, colors_vec){
  final_list <- list()
  for (a_name in names_vec) {
    if (a_name %in% names(colors_vec)) {
      final_list[a_name] = colors_vec[a_name]
    } else {
      final_list[a_name] = "grey"
    }
  }
  return(final_list)
}
extend_colors_genera <- function(names_vec) {
  extend_colors(names_vec, genusColors)
}
extend_colors_family <- function(names_vec) {
  extend_colors(names_vec, familyColors)
}
paginate_save <- function(plot, variable, plot_name, pass_nrow = 5, pass_ncol = 5, position = "top", pass_scales = "fixed"){
  for(i in 1:n_pages(plot)){
  plot +
    facet_wrap_paginate(as.formula(paste("~", variable)), ncol = pass_ncol, nrow = pass_nrow, strip.position = position, scales = pass_scales, page = i)
    ggsave(filename = paste0(strsplit(plot_name, ".pdf")[[1]][1], "_", i, '.pdf'))
  }
}
is_tsv_file <- function(file_path) {
    # Read the first few lines from the file
    lines <- readLines(file_path, n = 5)
    
    # Check if each line has consistent tab delimiter
    has_tab_delimiter <- sapply(lines, function(line) {
      count <- length(strsplit(line, "\t")[[1]])
      count > 1 && all(lengths(strsplit(line, "\t")) == count)
    })
    
    # Check if all lines have consistent tab delimiter
    is_tsv <- all(has_tab_delimiter)
    
    return(is_tsv)
  }
  get_read_count_from_flagstat <- function(flagstat, search_string = "with itself and mate mapped") {
    if (is_tsv_file(flagstat)) {
      for (line in readLines(flagstat)) {
        if (grepl(search_string, line)) {
          number <- strsplit(line, "\t")[[1]][1] %>% as.numeric()
          return(number)
        } else {
          number <- NA
        }
      }
    } else {
      for (line in readLines(flagstat)) {
        if (grepl(search_string, line)) {
          number <- strsplit(line, " +")[[1]][1] %>% as.numeric()
          return(number)
        } else {
          number <- NA
        }
      }
    }
    return(number)
  }
read_comma_numbers <- function(x) {
  as.numeric(gsub(",", "", x))
}
##############
# define important vectors (manually done)
##############

samples_to_exclude <- c("F4-5", "F5-1", "M6-2")

samples_IN <- c("M1-1", "M1-2", "M1-3", "M1-4", "M1-5",
              "C1-1", "C1-2", "C1-3", "C1-4", "C1-5",
              "C2-1", "C2-2", "C2-3", "C2-4", "C2-5",
              "C3-1", "C3-2", "C3-3", "C3-4", "C3-5",
              "D1-1","D1-2","D1-3","D1-4","D1-5",
              "D2-1","D2-2","D2-3","D2-4","D2-5",
              "D3-1","D3-2","D3-3","D3-4","D3-5",
              "F1-1","F1-2","F1-3","F1-4","F1-5",
              "F2-1","F2-2","F2-3","F2-4","F2-5",
              "F3-1","F3-2","F3-3","F3-4","F3-5"
            )
samples_IN <- samples_IN[!samples_IN %in% samples_to_exclude]
samples_KE <- c("AcCh01", "AcCh02", "AcCh03", "AcCh04", "AcCh05",
              "AcCh06", "AcCh07", "AcCh08", "AcCh09", "AcCh10",
              "AcKn01", "AcKn02", "AcKn03", "AcKn04", "AcKn05",
              "AcKn06", "AcKn07", "AcKn08", "AcKn09", "AcKn10",
              "AmAi01", "AmAi02", "AmAi03", "AmAi04", "AmAi05",
              "AmAi06", "AmAi07", "AmAi08", "AmAi09", "AmAi10",
              "AmIu01", "AmIu02", "AmIu03", "AmIu04", "AmIu05",
              "AmIu06", "AmIu07", "AmIu08", "AmIu09", "AmIu10",
              "DrY1_F1", "DrY1_F2", "DrY1_F3", "DrY1_F4", "DrY1_F5", "DrY1_F6",
              "DrY1_N1", "DrY1_N2", "DrY1_N3", "DrY1_N4", "DrY1_N5", "DrY1_N6",
              "DrY1_W1", "DrY1_W2", "DrY1_W3", "DrY1_W4", "DrY1_W5", "DrY1_W6",
              "DrY2_F1", "DrY2_F2", "DrY2_F3", "DrY2_F4", "DrY2_F5",
              "DrY2_N1", "DrY2_N2", "DrY2_N3", "DrY2_N4", "DrY2_N5", "DrY2_N6",
              "GrY2_F1", "GrY2_F2", "GrY2_F3", "GrY2_F4", "GrY2_F5", "GrY2_F6",
              "GrY2_N1", "GrY2_N2", "GrY2_N3", "GrY2_N4", "GrY2_N5", "GrY2_N6",
              "GrY2_W1", "GrY2_W2", "GrY2_W3", "GrY2_W4", "GrY2_W5", "GrY2_W6"
            )
samples_JP <- c("AmAi01", "AmAi02", "AmAi03", "AmAi04", "AmAi05",
              "AmAi06", "AmAi07", "AmAi08", "AmAi09", "AmAi10",
              "AmIu01", "AmIu02", "AmIu03", "AmIu04", "AmIu05",
              "AmIu06", "AmIu07", "AmIu08", "AmIu09", "AmIu10",
              "AcKn01", "AcKn02", "AcKn03", "AcKn04", "AcKn05",
              "AcKn06", "AcKn07", "AcKn08", "AcKn09", "AcKn10",
              "AcCh01", "AcCh02", "AcCh03", "AcCh04", "AcCh05",
              "AcCh06", "AcCh07", "AcCh08", "AcCh09", "AcCh10"
            )
samples_SW <- c("DrY1_F1", "DrY1_F2", "DrY1_F3", "DrY1_F4", "DrY1_F5", "DrY1_F6",
              "DrY1_N1", "DrY1_N2", "DrY1_N3", "DrY1_N4", "DrY1_N5", "DrY1_N6",
              "DrY1_W1", "DrY1_W2", "DrY1_W3", "DrY1_W4", "DrY1_W5", "DrY1_W6",
              "DrY2_F1", "DrY2_F2", "DrY2_F3", "DrY2_F4", "DrY2_F5",
              "DrY2_N1", "DrY2_N2", "DrY2_N3", "DrY2_N4", "DrY2_N5", "DrY2_N6",
              "GrY2_F1", "GrY2_F2", "GrY2_F3", "GrY2_F4", "GrY2_F5", "GrY2_F6",
              "GrY2_N1", "GrY2_N2", "GrY2_N3", "GrY2_N4", "GrY2_N5", "GrY2_N6",
              "GrY2_W1", "GrY2_W2", "GrY2_W3", "GrY2_W4", "GrY2_W5", "GrY2_W6"
            )
samples_MY <- c("M2-1", "M2-2", "M2-3", "M2-4", "M2-5",
              "M3-1", "M3-2", "M3-3", "M3-4", "M3-5",
              "M4-1", "M4-2", "M4-3", "M4-4", "M4-5",
              "M5-1", "M5-2", "M5-3", "M5-4", "M5-5",
              "M6-1", "M6-2", "M6-3", "M6-4", "M6-5",
              "M7-1", "M7-2", "M7-3", "M7-4", "M7-5",
              "C4-1", "C4-2", "C4-3", "C4-4", "C4-5",
              "C5-1", "C5-2", "C5-3", "C5-4", "C5-5",
              "C6-1", "C6-2", "C6-3", "C6-4", "C6-5",
              "C7-1", "C7-2", "C7-3", "C7-4", "C7-5",
              "C8-1", "C8-2", "C8-3", "C8-4", "C8-5",
              "C9-1", "C9-2", "C9-3", "C9-4", "C9-5",
              "D4-1","D4-2","D4-3","D4-4","D4-5",
              "D5-1","D5-2","D5-3","D5-4","D5-5",
              "D6-1","D6-2","D6-3","D6-4","D6-5",
              "D7-1","D7-2","D7-3","D7-4","D7-5",
              "D8-1","D8-2","D8-3","D8-4","D8-5",
              "D9-1","D9-2","D9-3","D9-4","D9-5",
              "F4-1","F4-2","F4-3","F4-4","F4-5",
              "F5-1","F5-2","F5-3","F5-4","F5-5",
              "F6-1","F6-2","F6-3","F6-4","F6-5",
              "F7-1","F7-2","F7-3","F7-4","F7-5",
              "F8-1","F8-2","F8-3","F8-4","F8-5",
              "F9-1","F9-2","F9-3","F9-4","F9-5",
              "A1-1","A1-2","A1-3","A1-4","A1-5",
              "A2-1","A2-2","A2-3","A2-4","A2-5",
              "A3-1","A3-2","A3-3","A3-4","A3-5",
              "A4-1","A4-2","A4-3","A4-4","A4-5",
              "A5-1","A5-2","A5-3","A5-4","A5-5",
              "A6-1","A6-2","A6-3","A6-4","A6-5"
            )
samples_MY <- samples_MY[!samples_MY %in% samples_to_exclude]
samples_AM <- c("M1-1", "M1-2", "M1-3", "M1-4", "M1-5",
              "M2-1", "M2-2", "M2-3", "M2-4", "M2-5",
              "M3-1", "M3-2", "M3-3", "M3-4", "M3-5",
              "M4-1", "M4-2", "M4-3", "M4-4", "M4-5",
              "M5-1", "M5-2", "M5-3", "M5-4", "M5-5",
              "M6-1", "M6-2", "M6-3", "M6-4", "M6-5",
              "M7-1", "M7-2", "M7-3", "M7-4", "M7-5",
              "AmAi01", "AmAi02", "AmAi03", "AmAi04", "AmAi05",
              "AmAi06", "AmAi07", "AmAi08", "AmAi09", "AmAi10",
              "AmIu01", "AmIu02", "AmIu03", "AmIu04", "AmIu05",
              "DrY1_W1", "DrY1_W2", "DrY1_W3", "DrY1_W4", "DrY1_W5", "DrY1_W6",
              "DrY1_F1", "DrY1_F2", "DrY1_F3", "DrY1_F4", "DrY1_F5", "DrY1_F6",
              "DrY1_N1", "DrY1_N2", "DrY1_N3", "DrY1_N4", "DrY1_N5", "DrY1_N6",
              "DrY2_F1", "DrY2_F2", "DrY2_F3", "DrY2_F4", "DrY2_F5",
              "DrY2_N1", "DrY2_N2", "DrY2_N3", "DrY2_N4", "DrY2_N5", "DrY2_N6",
              "GrY2_F1", "GrY2_F2", "GrY2_F3", "GrY2_F4", "GrY2_F5", "GrY2_F6",
              "GrY2_N1", "GrY2_N2", "GrY2_N3", "GrY2_N4", "GrY2_N5", "GrY2_N6",
              "GrY2_W1", "GrY2_W2", "GrY2_W3", "GrY2_W4", "GrY2_W5", "GrY2_W6"
            )
samples_AM <- samples_AM[!samples_AM %in% samples_to_exclude]
samples_AC <- c("C1-1", "C1-2", "C1-3", "C1-4", "C1-5",
              "C2-1", "C2-2", "C2-3", "C2-4", "C2-5",
              "C3-1", "C3-2", "C3-3", "C3-4", "C3-5",
              "C4-1", "C4-2", "C4-3", "C4-4", "C4-5",
              "C5-1", "C5-2", "C5-3", "C5-4", "C5-5",
              "C6-1", "C6-2", "C6-3", "C6-4", "C6-5",
              "C7-1", "C7-2", "C7-3", "C7-4", "C7-5",
              "C8-1", "C8-2", "C8-3", "C8-4", "C8-5",
              "C9-1", "C9-2", "C9-3", "C9-4", "C9-5",
              "AcKn01", "AcKn02", "AcKn03", "AcKn04", "AcKn05",
              "AcKn06", "AcKn07", "AcKn08", "AcKn09", "AcKn10",
              "AcCh01", "AcCh02", "AcCh03", "AcCh04", "AcCh05",
              "AcCh06", "AcCh07", "AcCh08", "AcCh09", "AcCh10"
            )
samples_AD <- c("D1-1","D1-2","D1-3","D1-4","D1-5",
              "D2-1","D2-2","D2-3","D2-4","D2-5",
              "D3-1","D3-2","D3-3","D3-4","D3-5",
              "D4-1","D4-2","D4-3","D4-4","D4-5",
              "D5-1","D5-2","D5-3","D5-4","D5-5",
              "D6-1","D6-2","D6-3","D6-4","D6-5",
              "D7-1","D7-2","D7-3","D7-4","D7-5",
              "D8-1","D8-2","D8-3","D8-4","D8-5",
              "D9-1","D9-2","D9-3","D9-4","D9-5"
            )
samples_AF <- c("F1-1","F1-2","F1-3","F1-4","F1-5",
              "F2-1","F2-2","F2-3","F2-4","F2-5",
              "F3-1","F3-2","F3-3","F3-4","F3-5",
              "F4-1","F4-2","F4-3","F4-4","F4-5",
              "F5-1","F5-2","F5-3","F5-4","F5-5",
              "F6-1","F6-2","F6-3","F6-4","F6-5",
              "F7-1","F7-2","F7-3","F7-4","F7-5",
              "F8-1","F8-2","F8-3","F8-4","F8-5",
              "F9-1","F9-2","F9-3","F9-4","F9-5"
            )
samples_AF <- samples_AF[!samples_AF %in% samples_to_exclude]
samples_AA <- c("A1-1","A1-2","A1-3","A1-4","A1-5",
              "A2-1","A2-2","A2-3","A2-4","A2-5",
              "A3-1","A3-2","A3-3","A3-4","A3-5",
              "A4-1","A4-2","A4-3","A4-4","A4-5",
              "A5-1","A5-2","A5-3","A5-4","A5-5",
              "A6-1","A6-2","A6-3","A6-4","A6-5"
            )
samples_IN_MY <- c("M1-1", "M1-2", "M1-3", "M1-4", "M1-5",
              "M2-1", "M2-2", "M2-3", "M2-4", "M2-5",
              "M3-1", "M3-2", "M3-3", "M3-4", "M3-5",
              "M4-1", "M4-2", "M4-3", "M4-4", "M4-5",
              "M5-1", "M5-2", "M5-3", "M5-4", "M5-5",
              "M6-1", "M6-2", "M6-3", "M6-4", "M6-5",
              "M7-1", "M7-2", "M7-3", "M7-4", "M7-5",
              "C1-1", "C1-2", "C1-3", "C1-4", "C1-5",
              "C2-1", "C2-2", "C2-3", "C2-4", "C2-5",
              "C3-1", "C3-2", "C3-3", "C3-4", "C3-5",
              "C4-1", "C4-2", "C4-3", "C4-4", "C4-5",
              "C5-1", "C5-2", "C5-3", "C5-4", "C5-5",
              "C6-1", "C6-2", "C6-3", "C6-4", "C6-5",
              "C7-1", "C7-2", "C7-3", "C7-4", "C7-5",
              "C8-1", "C8-2", "C8-3", "C8-4", "C8-5",
              "C9-1", "C9-2", "C9-3", "C9-4", "C9-5",
              "D1-1","D1-2","D1-3","D1-4","D1-5",
              "D2-1","D2-2","D2-3","D2-4","D2-5",
              "D3-1","D3-2","D3-3","D3-4","D3-5",
              "D4-1","D4-2","D4-3","D4-4","D4-5",
              "D5-1","D5-2","D5-3","D5-4","D5-5",
              "D6-1","D6-2","D6-3","D6-4","D6-5",
              "D7-1","D7-2","D7-3","D7-4","D7-5",
              "D8-1","D8-2","D8-3","D8-4","D8-5",
              "D9-1","D9-2","D9-3","D9-4","D9-5",
              "F1-1","F1-2","F1-3","F1-4","F1-5",
              "F2-1","F2-2","F2-3","F2-4","F2-5",
              "F3-1","F3-2","F3-3","F3-4","F3-5",
              "F4-1","F4-2","F4-3","F4-4","F4-5",
              "F5-1","F5-2","F5-3","F5-4","F5-5",
              "F6-1","F6-2","F6-3","F6-4","F6-5",
              "F7-1","F7-2","F7-3","F7-4","F7-5",
              "F8-1","F8-2","F8-3","F8-4","F8-5",
              "F9-1","F9-2","F9-3","F9-4","F9-5",
              "A1-1","A1-2","A1-3","A1-4","A1-5",
              "A2-1","A2-2","A2-3","A2-4","A2-5",
              "A3-1","A3-2","A3-3","A3-4","A3-5",
              "A4-1","A4-2","A4-3","A4-4","A4-5",
              "A5-1","A5-2","A5-3","A5-4","A5-5",
              "A6-1","A6-2","A6-3","A6-4","A6-5"
            )     
samples_IN_MY <- samples_IN_MY[!samples_IN_MY %in% samples_to_exclude]
samples <- c(samples_AM, samples_AC, samples_AD, samples_AF, samples_AA)
location_order <- c("Malaysia", "India", "Switzerland", "Japan")
location_order_color <- c("Malaysia" = "#1f78b4", "India" = "#33a02c", "Switzerland" = "#e31a1c", "Japan" = "#fb9a99")
host_order <- c("Apis mellifera", "Apis cerana", "Apis dorsata", "Apis florea", "Apis andreniformis")
host_order_color <- c("Apis mellifera" = brewer.pal(9, "Pastel1")[2], "Apis cerana" = brewer.pal(9, "Pastel1")[1], "Apis dorsata" = brewer.pal(9, "Pastel1")[4], "Apis florea" = brewer.pal(9, "Pastel1")[3], "Apis andreniformis" = brewer.pal(9, "Pastel1")[5])
host_order_color_dark <- c("Apis mellifera" = brewer.pal(9, "Set1")[2], "Apis cerana" = brewer.pal(9, "Set1")[1], "Apis dorsata" = brewer.pal(9, "Set1")[4], "Apis florea" = brewer.pal(9, "Set1")[3], "Apis andreniformis" = brewer.pal(9, "Set1")[5])
location_country_colors <- c("India" = brewer.pal(9, "Set3")[1], "Switzerland" = brewer.pal(9, "Pastel1")[8], "Malaysia" = brewer.pal(9, "Pastel1")[6], "Japan" = brewer.pal(9, "Pastel1")[7])
colony_order <- c("M_1", "M_Iu", "M_Ai", "M_DrY2_F", "C_1", "C_2", "C_3", "C_Kn", "C_Ch", "D_1", "D_2", "D_3", "F_1", "F_2", "F_3")
location_order <- c("AIST_Am", "UT_Am", "Bee park, GKVK_Am","Les Droites_Am",
                    "NCBS campus_Ac", "Bee park, GKVK_Ac", "Chiba_Ac", "Kanagawa_Ac",
                    "Biological sciences building, IISc_Ad","House near NCBS_Ad","Naideli hostel_Ad",
                    "Bangalore outskirts_Af")


genera <- c(
    "g__Lactobacillus",
    "g__Bifidobacterium",
    "g__Bombilactobacillus",
    "g__CALYQJ01",
    "g__Snodgrassella",
    "g__Gilliamella",
    "g__CALYQQ01",
    "g__Apibacter",
    "g__Bartonella_A",
    "g__Frischella",
    "g__Commensalibacter",
    "g__Apilactobacillus",
    "g__Bombella",
    "g__Dysgonomonas",
    "g__Enterobacter",
    "g__Pectinatus",
    "g__Spiroplasma",
    "g__Zymobacter",
    "g__Entomomonas",
    "g__Saezia",
    "g__Parolsenella",
    "g__WRHT01",
    "g__Cronobacter",
    "g__Melissococcus",
    "g__Rosenbergiella",
    "g__Pluralibacter",
    "g__KACC-21507",
    "g__"
  )         
genusColors <- list("g__Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                    "g__Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                    "g__CALYQJ01" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                    "g__Dysgonomonas" = brewer.pal(11, "Spectral")[2],
                    "g__Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                    "g__Apibacter" = brewer.pal(11, "Spectral")[4],
                    "g__Bombella" = brewer.pal(11, "Spectral")[5],
                    "g__Commensalibacter" = brewer.pal(11, "Spectral")[6],
                    "g__Bartonella_A" = brewer.pal(11, "Spectral")[7],
                    "g__Frischella" = brewer.pal(11, "Spectral")[8],
                    "g__Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                    "g__Snodgrassella" = brewer.pal(11, "Spectral")[10],
                    "g__Gilliamella" = brewer.pal(11, "Spectral")[11],
                    "g__CALYQQ01" = brewer.pal(11, "Spectral")[11],
                    "g__Spiroplasma" = brewer.pal(8, "Set1")[8],
                    "g__Pectinatus" = brewer.pal(8, "Dark2")[1],
                    "g__WRHT01" = brewer.pal(8, "Dark2")[3],
                    "g__Saezia" = brewer.pal(8, "Dark2")[6],
                    "g__Enterobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[1],
                    "g__Q3-R57-64" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[1],
                    "g__Zymobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[2],
                    "g__Entomomonas"= head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[4],
                    "g__Escherichia" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
                    "g__Klebsiella" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
                    "g__Parolsenella" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Cronobacter" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[2],
                    "g__Pluralibacter" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[4],
                    "g__Rosenbergiella" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[5],
                    "g__Melissococcus" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[6],
                    "g__KACC-21507" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[7],
                    "g__Pantoea" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Acinetobacter" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Corynebacterium" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Kosakonia" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__JAATFO01" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Hafnia" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Zymomonas" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Fructobacillus" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Floricoccus" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__" = "#000000"
)
# genusColors <- list("g__Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
#                     "g__Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
#                     "g__CALYQJ01" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
#                     "g__Dysgonomonas" = brewer.pal(11, "Spectral")[2],
#                     "g__Bifidobacterium" = brewer.pal(11, "Spectral")[3],
#                     "g__Apibacter" = brewer.pal(11, "Spectral")[4],
#                     "g__Bombella" = brewer.pal(11, "Spectral")[5],
#                     "g__Commensalibacter" = brewer.pal(11, "Spectral")[6],
#                     "g__Bartonella_A" = brewer.pal(11, "Spectral")[7],
#                     "g__Frischella" = brewer.pal(11, "Spectral")[8],
#                     "g__Apilactobacillus" = brewer.pal(11, "Spectral")[9],
#                     "g__Snodgrassella" = brewer.pal(11, "Spectral")[10],
#                     "g__Gilliamella" = brewer.pal(11, "Spectral")[11],
#                     "g__CALYQQ01" = brewer.pal(11, "Spectral")[11],
#                     "g__Spiroplasma" = brewer.pal(8, "Set1")[8],
#                     "g__Pectinatus" = brewer.pal(8, "Dark2")[1],
#                     "g__WRHT01" = brewer.pal(8, "Dark2")[3],
#                     "g__Saezia" = brewer.pal(8, "Dark2")[6],
#                     "g__Enterobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[1],
#                     "g__Q3-R57-64" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[1],
#                     "g__Zymobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[2],
#                     "g__Entomomonas"= head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[4],
#                     "g__Escherichia" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
#                     "g__Klebsiella" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
#                     "g__Parolsenella" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[8],
#                     "g__Cronobacter" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[2],
#                     "g__Pluralibacter" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[4],
#                     "g__Rosenbergiella" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[5],
#                     "g__Melissococcus" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[6],
#                     "g__KACC-21507" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[7],
#                     "g__" = "#000000"
# )


genusColors_char <- c("g__Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                    "g__Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                    "g__CALYQJ01" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                    "g__Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                    "g__Gilliamella" = brewer.pal(11, "Spectral")[11],
                    "g__CALYQQ01" = brewer.pal(11, "Spectral")[11],
                    "g__Frischella" = brewer.pal(11, "Spectral")[8],
                    "g__Bartonella_A" = brewer.pal(11, "Spectral")[7],
                    "g__Snodgrassella" = brewer.pal(11, "Spectral")[10],
                    "g__Apibacter" = brewer.pal(11, "Spectral")[4],
                    "g__Commensalibacter" = brewer.pal(11, "Spectral")[6],
                    "g__Bombella" = brewer.pal(11, "Spectral")[5],
                    "g__Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                    "g__Dysgonomonas" = brewer.pal(11, "Spectral")[2],
                    "g__Spiroplasma" = brewer.pal(8, "Set1")[8],
                    "g__WRHT01" = brewer.pal(8, "Dark2")[5],
                    "g__Pectinatus" = brewer.pal(8, "Dark2")[1],
                    "g__Enterobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[1],
                    "g__Q3-R57-64" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[1],
                    "g__Zymobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[2],
                    "g__Entomomonas"= head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[4],
                    "g__Escherichia" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
                    "g__Klebsiella" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
                    "g__Saezia" = brewer.pal(8, "Dark2")[6],
                    "g__Parolsenella" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__Cronobacter" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[2],
                    "g__Pluralibacter" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[4],
                    "g__Rosenbergiella" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[5],
                    "g__Melissococcus" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[6],
                    "g__KACC-21507" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[7],
                    "g__Pantoea" = head(colorRampPalette(c(brewer.pal(11, "PiYG")[2], "#FFFFFF"))(10), -1)[8],
                    "g__" = "#000000"
)

families <- c("f__Lactobacillaceae", "f__Bifidobacteriaceae", "f__Enterobacteriaceae", "f__Neisseriaceae", "f__Rhizobiaceae_A", "f__Selenomonadaceae", "f__Weeksellaceae", "f__Dysgonomonadaceae", "f__Mycoplasmataceae", "f__Halomonadaceae", "f__Pseudomonadaceae", "f__Burkholderiaceae", "f__Atopobiaceae", "f__Desulfovibrionaceae", "f__Acetobacteraceae", "f__", "f__Streptococcaceae")
familyColors <- list(
  "f__Lactobacillaceae" = brewer.pal(11, "Spectral")[1],
  "f__Bifidobacteriaceae" = brewer.pal(11, "Spectral")[3],
  "f__Enterobacteriaceae" = brewer.pal(11, "Spectral")[11],
  "f__Neisseriaceae" = brewer.pal(11, "Spectral")[10],
  "f__Rhizobiaceae_A" = brewer.pal(11, "Spectral")[7],
  "f__Weeksellaceae" = brewer.pal(11, "Spectral")[4],
  "f__Acetobacteraceae" = brewer.pal(11, "Spectral")[6],
  "f__Dysgonomonadaceae" = brewer.pal(11, "Spectral")[2],
  "f__Mycoplasmataceae" = brewer.pal(8, "Set1")[8],
  "f__Desulfovibrionaceae" = brewer.pal(8, "Dark2")[3],
  "f__Selenomonadaceae" = brewer.pal(8, "Dark2")[1],
  "f__Halomonadaceae" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[2],
  "f__Pseudomonadaceae" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[4],
  "f__Burkholderiaceae" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
  "f__Atopobiaceae" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[8],
  "f__Streptococcaceae" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[9],
  "f__" = "#000000"
)
#
motus <- c("Lactobacillus", "Bifidobacterium", "Gilliamella", "Frischella", "Snodgrassella", "Bartonella", "Bombella", "Acetobacteraceae", "Dysgonomonadaceae", "Spiroplasma", "Flavobacteriaceae", "Fructobacillus", "unassigned")
motuColors <- list(
  "Lactobacillus" = brewer.pal(11, "Spectral")[1],
  "Bifidobacterium" = brewer.pal(11, "Spectral")[3],
  "Flavobacteriaceae" = brewer.pal(11, "Spectral")[4],
  "Bombella" = brewer.pal(11, "Spectral")[5],
  "Acetobacteraceae" = brewer.pal(11, "Spectral")[6],
  "Bartonella" = brewer.pal(11, "Spectral")[7],
  "Frischella" = brewer.pal(11, "Spectral")[8],
  "Gilliamella" = brewer.pal(11, "Spectral")[11],
  "Snodgrassella" = brewer.pal(11, "Spectral")[10],
  "Dysgonomonadaceae" = brewer.pal(11, "Spectral")[2],
  "Fructobacillus" = brewer.pal(11, "Spectral")[9],
  "Spiroplasma" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[2], "#FFFFFF"))(10), -1)[8],
  "unassigned" = "black"
)
#
# This function is from https://github.com/Russel88/MicEco/blob/91f8e6f5d67e0bfd018dc3b46da3994b1eadeb46/R/adonis_OmegaSq.R#L10
#' Calculate (partial) Omega-squared (effect-size calculation) for PERMANOVA and add it to the input object
#'
#' @param adonisOutput An adonis object
#' @param partial Should partial omega-squared be calculated (sample size adjusted). Default TRUE
#' @return Original adonis object with the (partial) Omega-squared values added
#' @import vegan
#' @export
adonis_OmegaSq <- function(adonisOutput, partial = TRUE){
    if(!(is(adonisOutput, "adonis") || is(adonisOutput, "anova.cca")))
        stop("Input should be an adonis object")
    if (is(adonisOutput, "anova.cca")) {
        aov_tab <- adonisOutput
        aov_tab$MeanSqs <- aov_tab$SumOfSqs / aov_tab$Df
        aov_tab$MeanSqs[length(aov_tab$Df)] <- NA
    } else {
        aov_tab <- adonisOutput$aov.tab
    }
    heading <- attr(aov_tab, "heading")
    MS_res <- aov_tab[pmatch("Residual", rownames(aov_tab)), "MeanSqs"]
    SS_tot <- aov_tab[rownames(aov_tab) == "Total", "SumsOfSqs"]
    N <- aov_tab[rownames(aov_tab) == "Total", "Df"] + 1
    if(partial){
        omega <- apply(aov_tab, 1, function(x) (x["Df"]*(x["MeanSqs"]-MS_res))/(x["Df"]*x["MeanSqs"]+(N-x["Df"])*MS_res))
        aov_tab$parOmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    } else {
        omega <- apply(aov_tab, 1, function(x) (x["SumsOfSqs"]-x["Df"]*MS_res)/(SS_tot+MS_res))
        aov_tab$OmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    }
    if (is(adonisOutput, "adonis"))
        cn_order <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2",
                      if (partial) "parOmegaSq" else "OmegaSq", "Pr(>F)")
    else
        cn_order <- c("Df", "SumOfSqs", "F", if (partial) "parOmegaSq" else "OmegaSq",
                      "Pr(>F)")
    aov_tab <- aov_tab[, cn_order]
    attr(aov_tab, "names") <- cn_order
    attr(aov_tab, "heading") <- heading
    if (is(adonisOutput, "adonis"))
        adonisOutput$aov.tab <- aov_tab
    else
        adonisOutput <- aov_tab
    return(adonisOutput)
}
