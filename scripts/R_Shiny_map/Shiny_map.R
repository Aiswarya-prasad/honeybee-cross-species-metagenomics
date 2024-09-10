# install on conda using accompanying yaml file
# make sure that you have already run the two before
# conda config --add channels new_channel
# for channels bioconda and conda-forge
# use command,
# conda env create -f rmd-shiny-env.yaml
# or
# cd /...<project_path>.../20211018_aprasad_ApisCrossSpeciesAnalysis/DataCollection/R_Shiny_map
# if (!require("shiny")) install.packages("shiny")
# if (!require("dplyr")) install.packages("dplyr")
# if (!require("leaflet")) install.packages("leaflet")
# if (!require("DT")) install.packages("DT")
# if (!require("xlsx")) install.packages("xlsx")
# conda activate rmd-shiny-env

#####
# run using
# R -e "shiny::runApp('Shiny_map.R')"
#####

require(shiny)
require(dplyr)
require(leaflet)
require(DT)
require(xlsx)
require(RColorBrewer)

parse_ID <- function(ID){
  ID = as.character(ID)
  date = strsplit(ID, split = "_")[[1]][1]
  person = strsplit(ID, split = "_")[[1]][2]
  tmp = strsplit(ID, split = "_")[[1]][3]
  colony = strsplit(tmp, split = "-")[[1]][1]
  bee = strsplit(tmp, split = "-")[[1]][2]
  info = list(
    "date" = date,
    "person" = person,
    "colony" = colony,
    "bee" = bee
  )
  return(info)
}

transform_species_names <- function(long_name){
  if (grepl("mellifera", long_name)) {
    return("Apis mellifera")
  }
  if (grepl("A. m", long_name)) {
    return("Apis mellifera")
  }
  if (grepl("cerana", long_name)) {
    return("Apis cerana")
  }
  if (grepl("dorsata", long_name)) {
    return("Apis dorsata")
  }
  if (grepl("florea", long_name)) {
    return("Apis florea")
  }
  if (grepl("andreniformis", long_name)) {
    return("Apis andreniformis")
  }
  return(long_name)
}

transform_species_names_mel <- function(long_name){
  if (grepl("mellifera", long_name)) {
    return("Apis mellifera")
  }
  if (grepl("littorea", long_name)) {
    return("A. m. littorea")
  }
  if (grepl("jemenica", long_name)) {
    return("A. m. jemenetica")
  }
  if (grepl("monticola", long_name)) {
    return("A. m. monticola")
  }
  if (grepl("scutellata", long_name)) {
    return("A. m. scutellata")
  }
  if (grepl("simensis", long_name)) {
    return("A. m. simensis")
  }
  if (grepl("lamarkii", long_name)) {
    return("A. m. lamarkii")
  }
  if (grepl("cerana", long_name)) {
    return("Apis other")
  }
  if (grepl("dorsata", long_name)) {
    return("Apis other")
  }
  if (grepl("florea", long_name)) {
    return("Apis other")
  }
  if (grepl("andreniformis", long_name)) {
    return("Apis other")
  }
  return(long_name)
}

server <- function(input, output) {
  # Import Data and clean it
  metadata_shiny <- read.xlsx(paste0(getwd(), "/221220_metadata_compiled.xlsx"), 1)[-1, ]
  metadata_shiny <- metadata_shiny[!duplicated(metadata_shiny$Colony_ID),]
  metadata_shiny <- as.data.frame(metadata_shiny) %>%
                      mutate(Species_long = Species) %>%
                        # mutate(Species = Vectorize(transform_species_names_mel)(Species_long))
                        mutate(Species = Vectorize(transform_species_names)(Species_long))
  metadata_shiny$Location_latitude <- as.numeric(metadata_shiny$Location_latitude)
  metadata_shiny$Location_longitude <- as.numeric(metadata_shiny$Location_longitude)
  metadata_shiny <- metadata_shiny[which(!is.na(metadata_shiny$Location_latitude)), ]
  metadata_shiny <- metadata_shiny[which(!is.na(metadata_shiny$Location_longitude)), ]
  # modified (smaller popup label)
  metadata_shiny <- dplyr::mutate(metadata_shiny, cntnt=paste0('<strong>Name (eg. individual): </strong>',ID,
                                          '<br><strong>Species:</strong> ',Species,
                                          '<br><strong>Colony_type:</strong> ',Colony_type,
                                          '<br><strong>Location_name:</strong> ',Location_name,
                                          '<br><strong>Colony_ID:</strong> ', Colony_ID,
                                          '<br><strong>Altitude:</strong> ',Location_altitude,
                                          '<br><strong>Location type:</strong> ',Location_type,
                                          '<br><strong>Notes:</strong> ',Notes_1))
  # new column for the popup label
  # metadata_shiny <- dplyr::mutate(metadata_shiny, cntnt=paste0('<strong>Name (eg. individual): </strong>',ID,
  #                                         '<br><strong>Species:</strong> ',Species,
  #                                         '<br><strong>Colony_type:</strong> ',Colony_type,
  #                                         '<br><strong>Location_name:</strong> ',Location_name,
  #                                         '<br><strong>Collected by:</strong> ',parse_ID(ID)$person,
  #                                         '<br><strong>Collected on:</strong> ',paste0(Collection_date_d,"/", Collection_date_m,"/", Collection_date_y),
  #                                         '<br><strong>Colony_ID:</strong> ', Colony_ID,
  #                                         '<br><strong>Extraction:</strong> ',DNA_extracted,
  #                                         '<br><strong>Short_Name:</strong> ', Short_name,
  #                                         '<br><strong>Hive_description:</strong> ',Hive_description,
  #                                         '<br><strong>Altitude:</strong> ',Location_altitude,
  #                                         '<br><strong>Location type:</strong> ',Location_type,
  #                                         '<br><strong>Location latitude:</strong> ',Location_latitude,
  #                                         '<br><strong>Location longitude:</strong> ',Location_longitude,
  #                                         '<br><strong>Sequenced on:</strong> ',Sequenced_on,
  #                                         '<br><strong>Notes:</strong> ',Notes_1))
  pal <- colorFactor(pal = c("#ff8c00", "#ec008c", "#e41a1c", "#00188f", "#68217a", "#00bcf2", "#00b294", "#4daf4a", "#bad80a", "#377eb8"), domain = unique(metadata_shiny$Species))
  # only species not subspecies:
  pal <- colorFactor(pal = c("#ff8c00", "#ec008c", "#e41a1c", "#00188f", "#68217a", "#00bcf2", "#00b294", "#4daf4a", "#bad80a", "#377eb8"), domain = unique(metadata_shiny$Species))
  # pal <- colorFactor(pal = c("#fed9a6", "#ec008c", "#fbb4ae", "#00188f", "#decbe4", "#00bcf2", "#00b294", "#ccebc5", "#bad80a", "#b3cde3"), domain = unique(metadata_shiny$Species))
  # filter based on slider input
  my_icons <- iconList(
  "W" <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/black-circle-icon-23.png",
                          iconWidth = 18, iconHeight = 18),
  "M" <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/black-square-frame-23.png",
                          iconWidth = 18, iconHeight = 18),
  "SW" <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/triangle-png-28.png",
                            iconWidth = 18, iconHeight = 18)
      )
  # metadata_shiny_filt <- reactive({metadata_shiny})
  metadata_shiny_filt <- reactive({
    metadata_shiny$Collection_date_y <- as.numeric(metadata_shiny$Collection_date_y)
    df <- metadata_shiny[which(metadata_shiny$Collection_date_y >= input$Year[1] & metadata_shiny$Collection_date_y <= input$Year[2]),]
    year1 <- as.numeric(format(input$dateRange[1], "%y"))
    month1 <- as.numeric(format(input$dateRange[1], "%m"))
    day1 <- as.numeric(format(input$dateRange[1], "%d"))
    year2 <- as.numeric(format(input$dateRange[2], "%y"))
    month2 <- as.numeric(format(input$dateRange[2], "%m"))
    day2 <- as.numeric(format(input$dateRange[2], "%d"))
    start_date <- as.Date(paste0(day1, "/", month1, "/", year1), "%d/%m/%y")
    end_date <- as.Date(paste0(day2, "/", month2, "/", year2), "%d/%m/%y")
    df <- metadata_shiny %>%
            mutate(collection_date = paste0(Collection_date_d, "/", Collection_date_m, "/", Collection_date_y)) %>%
              filter(as.Date(collection_date, "%d/%m/%y") >= start_date & as.Date(collection_date, "%d/%m/%y") <= end_date)
    df
  })

  # create the leaflet map
  output$metamap <- renderLeaflet({
        leaflet(metadata_shiny_filt()) %>%
        addCircles(lng = ~Location_longitude, lat = ~Location_latitude) %>%
        addTiles() %>%
        addScaleBar() %>%
        addCircleMarkers(data = metadata_shiny_filt(), lat =  ~Location_latitude, lng = ~Location_longitude,
                         radius = 5, popup = ~as.character(cntnt),
                         color = ~pal(Species),
                         fill = ~pal(Species),
                         stroke = FALSE, fillOpacity = 0.8)%>%
        addLegend(pal=pal, values=metadata_shiny_filt()$Species, opacity=1, na.label = "Not Available")%>%
        addEasyButton(easyButton(
          icon="fa-crosshairs", title="ME",
          onClick=JS("function(btn, map){ map.locate({setView: true}); }")))
      })

  #create a data object to display data

  output$data <-DT::renderDataTable(datatable(
      subset(metadata_shiny_filt(), select = -c(cntnt)),
      extensions = "Scroller",
      filter = c("top"),
      # options=list(columnDefs = list(list(visible=FALSE, targets=c(29,30,31,32,33))))
  ))

  output$data_extra <-DT::renderDataTable(datatable(
      subset(metadata_shiny_filt(), select = -c(cntnt)),
      extensions = "Scroller",
      filter = c("top"),
      options=list(columnDefs = list(list(visible=FALSE, targets=3:20)))
  ))


}

ui = navbarPage("Metadata map", id="main",
           tabPanel("Map",
                sidebarLayout(
                       sidebarPanel(
                         dateRangeInput('dateRange',
                                         label = 'Filter by date',
                                         start = as.Date('2000-01-01') , end = as.Date(Sys.Date())
                          )
                        ),
                mainPanel(leafletOutput("metamap", height=1000))
              )
            ),
           tabPanel("Data", DT::dataTableOutput("data")),
           tabPanel("Notes", DT::dataTableOutput("data_extra"))
         )

shinyApp(ui = ui, server = server)
