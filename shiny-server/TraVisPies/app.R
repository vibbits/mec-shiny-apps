# Description ---------------------------------------------------
###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: UI for creating and saving tracer metabolomics pie chart plot 
# images proposed in the manuscript TraVis Pies: An intuitive visualization 
# approach for stable-isotope tracer metabolomics.
# Applies desired settings to a compound of
# choice to showcase their effects, then allows to user to apply these to all
# compounds in the database once they are satisfied, saving the resulting
# figures to a folder of choice

###Input
# Either the user clicks on the create new input file button, or the user 
# chooses a file to upload. The former has to be done when starting from new
# data, and will take 3 input .csv files respectively containing metadata, 
# abundance data and fractional contribution data and merge them to a single
# file to be used in TraVis pies, which can be saved to skip this step next
# time the same data is used by choosing the upload option.

###Output
# figures generated with desired settings from input data,
# saved locally to provided folder

#Developer settings
#Do not change below setting unless you know what you are doing!

#indicate whether local version of code should be used, otherwise uses web version,
#reason: change how output files should be saved (locally or with download)
local_version<-F      

#closes the R session when the application window is closed. Only required
#when running as a standalone app, could be annoying or troublesome otherwise
closeR_with_app<- F   


####TODO someday when time
#add tabs to switch between panels by turning off invisible mode, 
#double check logic to avoid errors:
#-switching in and out of inputcleaner should be fine
#-output tab should only be available if an output is generated in the visual
#module

# Functions and libraries ---------------------------------------------------------------
#libraries for shiny application
library(shiny)
library(shinyFeedback)   #for error messages on box
library(shinyjs)      #for show/hide button
library(shinyFiles)   #for saving locally from shiny application
library(DT)      #to display interactive data tables in shiny
library(colourpicker) #to allow colours to be chosen

#libraries for functions used
library(here)    #to make r source from file location instead of magic stuff
library(vroom)   #for error messages on box
library(forcats)      #for factor manipulation
library(dplyr)        #for faster.easier manipulation of data
library(tibble)       #for manipulating tibbles
library(readr)        #for writing .csv file of merged output
library(tidyr)        #for restructuring data tibbles
library(ggplot2)      #for generating the pie chart plots

#load functions to support the app 
source(here::here("Functions and modules/TraVis_Pies_functions.R"))


if (local_version) {
  library(extrafont)    #for using other fonts
  check_install_fonts()
} else {
  
  # specific for MEC linux server, install add extrafontdb folder to 
  #library paths
  .libPaths(c('r-lib', .libPaths()))
  library(extrafontdb) 
  library(extrafont)    #for using other fonts
  check_install_fonts()
} 

#Some operations only for windows systems
if (.Platform$OS.type=="windows") {
  #required for zipping output files in web tool with rtools zip tool:
  #check if rtools path is set in .Renviron file yet (if that file exists at all),
  #set it and reread Renviron file if not
  if(!file.exists("~/.Renviron")) {
    write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron",
          append = TRUE)
    readRenviron("~/.Renviron")
  } else {
    if (!grepl('PATH=${RTOOLS40_HOME}\\usr\\bin;${PATH}',
               read.delim("~/.Renviron"),fixed=T)) {
      write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron",
            append = TRUE)
    }
    readRenviron("~/.Renviron")
  }
}

#load the modules
source(here::here("Functions and modules/TraVis_InputCleaner.R"))
source(here::here("Functions and modules/TraVis_input.R"))
source(here::here("Functions and modules/TraVis_Visualisation.R"))
source(here::here("Functions and modules/TraVis_output.R"))

#Wrapper main app-------------------------------------------------------------------
#show inputmodule and vizmodule, unless the create new input button was pressed
#in the inputmodule, then show inputclean
ui <- fluidPage(
  #App explanation text
  fluidRow(
    column(
      4,titlePanel(h1("Travis Pies v1.31"))
    ),
    column(
      3,downloadButton(outputId= "down_manual",
                       label = "User manual",
                       style = "margin-top: 25px;margin-bottom: 25px;")
    )
  ),
  htmlOutput("text"),
  
  tabsetPanel(
    id = "wizard",
    type = "hidden",
    selected =  "Visualisation",
    tabPanel("InputClean", 
             travis_cleaner_ui("cleaner")
    ),
    tabPanel("Visualisation", 
             travis_input_ui("input"),
             travis_viz_ui("visual"),
             
             #click when obtained desired visual settings to apply
             fluidRow(column(5,
                             disabled(
                               actionButton("save_plots",
                                            "Save plots with these settings")))
             ),
             
    ),
    tabPanel("Output",
             if(local_version) {
               travis_outputlocal_ui("output")
             } else {
               travis_outputweb_ui("output")
             },
             HTML("<br/>"),
             actionButton("return_settings","Return to settings")
    )
  )
)

server <- function(input, output, session) {
  
  #Write text to put on top as explanation, need to use server output to be able 
  #to write multiple lines
  output$text <- renderText({
    paste0("<b>Create and save tracer metabolomics pie chart plots.</b><br/> ",
           "Allows creation of a large set tracer metabolomics pie charts ",
           "based on a compatible input file, that can be created or ",
           "uploaded.<br/> Based on this data, the visualisation will be ",
           "shown for a compound of choice, and the user can change the ",
           "settings until they are satisfied. <br/> Then these figures can be",
           " generated for all compounds in the input data and saved in a ",
           "folder of choice <br/>")
  })
  
  #Download manual
  output$down_manual = downloadHandler(
    filename = "TraVis_Pies_manual.docx",
    content = function(file){
      file.copy(here::here("TraVis_Pies_manual.docx"),
                file)
    }
  )
  
  #prepare reactive values
  input_v <- travis_input_server("input")
  create_v <- travis_cleaner_server("cleaner",local_version)
  v <- reactiveValues(active_tb=data.frame(NA))
  visual_v<-travis_viz_server("visual",reactive(v$active_tb))
  
  #Switch to inputclean module if create button is pressed in input module
  #switch back once inputclean is marked as finished
  observeEvent(input_v$create_new,{
    if (input_v$create_new) {
      updateTabsetPanel(inputId = "wizard", 
                        selected = "InputClean")
    }
  })
  
  ####setting active table 
  #below are two events that cause the active table to be set to whichever
  #of them occured the last
  
  #if the continue button is clicked in the clean input UI, set create_new
  #to false again, set the active data to the merged input data
  #switch back to the main UI and set create_v$finished back to false
  #so the user can modify the input by clicking create new again
  #NOTE: changing a reactive value also changes the reactive value that used
  #to be its input, so the value in the module is also reset!
  observeEvent(create_v$finished,{
    if (create_v$finished) {
      input_v$create_new<-F
      v$active_tb <- create_v$merged_tb
      updateTabsetPanel(inputId = "wizard", selected = "Visualisation")
      create_v$finished<-F
    }
  })
  
  #set active table to the value of the uploaded table
  observeEvent(input_v$checked_upload_tb,{
    v$active_tb <- input_v$checked_upload_tb
  })
  
  
  ### visual module output
  #get the output of the visual module once there is a valid input tibble
  observeEvent(v$active_tb,{
    req(!is.na(v$active_tb[1,1]))
    visual_v<-travis_viz_server("visual",reactive(v$active_tb))
  })
  
  #prepares output and enables continuing if all settings valid
  observe({
    disable("save_plots")
    req(visual_v$ready)
    
    #changes download version depending on whether code should use 
    #structure for saving locally on server or through download
    if(local_version) {
      travis_outputlocal_server("output",v_settings = visual_v,
                                tb = reactive(v$active_tb))
    } else {
      travis_outputweb_server("output",v_settings = visual_v,
                              tb = reactive(v$active_tb))
    }
    
    enable("save_plots")
  })
  
  ###output module interaction
  #if the continue button is clicked in the visualisation UI, 
  #switch to the output UI and load server with this input
  observeEvent(input$save_plots,{
    req(!is.na(visual_v$norm))
    updateTabsetPanel(inputId = "wizard", selected = "Output")
  })
  
  #return to visual tab if return button is pressed in output tab
  observeEvent(input$return_settings,{
    updateTabsetPanel(inputId = "wizard", selected = "Visualisation")
  })
  
  # close the R session when app window closes if local
  if (closeR_with_app) {
    session$onSessionEnded(function() {
      stopApp()
      q("no")
    })
  }
} 

shinyApp(ui, server)

