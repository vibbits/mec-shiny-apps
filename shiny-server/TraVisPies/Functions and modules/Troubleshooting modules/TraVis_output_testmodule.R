# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: Travis Pies shiny module that takes in desired settings for pie
# chart plots and applies them on a provided curated dataset to 
# generate figures with them.

###Input
# desired settings for pie chart plots and Travis curated tibble

###Output
# figures saved locally to provided folder

# Functions and libraries ---------------------------------------------------------------

#libraries for UI
library(shiny)
library(shinyFeedback)   #for error messages on box
library(shinyjs)      #for disable button
library(shinyFiles)   #for saving locally

#library for calling helper functions from r script
library(here)    #to make r source from file location instead of magic stuff
source(here::here("Functions and modules/TraVis_Pies_functions_beta.R"))




#Shiny Modules-------------------------------------------------------------------
travis_output_ui <- function(id) {
  fluidPage(
    #load shiny package elements
    
    #Button to start process
    fluidRow(
      column(2,
             actionButton(NS(id, "save_plots"),
                                   label = "Save plots to directory")
      )
    )
  )
}

travis_output_server <- function(id,v_settings,tb) {
  moduleServer(id,function(input, output, session) {

    #Write text to put on top as explanation, need to use server output to be able 
    #to write multiple lines
    output$text <- renderText({
      paste0("<b>Test action button sheninigans.</b><br/><br/>")

    })
    
    
    
    ## do the action
    observeEvent(input$save_plots,{
    print(Sys.time())
    })
  })
}

#Wrap modules in caller function to test-------------------------------------------------------------------
#to make correspond to input
example_tb<-tibble(Sample = c("S1","S2","S3","S4","S1","S2","S3","S4",
                              "S1","S2","S3","S4"),
                   Cohort=c("coh1","coh1","coh2","coh3",
                            "coh1","coh1","coh2","coh3",
                            "coh1","coh1","coh2","coh3"),
                   Mugwort2=c("Mugw1","Mugw1","Mugw2","Mugw2",
                              "Mugw1","Mugw1","Mugw2","Mugw2",
                              "Mugw1","Mugw1","Mugw2","Mugw2"),
                   datatype=c("Abund","Abund","Abund","Abund",
                              "FracCont","FracCont","FracCont","FracCont",
                              "NormAbund","NormAbund","NormAbund","NormAbund"),
                   Result = c(40005,45858,7000000,5000000,
                              .5,.453,.32,1,
                              400005,458058,700000,500000),
                   Result2 = c(40005,4585824,7552123,50000,
                               .8,.75,.32,0.3,
                               400005,45805824,755123,5000))
norm<-T

#example without normalized abundance
# example_tb<-tibble(Sample = c("S1","S2","S3","S4","S1","S2","S3","S4",
#                               "S1","S2","S3","S4"),
#                    Cohort=c("coh1","coh1","coh2","coh3",
#                             "coh1","coh1","coh2","coh3",
#                             "coh1","coh1","coh2","coh3"),
#                    Mugwort2=c("Mugw1","Mugw1","Mugw2","Mugw2",
#                               "Mugw1","Mugw1","Mugw2","Mugw2",
#                               "Mugw1","Mugw1","Mugw2","Mugw2"),
#                    datatype=c("Abund","Abund","Abund","Abund",
#                               "FracCont","FracCont","FracCont","FracCont",
#                               "Abund","Abund","Abund","Abund"),
#                    Result = c(40005,45858,7000000,5000000,
#                               .5,.453,.32,1,
#                               400005,458058,700000,500000),
#                    Result2 = c(40005,4585824,7552123,50000,
#                                .8,.75,.32,0.3,
#                                400005,45805824,755123,5000))
#norm<-F

datatype_index<-which(colnames(example_tb)=="datatype")

example_tb<-mutate(example_tb,Sample= as.character(Sample),
                   across(2:(datatype_index-1), as.factor))
factorder_input<-unique(as.character(pull(example_tb[,2])))
example_tb[,2]<-fct_relevel(pull(example_tb[,2]),factorder_input)

example_settings<-list(finish=F,fact_name = "Cohort",
                       fact_order=factorder_input,norm = norm,
                       label_decimals = 0,
                       min_lab_dist = 0.42,
                       percent_add = F,
                       FC_position = "center",
                       col_labeling = c("#bfbfbf","#ffd966"),
                       circlelinetypes = c(1,1,1,1),
                       circlelinecolor = "gray",
                       maxcol_facet= 4,
                       include_name = T,show_P=F,alpha=0.7,
                       font="sans",otherfontsize = 16,
                       legendtitlesize =16,
                       cohortsize = 18,
                       include_legend = T,
                       show_P=T)


travis_output_app<- function() {
  ui <- fluidPage(
    travis_output_ui("output")
  )
  

  server <- function(input, output, session) {
    travis_output_server("output",v_settings = 
                           do.call("reactiveValues",example_settings),
                        tb = reactive(example_tb))
  }
  
  shinyApp(ui, server)  
}

travis_output_app()


