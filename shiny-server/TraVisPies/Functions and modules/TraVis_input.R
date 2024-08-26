# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: Input module for the TraVis Pies ui. Give choice of creating a new
# TraVis input file from three raw input files by opening the TraVis input 
# cleaner module, or simply upload an already saved TraVis input file after
# previously cleaning the raw input. In case of the latter, this file is 
# checked to see if it conforms to the required input format, then formatted
# further to be identical to a freshly created file

###Input
# Either the user clicks on the create new input file button, or the user 
# chooses a file to upload. By writing the output value checking if the create
# button was clicked to a reactive value in the parent application
# and setting that to false afterwards there, the create file can be repeated

###Output
# The module outputs a reactiveValues v that contains two outputs
# v$create_new is False by default but become TRUE when the user clicks on
# the create new input button
# v$checked_upload_tb is data.frame(NA) by default but is replaced by
# a tibble a file with the correct format was uploaded.


# Functions and libraries ---------------------------------------------------------------
library(shiny)
library(vroom)   #for error messages on box
library(forcats)      #for factor manipulation
library(dplyr)        #for faster.easier manipulation of data as tibbles
library(here)    #to find example datafiles


#helper function for input module to check uploaded merged file for validity
check_uploaded_table<-function(tb){
  errormessage<-paste0("Uploaded file contents do not match TraVis Input ",
                       "file format. Please select the right file or ",
                       "create a new TraVis Input file and do not manually ",
                       "modify this file")
  
  #check required columns Sample and datatype
  if (!all(c("Sample","datatype") %in% colnames(tb))) 
    validate(errormessage)
  
  
  #check two columns before datatype (Sample and factor)
  datatype_index<-which(colnames(tb)=="datatype")
  if (!datatype_index>2) validate(errormessage)
  
  
  #check if datatype contains Fraccont and Abundance, and only allow but don't 
  #require NormAbund as a third possibility
  if (! all(c("Abund","FracCont") %in% tb$datatype) |
      ! all(tb$datatype %in% c("Abund","FracCont","NormAbund","Isotopologues")))
    validate(errormessage)
  
  return("")
}

#helper function for input module to format uploaded merged file
format_uploaded_table<-function(tb){
  
  #Make sure table has right formatting
  datatype_index<-which(colnames(tb)=="datatype")
  tb<-mutate(tb,across(1:(datatype_index-1), as.character))
  
  return(tb)
}


#Input module UI and server-------------------------------------------------------------------
travis_input_ui<- function(id) {
  fluidPage(
    #set color for error messages
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
      color: red;
      }
      "))
    ),
    
    #Input selection: upload existing merged file or create new
    div(style="display:inline-block",h3(paste0("Choose input method"))),
    div(style="display:inline-block; ",
        actionButton(NS(id,"ui_showhide"), label = "show / hide")),
    
    htmltools::div(
      id=NS(id,"Input_options"),
      fluidRow(
        column(
          4,
          actionButton(NS(id,"create_input"),
                         label = "Create TraVis input file from curated data",
                         style = "margin-top: 10px;margin-bottom: 10px;")
        ),
      ),  
      fluidRow(
        column(
          4,
          fileInput(NS(id,"upload_file"),
                    "Upload previously created standardized TraVis input file",
                    accept = ".csv",width = "100%")
        ),
        column(
          3,
          actionButton(NS(id,"demo_input"),
                       label = "Use demo input",
                       style = "margin-top: 25px;margin-bottom: 25px;"),
        ),
      ),  
      textOutput(NS(id,"input_valid")),
      verbatimTextOutput("check")
      
      # tableOutput(NS(id,"inputtable"))
    )
  )
}


travis_input_server <- function(id) {
  moduleServer(id,function(input, output, session) {
    #Reactivevalues with output
    #exist_tb is the input uploaded in the ui of this module
    v<-reactiveValues(create_new=F,
                      checked_upload_tb=data.frame(NA))
    
    ###Show/hide input options
    observeEvent(input$ui_showhide, {
      if(input$ui_showhide %% 2 == 0){
        shinyjs::show(id = "Input_options")
      }else{
        shinyjs::hide(id = "Input_options")
      }
    })
    
    #Fileloading input data
    upload_tb<-reactive({
      req(input$upload_file$datapath)
      vroom::vroom(input$upload_file$datapath,delim = ",",
                   show_col_types = FALSE)
    })
    
    #Change reactive value create_new to true if creation of new file is
    #requested
    observeEvent(input$create_input,{
      v$create_new<-T
    })
    
    #Use demo data if requested
    observeEvent(input$demo_input,{
      v$checked_upload_tb<-format_uploaded_table(
        read_csv(here::here(
          paste0("Example_data/Standardized input/Input_Example_standardized ",
          "w isotopologues.csv")))
      )
    })
  
    #Check if input file is valid
    output$input_valid<-renderText({
      req(nrow(upload_tb())>0)
      
      #perform check, forward validate if fails
      validstring<-check_uploaded_table(upload_tb())
      
      #if check succeeded update reactivevalue, save correctly formatted
      #data table correctly
      v$checked_upload_tb<-format_uploaded_table(upload_tb())
      
      return(validstring)
    })
    
    
    
    return(v)
  })
}