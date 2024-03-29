# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: Shiny module that cleans inputdata for TraVis Pies Shiny v1.0, 
# entering dummy columns, removing superficial columns and standardizing naming
# and order to ensure that the resulting merged tibble can be reliably used by
# the other modules of the Travis Pies shiny app.

###Input
# UI resuests 3 input files: metadata, raw abundance and fractional
# contribution. Several other tidyverse packages are used to get the data in the 
# right format to generate the required figures.

###Output
# A curated TraVis input tibble that should always be compatible with downstream
# TraVis Pies apps. If input data is not internally consistent, warnings or 
# errors willguide the user to either make it so in the app, or to modify the 
# input data and try again. Contains Sample cohort and compounds columns, and 
# for each compound in each sample an abundance and fractional contribution. If 
# a  normalisation variable was specified, in addition there will also be a 
# normalized abundance for each compound in each sample

###Input file tips
# The samples considered will be based on the metadata file provided. Samples 
# in other files that are not in the metadatafile will not be included in 
# the final tibble.


# Functions and libraries ---------------------------------------------------------------
#libraries for UI
library(shiny)
library(shinyFeedback)   #for error messages on box
library(shinyjs)      #for show/hide button
library(shinyFiles)   #for saving locally
library(here)    #to make r source from file location instead of magic stuff
library(DT)      #to display interactive data tables in shiny

#library for calling helper functions from r script
source(here::here("Functions and modules/TraVis_Pies_functions.R"))

#updates inputbox choices with choice list, and add an option to pick none
#specified in noPick, then selects the default option based other input:
#select by preference the first choice that contains a preference string
#if no match, will pick the first available choice by default unless 
#prefOrNo is set to TRUE in which case it will pick the noPick option
upd_selInp_guess_default<-function(inputId,choices,prefString,noPick="None",
                                   prefOrNo=F){
  allchoices <- c(choices,noPick)
  choice_selected_index<-grep(prefString,tolower(allchoices))
  if (length(choice_selected_index) == 0) {
    choice_selected <- if_else(prefOrNo,allchoices[length(allchoices)],
                               allchoices[1])
  } else {
    choice_selected <- allchoices[choice_selected_index[1]]
  }
  
  updateSelectInput(inputId = inputId, choices = allchoices,
                    selected=choice_selected) 
}

#Explanation texts-------------------------------------------------------------------
text_exp_3file_FC<-
  paste0("See Example 3-file input FC",
         " above. All files contain an identically named column",
         " with sample names.<br/>Metadata optionally contains cohorts and ",
         "normalisation columns. The other two contain only the sample ",
         "name column and their respective data columns named after the ",
         "compounds which they represent.<br/>Fractional contribution ",
         "data can be provided as a fraction between slightly less than 0 ",
         "and 1 or a percentage between slightly less than 0% and 100%. ",
         "12C controls should be in a separate cohort.",
         "<br/> Normalisation happens by dividing abundances by the ",
         "normalisation value provided.")

text_exp_3file_iso<-
  paste0("See Example 3-file input isotopologues",
         " above. All files contain an identically named column",
         " with sample names.<br/>Metadata optionally contains cohorts and ",
         "normalisation columns. The other two contain only the sample ",
         "name column and their respective data columns named after the ",
         "compounds which they represent.<br/>Isotopologue ",
         "data can be provided as a fraction or a percentage, with negative ",
         "values left as is or set to 0. ",
         "12C controls should be in a separate cohort.",
         "<br/> Normalisation happens by dividing abundances by the ",
         "normalisation value provided.")

text_exp_2file_et<-
  paste0("See Example 2-file Escher-Trace input above. ",
         "Metadata contains sample names and optionally cohorts ",
         "and normalisation columns. 12C controls should be in a separate ",
         "cohort, if provided at all.",
         "<br/>The isotopologue input file is the same as ",
         "the corrected isotopologue input file for Escher trace. The ", 
         "first row of the CSV must include the following headings in ",
         "order: Metabolite, Fragment followed by the sample names.",
         "<br/> Sample names in the columns of this file should correspond to ",
         "sample names in the rows of the metadata file. ",
         "More details on how to organise the input in the manual ",
         "that can be downloaded by clicking the button at the top.",
         "<br/> Normalisation happens by dividing abundances by the ",
         "normalisation value provided.")
  

#Local output module-------------------------------------------------------------------
tibble_localoutput_ui <- function(id) {
  strong(paste0("Optional: save the merged table. Copypaste the folderpath ",
                "in textbox below and save to specified folder for faster saving"),
         fluidRow(
           column(4,
                  textInput(NS(id,"userdir"),"",)
           ),
           column(3,
                  shinySaveButton(NS(id,"saveuser"), "Save to specified directory",
                                  "Save file as ...", filetype=list(csv="csv"),
                                  style = "margin-top: 20px;margin-bottom: 25px;",
                                  filename = "merged data")
           ),
           column(3,
                  shinySaveButton(NS(id,"saveAll"), 
                                  "Save to any directory",
                                  "Save file as ...", filetype=list(csv="csv"),
                                  style = "margin-top: 20px;margin-bottom: 25px;",
                                  filename = "merged data")
           )),
  )
}

tibble_localoutput_server <- function(id,tb) {
  moduleServer(id,function(input, output, session) {
    
    ## Allow user to save to all local folders 
    observe({
      volumes <- getVolumes() # this makes the directory at the base of your computer.
      shinyFileSave(input, "saveAll", roots=volumes(), session=session)
      fileinfo <- parseSavePath(volumes, input$saveAll)
      if (nrow(fileinfo) > 0) {
        write_csv(tb, as.character(fileinfo$datapath))
      }
    })
    
    ## Allow user to save to specified folder
    observe({
      usersavepath<-trimws(gsub("\\\\", "/",input$userdir)) #get valid path
      volumes<- c(user=usersavepath)
      shinyFileSave(input, "saveuser", roots=volumes, session=session)
      fileinfo <- parseSavePath(volumes, input$saveuser)
      if (nrow(fileinfo) > 0) {
        write_csv(tb, as.character(fileinfo$datapath))
      }
    })
    
    # # close the R session when Chrome closes
    # session$onSessionEnded(function() { 
    #   stopApp()
    #   q("no") 
    # })
  })
}

#Web output module-------------------------------------------------------------------
tibble_weboutput_ui <- function(id) {
  fluidRow(
    column(2,
           downloadButton(NS(id, "download_merged"),
                          label = "Download merged data")
    )
  )
}

tibble_weboutput_server <- function(id,tb) {
  moduleServer(id,function(input, output, session) {
    
    # Create and download the plots
    output$download_merged = downloadHandler(
      filename = 'merged data.csv',
      content = function(file){
        # Set temporary working directory
        owd <- setwd( tempdir())
        on.exit( {
          setwd( owd)
        })
        write_csv(tb,file)
        
      }
    )
    
  })
}



#Cleaner module-------------------------------------------------------------------
travis_cleaner_ui <- function(id) {
  fluidPage(
    #load shiny package elements
    shinyFeedback::useShinyFeedback(),
    useShinyjs(),    ## IMPORTANT: so shiny knows to use the shinyjs library
    
    #set color for error messages
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
      color: red;
      }
      "))
    ),
    
    #App explanation text
    fluidRow(
      column(
        4, titlePanel(h2("Input Standardizer"))
      )
    ),
    fluidRow(h3(),
      column(
        7,radioButtons(inputId = (NS(id,"input_type")),
                     label = "Choose input type",
                     choiceNames = c(
                       paste0("Upload 3 files with metadata, abundance data ",
                       "and fractional contribution data"),
                       paste0("Upload 3 files with metadata, abundance data ",
                       "and isotopologue data"),
                       paste0("Upload 2 files with metadata and abundance/",
                              "isotopologue data (Escher-Trace input)")
                       ),
                     choiceValues = c("3file_FC","3file_iso","2file"),
                     width="100%"
                    )
      ),
      column(
        3,downloadButton(outputId = (NS(id,"down_3example_FC")),
                         label = "Example 3-file input FC",
                         style = "margin-top: 10px;margin-bottom: 0px;"),
        downloadButton(outputId = (NS(id,"down_3example_iso")),
                       label = "Example 3-file input isotopologues",
                       style = "margin-top: 0px;margin-bottom: 0px;"),
        downloadButton(outputId = (NS(id,"down_2example")),
                       label = "Example 2-file input Escher-Trace",
                       style = "margin-top: 0px;margin-bottom: 5px;")
      )
    ),
    
    div(style="display:inline-block; ",
        actionButton(NS(id,"input_showhide_explanation"), 
                     label = "Show / hide explanation")),
    
    shinyjs::hidden(
      htmltools::div(
        id=NS(id,"explanation"),
        htmlOutput(NS(id,"explanation_text"))
      )
    ),
    
    h2(paste0("Upload all input files below to continue, uploaded",
              " data can be viewed below")
    ),
    
    tabsetPanel(
      id = (NS(id,"input_panel")),
      type = "hidden",
      selected =  "3file_FC_input",
      tabPanel("3file_FC_input", 
               #input files
               fluidRow(
                 tags$style(HTML("div.input-group {margin-bottom: -30px;}")),
                 column(
                   3,fileInput(NS(id,"meta3_FC_file"), "Metadata",accept = ".csv")),
                 column(
                   3,fileInput(NS(id,"abund_FC_file"), "Abundance data",accept = ".csv")),
                 column(
                   3,fileInput(NS(id,"frac_file"), "Fractional contribution data",accept = ".csv"))
               ),
      ),
      tabPanel("3file_iso_input", 
               #input files
               fluidRow(
                 tags$style(HTML("div.input-group {margin-bottom: -30px;}")),
                 column(
                   3,fileInput(NS(id,"meta3_iso_file"), "Metadata",
                               accept = ".csv")),
                 column(
                   3,fileInput(NS(id,"abund_iso_file"), "Abundance data",
                               accept = ".csv")),
                 column(
                   4,fileInput(NS(id,"iso_col_file"),"Corrected isotopologue data",
                               accept = ".csv",width = "100%"))
                 )
      ),
      tabPanel("2file_input", 
               #input files
               fluidRow(
                 tags$style(HTML("div.input-group {margin-bottom: -30px;}")),
                 column(
                   3,fileInput(NS(id,"meta2_file"), "Metadata",accept = ".csv")),
                 column(
                   6,fileInput(NS(id,"iso_et_file"), 
                               paste0("Corrected isotopologue and abundance",
                               "(Escher-Trace input)"),width = "100%",
                               accept = ".csv"))
                 )
      )
    ),

    HTML(paste0("<b>Samples present in metadata required in the other files, ",
             "samples not present in metadata will be removed</b><br/><br/>")),
    textOutput(NS(id,"input_valid")),
    
    #Finetune selection of metadata columns, select compounds
    conditionalPanel(
      "output.input_valid === ''",
      ns = NS(id),
      h3("Adapt the selected input settings if necessary or desired"),
      fluidRow(
        column(4,
               #for some reason shinyfeedback only works well if selectize is off
               selectInput(NS(id,"sample_column"), label= "Select sample column",
                           choices = NULL,selectize = FALSE)
        ),
        column(7,
               selectInput(NS(id,"norm_column"),
                           label= paste0("Select normalisation column. ",
                                         "Normalisation can still be ",
                                         "toggled further."),
                           choices = NULL,
                           width="100%")
        )
      ),
      fluidRow(
        column(3,
               selectInput(NS(id,"cohort_column"), 
                           label= "Select cohort column if present",
                           choices = NULL)
        )
      ),
      strong("Select compounds to be analysed, by default all in abundance file"),
      actionButton(NS(id,"comp_showhide"), label = "show / hide"),
      fluidRow(
        column(10,
               selectInput(NS(id,"compounds"),label=NULL,choices= NULL, multiple = T,
                           width="100%")
        )
      ),
      htmlOutput(NS(id,"sample_valid")),
      
      
      #buttons to generate merged file and continue
      h3("Generate merged table and continue"),
      
      fluidRow(
        column(2,
               disabled(actionButton(NS(id, "merge"),
                                     label = "Merge inputdata"))
        ),
        column(3,
               disabled(actionButton(NS(id, "finish"), 
                            label = "Continue with this data")
        )),
        
      ),
      
        tabsetPanel(
          id = NS(id,"downoptions"),
          type = "hidden",
          selected =  "None",
          tabPanel("None",
          ),
          tabPanel("local",
                   tibble_localoutput_ui(NS(id,"localout"))
          ),
          tabPanel("web",
                   tibble_weboutput_ui(NS(id,"webout"))
          )
        )
    ),  
      
    
    #table output
    conditionalPanel(
      "output.any_upload == true",
      ns = NS(id),
      h3("Select a dataset to view"),
      fluidRow(column(
        12,selectInput(NS(id,"view_table"), label= NULL,choices = NULL)
      )
      ),
      fluidRow(column(
        12,div(DT::dataTableOutput(NS(id,"table")))
      )
      )
    ),
    
    
    
    
    #other page descriptors
    title = "Travis Pies Input check v1.0"
  )
}

travis_cleaner_server <- function(id,local_version=T) {
  moduleServer(id,function(input, output, session) {
    
    
    #reactive value prepared for storing merged_tb, initialized as NA
    #for tests until a table is generated
    v<-reactiveValues(merged_tb=data.frame(NA),
                      finished=F)
    
    #Download example input files
    output$down_3example_FC = downloadHandler(
      filename = 'Input_meta_abund_FC_Example.zip',
      content = function(file){
        file.copy(here::here("Example_data/Input_meta_abund_FC_Example.zip"),
                  file)
      }
    )
    output$down_3example_iso = downloadHandler(
      filename = 'Input_meta_abund_iso_Example.zip',
      content = function(file){
        file.copy(here::here("Example_data/Input_meta_abund_iso_Example.zip"),
                  file)
      }
    )
    output$down_2example = downloadHandler(
      filename = 'Input_meta_Escher-Trace_Example.zip',
      content = function(file){
        file.copy(here::here("Example_data/Input_meta_Escher-Trace_Example.zip"),
                  file)
      }
    )
    
    #return text depending on input option selected
    output$explanation_text <- renderText({
      if (input$input_type=="3file_FC") {
       return(text_exp_3file_FC)
      }
      
      if (input$input_type=="3file_iso") {
        return(text_exp_3file_iso)
        
      }
      
      if (input$input_type=="2file") {
        return(text_exp_2file_et)
        
      }
      
    })
    
    ###Show/hide texts options
    observeEvent(input$input_showhide_explanation, {
      if(input$input_showhide_explanation %% 2 == 1){
        shinyjs::show(id = "explanation")
      }else{
        shinyjs::hide(id = "explanation")
      }
    })
    
    #Switch to correct input file panel when type selection radiobutton is
    #pressed
    observeEvent(input$input_type,{
      if (input$input_type=="3file_FC") {
        updateTabsetPanel(inputId = "input_panel",
                          selected = "3file_FC_input")
      }
      if (input$input_type=="3file_iso") {
        updateTabsetPanel(inputId = "input_panel",
                          selected = "3file_iso_input")
      }
      if (input$input_type=="2file") {
        updateTabsetPanel(inputId = "input_panel",
                          selected = "2file_input")
      }
    })
    
    ###Loading input data
    #prepare data reactive variable to change depending on input
    meta_tb<-reactiveVal(tibble(NA))
    abund_tb<-reactiveVal(tibble(NA))
    frac_tb<-reactiveVal(tibble(NA))
    iso_col_tb<-reactiveVal(tibble(NA))
    iso_et_tb<-reactiveVal(tibble(NA))
    iso_tb<-reactiveVal(tibble(NA))
    
    

    #set to tibbles to newly uploaded corresponding data
    observeEvent(input$meta3_FC_file,{
      meta_tb(read_csv_clean(input$meta3_FC_file$datapath,remove_empty = T))
      }
    )
    
    observeEvent(input$meta3_iso_file,{
      meta_tb(read_csv_clean(input$meta3_iso_file$datapath,remove_empty = T))
    })

    observeEvent(input$meta2_file,{
      meta_tb(read_csv_clean(input$meta2_file$datapath,remove_empty = T))
    })
    
    observeEvent(input$abund_FC_file,{
      abund_tb(read_csv_clean(input$abund_FC_file$datapath,remove_empty = T))
    })
    
    observeEvent(input$abund_iso_file,{
      abund_tb(read_csv_clean(input$abund_iso_file$datapath,
                              remove_empty = T))
    })
    
    observeEvent(input$frac_file,{
      frac_tb(read_csv_clean(input$frac_file$datapath,remove_empty = T,
                             perc_to_num = T))
    })
    
    observeEvent(input$iso_col_file,{
      #Load in slightly cleaned iso data for checks
      iso_col_tb(read_csv_clean(input$iso_col_file$datapath,remove_empty = T,
             perc_to_num = T))

      #Transform columnwise iso data to rowwise iso data for easier further 
      #calculations. Assumes iso label separator suffix is "_"
      iso_tb(extract_col_isotopologues(iso_col_tb(),iso_suffix_sep = "_"))
      print(iso_tb)
      
      #calculate fractional contribution from isotopologue data
      frac_tb(calculate_FC(iso_tb()))
    })
    
    observeEvent(input$iso_et_file,{
      iso_et_tb(read_csv_clean(input$iso_et_file$datapath,remove_empty = T))
      
      #Calculate abundances from escher trace input
      abund_tb(extract_et_abund(iso_et_tb()))
      
      #calculate isotopologue tibble from escher trace input
      iso_tb(extract_et_isotopologues(iso_et_tb()))
      
      #calculate fractional contribution from isotopologue data
      frac_tb(calculate_FC(iso_tb()))
    })

    
    #set to currently uploaded variables in new input option when switching
    observeEvent(input$input_type,{
      if (input$input_type=="3file_FC") {
        #only set if files are uploaded
        req(length(input$meta3_FC_file)*
              length(input$abund_FC_file)* length(input$frac_file)>0)
        meta_tb(read_csv_clean(input$meta3_FC_file$datapath,remove_empty = T))
        abund_tb(read_csv_clean(input$abund_FC_file$datapath,remove_empty = T))
        frac_tb(read_csv_clean(input$frac_file$datapath,remove_empty = T))
        iso_tb(tibble(NA))
      }
      
      if (input$input_type=="3file_iso") {
        #only set if files are uploaded
        req(length(input$meta3_iso_file)*
              length(input$abund_iso_file)* length(input$iso_col_file)>0)
        meta_tb(read_csv_clean(input$meta3_iso_file$datapath,
                               remove_empty = T))
        abund_tb(read_csv_clean(input$abund_iso_file$datapath,
                                remove_empty = T))
        iso_col_tb(read_csv_clean(input$iso_col_file$datapath,remove_empty = T,
                                  perc_to_num = T))
        frac_tb(tibble(NA))
        
        #Only if columnwise iso data loaded
        #Transform columnwise iso data to rowwise iso data for easier further 
        #calculations. Assumes iso label separator suffix is "_"
        req(nrow(iso_col_tb())*ncol(iso_col_tb())>0)
        iso_tb(extract_col_isotopologues(iso_col_tb(),iso_suffix_sep = "_"))
        
        #calculate fractional contribution from isotopologue data
        frac_tb(calculate_FC(iso_tb()))
      }
      
      if (input$input_type=="2file") {
        req(length(input$meta2_file)*
              length(input$iso_et_file)>0)
        req(nrow(meta_tb())*ncol(meta_tb())>0)
        meta_tb(read_csv_clean(input$meta2_file$datapath,remove_empty = T))
        iso_et_tb(read_csv_clean(input$iso_et_file$datapath,remove_empty = T,
                                  perc_to_num = T))
        abund_tb(tibble(NA))
        frac_tb(tibble(NA))
        
        #Only if eschter trace abundance + iso data loaded
        #Calculate abundances from escher trace input
        #calculate isotopologue tibble from escher trace input
        #calculate fractional contribution from isotopologue data
        req(nrow(iso_et_tb())*ncol(iso_et_tb())>0)
        abund_tb(extract_et_abund(iso_et_tb()))
        iso_tb(extract_et_isotopologues(iso_et_tb()))
        frac_tb(calculate_FC(iso_tb()))
      }
    })
    
    #add a check to see uploaded data is valid
    output$input_valid <- renderText({
      if (input$input_type=="3file_FC") {
        #only evaluate once all files loaded, standard valid is FALSE
        req(nrow(meta_tb())*ncol(meta_tb())*nrow(abund_tb())*ncol(abund_tb())*
              nrow(frac_tb())*ncol(frac_tb())>0)
        #check if at least one column name appears in all input datasets, required
        #for sample column. Do not continue
        common_column<- 0 < length(get_common_elements(colnames(meta_tb()),
                                                       colnames(abund_tb()),
                                                       colnames(frac_tb())))
        if (!common_column) {
          validate(
            paste0("The sample column needs to have the same name in each file. ",
                   "Yet there was not a single common column name in the files. ",
                   "Please check the names and make sure they have the same ",
                   "capitalisation.")
          )
        }
        
        #return empty text if all checks ok, conditionalpanel requires this empty
        #text!
        return("")
      }
      
      if (input$input_type=="3file_iso") {
        #only evaluate once all files loaded, standard valid is FALSE
        req(nrow(meta_tb())*ncol(meta_tb())>0,
            nrow(abund_tb())*ncol(abund_tb())>1,
            nrow(frac_tb())*ncol(frac_tb())>1)
        
        #check if at least one column name appears in all input datasets, required
        #for sample column. Do not continue
        common_column<- 0 < length(get_common_elements(colnames(meta_tb()),
                                                       colnames(abund_tb()),
                                                       colnames(iso_col_tb())))
        if (!common_column) {
          validate(
            paste0("The sample column needs to have the same name in each file. ",
                   "Yet there was not a single common column name in the files. ",
                   "Please check the names and make sure they have the same ",
                   "capitalisation.")
          )
        }
        #return empty text if all checks ok, conditionalpanel requires this empty
        #text!
        return("")
      }
      
      if (input$input_type=="2file") {
        #only evaluate once all files loaded, standard valid is FALSE
        req(nrow(meta_tb())*ncol(meta_tb())>0,input$iso_et_file)

        #check corrected isotopologue file validity
        if (!check_iso_input(iso_et_tb())=="OK") {
          validate(
            check_iso_input(iso_et_tb())
          )
        }
        
        #return empty text if all checks ok, conditionalpanel requires this empty
        #text!
        return("")
      }

    })
    outputOptions(output, 'input_valid', suspendWhenHidden=FALSE)
    
    
    ###metadata Column assignment
    #compound choices are all columns in abundance except sample
    #need return, otherwise nothing will be output except if the last if 
    #statement is true
    choises_sample<-reactive({
      if (input$input_type=="3file_FC") {
        return(get_common_elements(colnames(meta_tb()),
                            colnames(abund_tb()),
                            colnames(frac_tb())))
      }
      
      if (input$input_type=="3file_iso") {
        return(get_common_elements(colnames(meta_tb()),
                                   colnames(abund_tb()),
                                   colnames(iso_col_tb())))
      }
      
      if (input$input_type=="2file") {
        return(colnames(meta_tb()))
      }
      
    })
    
    #update sample choices when they change, select first option by default
    #close compounds box (makes it closed on first load
    #and all compounds will be reselected)
    observeEvent(choises_sample(), {
      shinyjs::hide(id = "compounds")
      updateSelectInput(inputId = "sample_column", choices = choises_sample(),
                        selected=choises_sample()[1]) 
    })
    
    #update derived abundance and fractional contribution column names
    observeEvent(input$sample_column,{
      if (input$input_type=="2file") {
        #only apply if data already uploaded
        req(nrow(frac_tb())*ncol(frac_tb())>1)
        
        #change sample column names (first column by design)
        abund_temp<-abund_tb()
        frac_temp<-frac_tb()
        colnames(abund_temp)[1]<-colnames(frac_temp)[1]<-input$sample_column
        abund_tb(abund_temp)
        frac_tb(frac_temp)
      }
    })
    
    #Normalisation choices are all metadata columns not selected for sample
    choises_norm<-reactive({
      colnames(meta_tb())[-which(colnames(meta_tb())==input$sample_column)]
    })
    
    #update normalisation column choices
    #select the first column name containing "norm" by default, 
    #if none present select "None" by default
    observeEvent(choises_norm(), {
      upd_selInp_guess_default(inputId = "norm_column",choices = choises_norm(),
                               prefString = "norm",noPick = "None",prefOrNo = T)
    })
    
    
    #Cohort column choices are all metadata columns not selected for sample
    #or normalisation
    choises_cohort<-reactive({
       colnames(meta_tb())[-which(colnames(meta_tb()) %in% c(input$sample_column,
                                                            input$norm_column))]
    })
    
    #update cohort column choices, making sure none is selectable
    #select the first column name containing "cohort" by default, 
    #if none present select first available column name by default
    #if no columns available at all will choose "None"
    observeEvent(choises_cohort(), {
      upd_selInp_guess_default(inputId = "cohort_column",
                               choices = choises_cohort(),prefString = "cohort",
                               noPick = "None",prefOrNo = F)
    })
    
    
    #Compound choices (all abundance columns minus sample column)
    choises_comp<-reactive({
      req(nrow(abund_tb())*ncol(abund_tb())>0)
      colnames(abund_tb())[-which(colnames(abund_tb())==
                                    input$sample_column)]
    })
    
    #update compound
    #select all options by default
    observeEvent(choises_comp(), {
      updateSelectInput(inputId = "compounds", choices = choises_comp(),
                        selected=choises_comp()) 
    })
    
    ## allow show/hide button to show/hide compounds dropdown list
    observeEvent(input$comp_showhide, {
      if(input$comp_showhide %% 2 == 1){
        shinyjs::show(id = "compounds")
      }else{
        shinyjs::hide(id = "compounds")
      }
    })
  
    #TODO/ the below didn't work to avoid an error coming from applying previous setting to new data before new settings are offered
    #try something
    #prepare reactive value that is FALSE right after input was changed
    #and only turns to TRUE after a short delay
    #when requiring a function to wait till this is true, this makes
    #sure downstream input options derived from earlier input
    #that then get applied to that input can also change before being applied
    #again
    # values <- reactiveValues(starting = TRUE)
    # session$onFlushed(function() {
    #   values$starting <- FALSE
    # })
    
    ###add a check to crosscheck samples and compounds between files and give 
    #appropriate errors or warnings
    output$sample_valid <- renderText({
      #requires sample column to be set first, because this will execute faster
      #than the finding possible sample column names and give an error otherwise
      #also requires a fractional contribution table to be calculated for this 
      #input type
      req(input$sample_column,nrow(frac_tb())*ncol(frac_tb())>1)

      # req(values$starting)
      
      #disables the merge button by default
      disable("merge")
      
      #generate error or warning messages if any
      check_output<-check_samples_compounds(
        meta_tb = meta_tb(),abund_tb = abund_tb(),frac_tb = frac_tb(),
        sample_column = input$sample_column,norm_column = input$norm_column)
      
      if (check_output$error) {
        validate(check_output$message)
      } else {
        outputtext<-check_output$message
      }
      
      #enables merge button if no validate errors from the function above
      enable("merge")
      
      return(outputtext)
    })
    outputOptions(output, 'sample_valid', suspendWhenHidden=FALSE)
    
    ## allow show/hide button to show/hide compounds dropdown list
    observeEvent(input$comp_showhide, {
      if(input$comp_showhide %% 2 == 1){
        shinyjs::show(id = "compounds")
      }else{
        shinyjs::hide(id = "compounds")
      }
      
    })

    ###construct merged table with button press
    #button press was required to avoid recalculation when updating eg. norm
    #column before cohort_column got updated, resulting in errors and crashes.
    #could not save as reactive as reactives always exists, and needed to test
    #when data is merged to only show table then. Can test reactivevalue, see 
    #tableviewer section
    observeEvent(input$merge,{
      meta_formatted_tb<-format_metadata(meta_tb = meta_tb(),
                                         sample_column = input$sample_column,
                                         factor_column = input$cohort_column,
                                         norm_column = input$norm_column)

      if (input$input_type=="3file_FC") {
        v$merged_tb<-merge_input(meta_tb = meta_formatted_tb,
                                 abund_tb = abund_tb(),
                                 frac_tb = frac_tb(),
                                 sample_col = input$sample_column,
                                 compounds = input$compounds)
      } else {
        print(iso_tb)
        v$merged_tb<-merge_input(meta_tb = meta_formatted_tb,
                                 abund_tb = abund_tb(),
                                 frac_tb = frac_tb(),
                                 iso_tb=iso_tb(),
                                 sample_col = input$sample_column,
                                 compounds = input$compounds)
      }
      
      enable("finish")
    })
    
    #reset merged table and disable continue button when input was changed
    #(until merged table calculated again)
    observe({
      input$input_type
      meta_tb()
      abund_tb()
      frac_tb()
      iso_tb()
      input$sample_column
      input$cohort_column
      input$norm_column
      input$compounds
      
      v$merged_tb<-tibble(NA)
      disable("finish")
      updateTabsetPanel(inputId = "downoptions", selected = "None")
    })
      
    ###give signal finished when continue button pressed
    observeEvent(input$finish,{
      v$finished<-T
    })
    
    #set right download option to show in ui
    observeEvent(v$merged_tb,{
      updateTabsetPanel(inputId = "downoptions", selected = "None")
      if (nrow(v$merged_tb)*ncol(v$merged_tb)>1) {
        if (local_version) {
          updateTabsetPanel(inputId = "downoptions", selected = "local")
          tibble_localoutput_server("localout",v$merged_tb)
        } else {
          updateTabsetPanel(inputId = "downoptions", selected = "web")
          tibble_weboutput_server("webout",v$merged_tb)
        }
      }
    })
    
    ###tableviewer
    #check whether any data uploaded to show table viewer
    output$any_upload <- reactive ({
      return(any(nrow(meta_tb())*ncol(meta_tb())>1,
                 nrow(abund_tb())*ncol(abund_tb())>1,
                 nrow(frac_tb())*ncol(frac_tb())>1))
    })
    outputOptions(output, 'any_upload', suspendWhenHidden=FALSE)
    
    #save available viewer choices depending on uploaded data
    #needed different approach for newly generated data, as a reactive always
    #exists used reactiveValues instead, NA before and dataframe after calc
    tablechoices<-reactive({
      choices <- NULL
      if (nrow(meta_tb())*ncol(meta_tb())>1) choices<-c(choices,"Metadata")
      if (nrow(abund_tb())*ncol(abund_tb())>1) choices<-c(choices,
                                                          "Abundance data")
      if (nrow(frac_tb())*ncol(frac_tb())>1) choices<-c(choices,
                                                "Fractional contribution data")
      if (nrow(iso_tb())*ncol(iso_tb())>1) choices<-c(choices,
                                                        "Isotopologue data")
      if (nrow(v$merged_tb)*ncol(v$merged_tb)>1) choices <- 
          c(choices,"Merged data")
      
      return(choices)
    })
    
    #update viewer choices if more is available, 
    #set to last processed in reactive. Don't support always showing last added
    #too much hassle
    observeEvent(tablechoices(), {
      updateSelectInput(inputId = "view_table", choices = tablechoices(),
                        selected=tablechoices()[length(tablechoices())]) 
    })
    
    #create table for output, depending on which table is to be shown
    output$table <- DT::renderDataTable({
      if (input$view_table == "Metadata") return(meta_tb())
      if (input$view_table == "Abundance data") return(abund_tb())
      if (input$view_table == "Fractional contribution data") return(frac_tb())
      if (input$view_table == "Isotopologue data") return(iso_tb())
      
      #need additional option to not load when variable is NA to avoid warning 
      #when resetting v$merged_tb due to changed input (table loads faster than
      #viewer option disappears)
      if (input$view_table == "Merged data" &
          !is.na(v$merged_tb[1,1])) return(v$merged_tb)
      
    },
    options = list(
      scrollX = TRUE )
    )
    
    # 
    # # close the R session when Chrome closes
    # session$onSessionEnded(function() { 
    #   stopApp()
    #   q("no") 
    # })
    
    return(v)
  })
}

#Wrap modules in caller function for testing-------------------------------------------------------------------
travis_inputcleanerApp<- function() {
  ui <- fluidPage(
    travis_cleaner_ui("InputClean")
  )
  
  server <- function(input, output, session) {
    travis_cleaner_server("InputClean",local_version = local_version)
  }
  
  shinyApp(ui, server)  
}

travis_inputcleanerApp()
