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

#function to check which elements are common among all vectors
get_common_elements<-function(...){
  #list all input objects, check if all elements are vectors
  vectorlist<-list(...)
  for(i in vectorlist) {
    if (!is.atomic(i)) stop(paste0(i," is not an atomic vector. "))
    if (is.matrix((i))) stop(paste0(i," is a matrix, not an atomic vector."))
  }
  
  #obtain elements common to all vectors in list
  Reduce(intersect, vectorlist)
}

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

#helper function for server to check samples across input files and check
#compounds present to generate informative errors and warnings where
#appropriate
check_samples_compounds<-function(meta_tb,abund_tb,frac_tb,sample_column,
                                        norm_column){
  sample_warnings<-NULL
  
  #Check if all samples in meta table are present abund table
  samples_miss_abund<-!all(pull(meta_tb,sample_column) %in%
                             pull(abund_tb,sample_column))
  if (samples_miss_abund) {
    validate(paste0("Sample from metadata file missing in abund file. If ",
                    "the correct sample column is chosen, verify samples ",
                    "and sample names in both files."))
  }
  
  #Check if all samples in meta table are present frac table
  samples_miss_frac<-!all(pull(meta_tb,sample_column) %in%
                            pull(frac_tb,sample_column))
  if (samples_miss_frac) {
    validate(paste0("Sample from metadata file missing in frac file. If ",
                    "the correct sample column is chosen, verify samples ",
                    "and sample names in both files."))
  }
  
  #check if normalisation column is numeric if there is a column specified
  if (norm_column != "None") {
    if (!is.numeric(pull(meta_tb,norm_column))) {
      validate(paste0("The chosen normalisation column does not contain ",
                      "numbers. Please pick the right column or check the input if this is it."))
    }
  }
  
  
  #Warn if more samples present in abund or frac file than in meta
  samples_ignored<-!all( abund_tb[,sample_column] %in%
                           meta_tb[,sample_column],
                         frac_tb[,sample_column] %in%
                           meta_tb[,sample_column])
  if (samples_ignored) {
    sample_warnings<-paste0(
      c(sample_warnings,
        paste0("Samples from abund and or frac file missing in ",
               "metadatafile. These samples will be removed from the ",
               "analysis."))
    )
  }
  #Warn if compounds present in abund file not frac file 
  #will be 100% unlabeled
  comp_ab_only<-colnames(abund_tb)[which(!colnames(abund_tb)%in%
                                           colnames(frac_tb))]
  if (length(comp_ab_only)>0) {
    sample_warnings<-
      c(sample_warnings,
        paste0("Following compounds only in abundance file, will be ",
               "considered fully unlabeled: ",
               paste(comp_ab_only,collapse = ", ")))
  }
  
  #Warn if compounds present in frac file not abund file
  #will be removed
  comp_fc_only<-colnames(frac_tb)[which(!colnames(frac_tb)%in%
                                          colnames(abund_tb))]
  if (length(comp_fc_only)>0) {
    sample_warnings<-
      c(sample_warnings,
        paste0("Following compounds only in fractional contribution file, ",
               " will be removed: ",
               paste(comp_fc_only,collapse = ", ")))
  }
  
  #Warn if compounds have 0 abundance in every sample, they will be dropped
  compounds_notdetected<-colnames(
    select(abund_tb,-where(has_nonzero))
  )
  
  if (length(compounds_notdetected)>0) {
    sample_warnings<-
      c(sample_warnings,
        paste0("Following compounds are never detected (abundance always 0), ",
               " and will be removed: ",
               paste(compounds_notdetected,collapse = ", ")))
  }
  
  #return empty text if no warnings, else give them in orange text (html)
  if (length(sample_warnings)>0) {
    return(paste("<b><p style='color:orange'>Warning: </b>",sample_warnings,
                 "</p>", sep = "<br/>"))
  } else {
    return("")
  }
}

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
    titlePanel(h1("Travis Pies v1.0 Input creator")),
    htmlOutput(NS(id,"text")),
    
    #input files
    h3(paste0("Upload all 3 input files below to continue, uploaded data can be ",
              "viewed below")),
    fluidRow(
      tags$style(HTML("div.input-group {margin-bottom: -30px;}")),
      column(
        3,fileInput(NS(id,"meta_file"), "Metadata",accept = ".csv")),
      column(
        3,fileInput(NS(id,"abund_file"), "Abundance data",accept = ".csv")),
      column(
        3,fileInput(NS(id,"frac_file"), "Fractional contribution data",accept = ".csv"))
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
      
      #how to optionally store merged file for later us
      # conditionalPanel(
      #   "output.merge_done == true",
      #   ns=NS(id),
      #   tabsetPanel(
      #     id = "downoptions",
      #     type = "hidden",
      #     selected =  "web",
      #     tabPanel("local", 
      #              tibble_localoutput_ui(NS(id,"localout"))
      #     ),
      #     tabPanel("web", 
      #              tibble_weboutput_ui(NS(id,"webout"))
      #     )
      #   )
      # )
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
    
    #Write text to put on top as explanation, need to use server output to be able 
    #to write multiple lines
    output$text <- renderText({
      paste0("<b>Upload 3 input files, metadata, abundance data and fractional ",
             "contribution data.</b><br/>See example input files provided with",
             " this application. All contain an identically named column",
             " with sample names. Metadata optionally contains cohorts and ",
             "normalisation columns. The other two contain only the sample ",
             "name column and their respective data columns named after the ",
             "compounds which they represent.<br/>Fractional contribution ",
             "data can be provided as a fraction between slightly less than 0 ",
             "and 1 or a percentage between slightly less than 0% and 100%.",
             "<br/> Normalisation happens by dividing abundances by the ",
             "normalisation value provided.")
    })
    
    #Fileloading input data
    meta_tb<-reactive({
      read_csv_clean(input$meta_file$datapath,remove_empty = T)
    })
    
    abund_tb<-reactive({
      read_csv_clean(input$abund_file$datapath,remove_empty = T)
    })
    frac_tb<-reactive({
      read_csv_clean(input$frac_file$datapath,remove_empty = T)
    })
    
    #add a check to see uploaded data is valid
    output$input_valid <- renderText({
      #only evaluate once all files loaded, standard valid is FALSE
      req(input$meta_file,input$abund_file,input$frac_file)
      
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
      ""
    })
    outputOptions(output, 'input_valid', suspendWhenHidden=FALSE)
    
    
    
    ###metadata Column assignment
    #compound choices are all columns in abundance except sample
    choises_sample<-reactive({
      get_common_elements(colnames(meta_tb()),
                          colnames(abund_tb()),
                          colnames(frac_tb()))
      
    })
    
    #update sample choices when they change, select first option by default
    #close compounds box (makes it closed on first load
    #and all compounds will be reselected)
    observeEvent(choises_sample(), {
      shinyjs::hide(id = "compounds")
      updateSelectInput(inputId = "sample_column", choices = choises_sample(),
                        selected=choises_sample()[1]) 
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
      colnames(abund_tb())[-which(colnames(abund_tb())==input$sample_column)]
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
    
    ###add a check to crosscheck samples and compounds between files and give 
    #appropriate errors or warnings
    output$sample_valid <- renderText({
      #requires sample column to be set first, because this will execute faster
      #than the finding possible sample column names and give an error otherwise
      req(input$sample_column)
      #disables the merge button by default
      disable("merge")
      
      outputtext<-check_samples_compounds(meta_tb = meta_tb(),abund_tb = abund_tb(),
                              frac_tb = frac_tb(),
                              sample_column = input$sample_column,
                              norm_column = input$norm_column)
      
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

      v$merged_tb<-merge_input(meta_tb = meta_formatted_tb,abund_tb = abund_tb(),
                               fraccon_tb = frac_tb(),sample_col = input$sample_column,
                               compounds = input$compounds)
      enable("finish")
    })
    
    #reset merged table and disable continue button when input was changed
    #(until merged table calculated again)
    observe({
      input$meta_file
      input$abund_file
      input$frac_file
      input$sample_column
      input$cohort_column
      input$norm_column
      input$compounds
      
      v$merged_tb<-data.frame(NA)
      disable("finish")
    })
      
    ###give signal finished when continue button pressed
    observeEvent(input$finish,{
      v$finished<-T
    })
    
    # #check to see if merged table is ready, to start download server
    # #and show download UI only then
    # output$merge_done <- reactive({
    #   #add ==true explicitly because if it returns a 1 element tibble with
    #   if (!is.na(v$merged_tb[1,1])) {
    #     return(TRUE)
    #   }
    # })
    # outputOptions(output, 'merge_done', suspendWhenHidden=FALSE)
    
    #set right download option to show in ui
    observeEvent(v$merged_tb,{
      updateTabsetPanel(inputId = "downoptions", selected = "none")
      if (!is.na(v$merged_tb[1,1])) {
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
      return(any(!is.null(input$meta_file),!is.null(input$abund_file),
                 !is.null(input$frac_file)))
    })
    outputOptions(output, 'any_upload', suspendWhenHidden=FALSE)
    
    #save available viewer choices depending on uploaded data
    #needed different approach for newly generated data, as a reactive always
    #exists used reactiveValues instead, NA before and dataframe after calc
    tablechoices<-reactive({
      choices <- NULL
      if (!is.null(input$meta_file)) choices<-c(choices,"Metadata")
      if (!is.null(input$abund_file)) choices<-c(choices,"Abundance data")
      if (!is.null(input$frac_file)) choices<-c(choices,
                                                "Fractional contribution data")
      if (!is.na(v$merged_tb[1,1])) choices <- c(choices,"Merged data")
      
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
    travis_cleaner_server("InputClean",local_version = F)
  }
  
  shinyApp(ui, server)  
}

travis_inputcleanerApp()
