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
source(here::here("Functions and modules/TraVis_Pies_functions_V1.0.R"))


# Outputsettings Module-------------------------------------------------------------------
travis_outputsettings_ui <- function(id) {
  fluidPage(
    #App explanation text
    titlePanel(h1("Travis Pies v1.0 Output")),
    htmlOutput(NS(id,"text")),

    #Saving options
    fluidRow(
      column(4,
             radioButtons(NS(id,"plottype"),
                         label=paste0("Which plot types should be saved?"),
                         choices= c("Stand-alone","Pathway-compatible","Both"))
      ),
      column(4,
             radioButtons(NS(id,"format"),
                          label=paste0("Which format should be used?"),
                          choices= c("png","tiff"))
      )
    ),

    fluidRow(
      column(12,
             selectInput(NS(id,"compounds"),
                         label=paste0("Select compounds to be plotted, by ",
                         "default all in input"),
                         choices= NULL, multiple = T,width="100%")
      )
    )
  )
}

travis_outputsettings_server <- function(id,tb) {
  moduleServer(id,function(input, output, session) {

    #Write text to put on top as explanation, need to use server output to be able
    #to write multiple lines
    output$text <- renderText({
      paste0("<b>Apply settings from visualisation screen to all selected ",
             "compounds output a plot for each.</b><br/><br/>")

    })

    ###Obtain compound choices
    choises_compound<-reactive({
      compnames<-colnames(tb())[(which(colnames(tb())=="datatype")+1):
                                  ncol(tb())]
      return(compnames)
    })

    observeEvent(choises_compound(), {
      updateSelectInput(inputId = "compounds", choices = choises_compound(),
                        selected=choises_compound())
    })

    return(
      list(
        plottype = reactive({ input$plottype }),
        format = reactive({ input$format }),
        compounds = reactive({ input$compounds })

      )
    )
  })
}

#Local output Module-------------------------------------------------------------------
travis_outputlocal_ui <- function(id) {
  fluidPage(
    #load shiny package elements
    shinyFeedback::useShinyFeedback(),
    useShinyjs(),    ## IMPORTANT: so shiny knows to use the shinyjs library
    
    #load settings ui, need to use NS(ID,...) in nested modules!
    travis_outputsettings_ui(NS(id,"outputsettings")), 
    
    #directory selection
    strong(paste0("Select a directory to save to or copypaste the ",
                  "desired folderpath in textbox below")),
    HTML(paste0("<p style='color:orange'>Warning: Will overwrite existing ",
                "files with the same name and extension!</p>")),
    fluidRow(
      column(3,
             shinyDirButton(NS(id,"directory_button"),
                            "Pick directory",
                            "Save plots to ...",
                            style = "margin-top: 20px;margin-bottom: 25px;")
      ),
      
      column(4,
             textInput(NS(id,"directory_text"),"","",width="100%")
      )
    ),
    
    #Button to start saving process
    fluidRow(
      column(2,
             disabled(actionButton(NS(id, "save_plots"),
                                   label = "Save plots to directory"))
      )
    )
  )
}

travis_outputlocal_server <- function(id,v_settings,tb) {
  moduleServer(id,function(input, output, session) {
    
    out_settings<-travis_outputsettings_server("outputsettings",
                                               tb = tb)
    
    ###Directory selection
    
    #Allow choosing any directory with the button
    # Set directory text to button directory when picked
    observeEvent(input$directory_button,{
      disable("save_plots")
      volumes<-getVolumes() #get all volumes of pc
      shinyDirChoose(input, "directory_button", roots=volumes(),session=session)
      savepath<-parseDirPath(volumes, input$directory_button)
      updateTextInput(inputId = "directory_text",
                      value = savepath)
      req(dir.exists(savepath))
      enable("save_plots")
    })
    
    ## Save directorytext, convert to R string if required
    #and check if the directory exists
    observeEvent(input$directory_text,{
      disable("save_plots")
      req(input$directory_text)
      target_savepath<-trimws(gsub("\\\\", "/",input$directory_text))
      
      #show feedback and disable plot button if directory does not exist
      shinyFeedback::feedbackDanger("directory_text",
                                    !dir.exists(target_savepath),
                                    "This directory doesn't exist.")
      
      req(dir.exists(target_savepath))
      enable("save_plots")
    })
    
    ## Generate pies to destination folder
    #load in all input!
    observeEvent(input$save_plots,{
      #disable save button while still calculating
      disable("save_plots")

      #get path to save to, corrected if required
      target_savepath<-trimws(gsub("\\\\", "/",input$directory_text))

      
      #set which charts to generate depending on requested plottype
      detail_charts<-F
      pathway_charts<-F

            if (out_settings$plottype() %in% c("Stand-alone","Both")) {
        detail_charts<-T
      }
      if (out_settings$plottype() %in% c("Pathway-compatible","Both")) {
        pathway_charts<-T
      }

      withProgress(message= "Saving plots", value=0, {
        
        #loop over each compound in input tibble
        for (compound in out_settings$compounds()) {
          printmessage<-paste0("Compound ",
                               which(out_settings$compounds()==compound),
                               " of ",length(out_settings$compounds()))
          print(printmessage)
          
          # Increment the progress bar, and update the detail text.
          incProgress(1/length(out_settings$compounds()), 
                      detail = printmessage)
          
          
          generate_pie(tb(),compound=compound,detail_charts=detail_charts,
                       pathway_charts=pathway_charts,savepath=target_savepath,
                       normalize=v_settings$norm,fact_name=v_settings$fact_name,
                       fact_order=v_settings$fact_order,
                       label_decimals=v_settings$label_decimals,
                       percent_add = v_settings$percent_add ,
                       FC_position=v_settings$FC_position,
                       min_lab_dist =v_settings$min_lab_dist,
                       circlelinecolor=v_settings$circlelinecolor,
                       circlelinetypes=v_settings$circlelinetypes,
                       maxcol_facet=v_settings$maxcol_facet,
                       include_name=v_settings$include_name,
                       show_P=v_settings$show_P,
                       col_labeling=v_settings$col_labeling,
                       alpha=v_settings$alpha,
                       otherfontsize=v_settings$otherfontsize,
                       font=v_settings$font,
                       legendtitlesize=v_settings$legendtitlesize,
                       cohortsize=v_settings$cohortsize,
                       include_legend=v_settings$include_legend,
                       format=out_settings$format())
          
          
        }  
      })
    
    print("Finished")
    enable("save_plots")
      
    })
  })
}


#Web output Module-------------------------------------------------------------------
travis_outputweb_ui <- function(id) {
  fluidPage(
    #load shiny package elements
    shinyFeedback::useShinyFeedback(),
    useShinyjs(),    ## IMPORTANT: so shiny knows to use the shinyjs library
    
    #load settings ui, need to use NS(ID,...) in nested modules!
    travis_outputsettings_ui(NS(id,"outputsettings")), 
    
    #Button to start making plots
    fluidRow(
      column(2,
             downloadButton(NS(id, "create_download_plots"),
                            label = "Download plots")
      )
    )
  )
}

travis_outputweb_server <- function(id,v_settings,tb) {
  moduleServer(id,function(input, output, session) {
    
    out_settings<-travis_outputsettings_server("outputsettings",
                                               tb = tb)
    enable("create_download_plots")
    
    # Create and download the plots
    output$create_download_plots = downloadHandler(
      filename = 'pies.zip',
      content = function(file){
        filelist<-NULL
        
        # Set temporary working directory
        owd <- setwd( tempdir())
        on.exit( {
          setwd( owd)
        })
        
        #see where temp files are saved in print console
        print(paste0("Saving to ",getwd()))
        

        #disable create button while still calculating
        disable("create_download_plots")
        
        #set which charts to generate depending on requested plottype
        detail_charts<-F
        pathway_charts<-F

        if (out_settings$plottype() %in% c("Stand-alone","Both")) {
          detail_charts<-T
          filelist<-c(filelist,"Pie charts/")
        }
        if (out_settings$plottype() %in% c("Pathway-compatible","Both")) {
          pathway_charts<-T
          filelist<-c(filelist,"Pie charts pathway/")
        }
        
        withProgress(message= "Making plots", value=0, {
          
          #loop over each compound in input tibble
          for (compound in out_settings$compounds()) {
            printmessage<-paste0("Compound ",
                                 which(out_settings$compounds()==compound),
                                 " of ",length(out_settings$compounds()))
            print(printmessage)
            
            # Increment the progress bar, and update the detail text.
            incProgress(1/length(out_settings$compounds()), 
                        detail = printmessage)
            
            
            generate_pie(tb(),compound=compound,detail_charts=detail_charts,
                         pathway_charts=pathway_charts,savepath=getwd(),
                         normalize=v_settings$norm,fact_name=v_settings$fact_name,
                         fact_order=v_settings$fact_order,
                         label_decimals=v_settings$label_decimals,
                         percent_add = v_settings$percent_add ,
                         FC_position=v_settings$FC_position,
                         min_lab_dist =v_settings$min_lab_dist,
                         circlelinecolor=v_settings$circlelinecolor,
                         circlelinetypes=v_settings$circlelinetypes,
                         maxcol_facet=v_settings$maxcol_facet,
                         include_name=v_settings$include_name,
                         show_P=v_settings$show_P,
                         col_labeling=v_settings$col_labeling,
                         alpha=v_settings$alpha,
                         otherfontsize=v_settings$otherfontsize,
                         font=v_settings$font,
                         legendtitlesize=v_settings$legendtitlesize,
                         cohortsize=v_settings$cohortsize,
                         include_legend=v_settings$include_legend,
                         format=out_settings$format())
          }  
        })
        
        print("Finished")
        enable("create_download_plots")
        
        # Zip them up
        zip( file, filelist)
      }
    )
    
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


travis_outputlocal_app<- function() {
  ui <- fluidPage(
    travis_outputlocal_ui("output")
  )
  

  server <- function(input, output, session) {
    travis_outputlocal_server("output",v_settings = 
                           do.call("reactiveValues",example_settings),
                        tb = reactive(example_tb))
  }
  
  shinyApp(ui, server)  
}

travis_outputweb_app<- function() {
  ui <- fluidPage(
    travis_outputweb_ui("output")
  )
  
  
  server <- function(input, output, session) {
    travis_outputweb_server("output",v_settings = 
                                do.call("reactiveValues",example_settings),
                              tb = reactive(example_tb))
  }
  
  shinyApp(ui, server)  
}

travis_outputweb_app()


