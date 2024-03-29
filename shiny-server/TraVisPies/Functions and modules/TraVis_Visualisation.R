# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: This module shows the effect of settings on the pie chart plot for
# a compound from a provided data set. The user can interactively select the 
# desired settings, which will be saved in order to transmit them to the
# output module where they will be applied to generate similar plots for all 
# desired compounds from thedataset

###Input
# A curated TraVis input tibble as generated by the InputCleaner module

###Output
# ReactiveValues object with selected settings for generating pie chart plots.
# This object also includes an element called ready that is FALSE until 
# valid settings were chosen

# Functions and libraries ---------------------------------------------------------------
library(here)    #to make r source from file location instead of magic stuff
library(shiny)
library(shinyjs)      #for show/hide button
library(colourpicker) #to allow colours to be chosen

#library for calling helper functions from r script
library(here)    #to make r source from file location instead of magic stuff
source(here::here("Functions and modules/TraVis_Pies_functions.R"))

#load font library. For windows only it loads these fonts for bitmap output
#as well, not required for other operating systems. Don't  import library 
#unless never done before on this system or unless 
if (.Platform$OS.type=="windows") {
  library(extrafont)
  loadfonts(device="win",quiet = T)
}

#check validity of numeric inputs
valid_num_input<-function(input,min=NULL,max=NULL,label=input){
  #check for both null and NA to be sure (normally NA it seems)
  if (is.null(input)|is.na(input)) {
    validate(
      paste0(label,": A value needs to be provided")
    )
  } else {
    if(!is.null(min)){
      if (input<min) {
        validate(
          paste0(label,": The value needs to be at least ",min)
        )
      }
    }
    
    if(!is.null(max)){
      if (input>max) {
        validate(
          paste0(label,": The value needs to be at most ",max)
        )
      }
    }
  }
}

#Visualisation module UI and server-------------------------------------------------------------------
travis_viz_ui <- function(id) {
  #predefine named linetype list for selecting
  linetype_names<-list("None","Solid","dashed","dotted","dot-dashed",
                      "doubledashed")
  linetype_values<-list(0,1,2,3,4,6)
  names(linetype_values)<-linetype_names
  
  #actual UI
  fluidPage(
    #load shiny package elements
    useShinyjs(),    ## IMPORTANT: so shiny knows to use the shinyjs library
    
    #set color for error messages, and adapt standard spacing
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
      color: red;
      }
      div.checkbox {margin-top: 0px;margin-bottom: -20px;}
      div.selectize-input {margin-bottom: -15px;}
      "))
    ),

    #Always shown crucial inputs for plot
    h3(paste0("Modify compound, cohort variable,cohort order and ",
    "normalisation")),
    fluidRow(
      column(5,
             selectInput(NS(id,"compound"),label = "Compound",choices=NULL)
      ),
      column(2,
             selectInput(NS(id,"factor"),label = "Cohort variable",choices=NULL)
      ),
      column(4,
             selectInput(NS(id,"fac_levels"),label = "Cohorts",choices=NULL,
                         multiple = T)
      )
    ),
    fluidRow(
      column(12,
             disabled(
               checkboxInput(NS(id,"norm"),
                             label = paste0("Use normalized abundances, ",
                             "disabled if no normalisation column was ",
                             "provided or present"),
                             value=T,
                             width = "100%")
               )
      ),
      column(12,
             disabled(
               checkboxInput(NS(id,"P_isotopologues"),
                             label = paste0("Show * in cohort name if ",
                             "significant isotopologue difference, ",
                             "disabled if no isotopologue data is ",
                             "present"),
                             value=F,
                             width = "100%")
             )
      )
    ),
    
    #Chart layout option, not shown by default. inline-block for displaying
    # show/hide button next to text
    div(style="display:inline-block",
        h3(paste0("Chart layout options"))),
    div(style="display:inline-block; ",
        actionButton(NS(id,"layout_showhide"), label = "show / hide")),
    
    shinyjs::hidden(
      htmltools::div(
        id=NS(id,"Layout_options"),
        fluidRow(
          column(3,
                 numericInput(NS(id,"maxcol_facet"),"Maximum charts per row",4,
                              min=1,step = 1)
          ),
          column(2,
                 radioButtons(NS(id,"FC_position"),"Position FC label",
                              c("center","slice"))
          ),
          column(4,
                 disabled(
                   numericInput(NS(id,"min_lab_dist"),
                                "Slice label distance from center",
                                0.42,min=0,max=1,step=0.01)
                  )

          )
        ),
        fluidRow(
          column(3,
                 numericInput(NS(id,"label_decimals"),
                              "FC label decimals",
                              0,min=0,step=1)
          ),
          column(3,
                 checkboxInput(NS(id,"percent_add"),
                               label = "Add '%' to FC label"),
                 checkboxInput(NS(id,"include_legend"),label = "Add legend",
                               value = T),
                 checkboxInput(NS(id,"log_abund"),
                               label = "Use base 10 logarithmic scale for abundance",
                               value = F)
                 
          ),
          column(3,
                 checkboxInput(NS(id,"show_P"),
                               label = "Show P values",
                               value = T),
                 checkboxInput(NS(id,"include_name"),
                               label = "Add compound name as title",
                               value = T)
          )
        ),
        
        fluidRow(
          column(4,
                 colourInput(NS(id,"colour_lab"),
                             "Pick colour labeled fraction",
                             value = "#ffd966")
          ),
          column(4,
                 colourInput(NS(id,"colour_unlab"),
                             "Pick colour unlabeled fraction",
                             value = "#bfbfbf")
          ),
          column(4,
                 colourInput(NS(id,"colour_circle"),
                             "Pick circle line colour",
                             value = "gray")
          )
        ),
        fluidRow(
          column(4,
                 sliderInput(NS(id,"alpha"),
                             "Set opacity (arrow buttons to finetune)",
                             value = 0.7, min=0, max = 1, step = 0.01)
          ),
          column(8,
                 div(style="margin-bottom: 15px;",
                     (paste0("Tip: copy and save the colour string to ",
                             "easily replicate your favorite colours")))
          )
        ),
        fluidRow(
          column(3,
                 selectInput(NS(id,"type_circle_100"),
                              "Largest circle linetype",
                              choices = linetype_values,selected = 1)
          ),
          column(3,
                 selectInput(NS(id,"type_circle_75"),
                              "2nd largest circle linetype",
                              choices = linetype_values,selected = 1)
          ),
          column(3,
                 selectInput(NS(id,"type_circle_50"),
                              "2nd smallest circle linetype",
                              choices = linetype_values,selected = 1)
          ),
          column(3,
                 selectInput(NS(id,"type_circle_25"),
                              "Smalles circle linetype",
                              choices = linetype_values,selected = 1)
          )
        ),
      )
    ),
    
    #Chart font options, not shown by default. inline-block for displaying
    # show/hide button next to text
    div(style="display:inline-block",
        h3(paste0("Chart font options"))),
    div(style="display:inline-block; ",
        actionButton(NS(id,"font_showhide"), label = "show / hide")),
    shinyjs::hidden(
      htmltools::div(
        id=NS(id,"font_options"),
        fluidRow(
          column(3,
                 selectInput(NS(id,"font"),
                             label = "Select font if available",
                             choices=unique(c(fonts(),"sans")),
                             selected = "sans")
          ),
          column(3,
                 numericInput(NS(id,"fontsize_cohort"),
                              "Fontsize cohort labels",
                              18,min=0,step=1)
          ),column(3,
                   numericInput(NS(id,"fontsize_legend"),
                                "Fontsize legend title",
                                16,min=0,step=1)
          ),column(3,
                   numericInput(NS(id,"fontsize_other"),
                                "Fontsize other text",
                                16,min=0,step=1)
          )
          
        )
      )
    ),
    #button to reset to default settings
    div(style="display:inline-block; ",
        actionButton(NS(id,"reset"), label = "Reset all options")),
 
    
    #output image, specify width and height to make sure other 
    #elements don't overlap
    imageOutput(NS(id,"pie_image"),width = "100%",height="100%"),
    
    #output caption
    div(style="display:inline-block; ",
        actionButton(NS(id,"caption_showhide"), label = "show / hide caption")),
    shinyjs::hidden(
      textOutput(NS(id,"caption"))
    )
  )
}

travis_viz_server <- function(id,tb) {
  moduleServer(id,function(input, output, session) {
    #prepare reactivevalues for passing settings out of function
    v_settings<-reactiveValues(ready=F,norm=NA,P_isotopologues = NA,
                      fact_name = NA,fact_order=NA,normalize = NA,
                      label_decimals = NA,
                      min_lab_dist = NA,
                      percent_add = NA,
                      FC_position = NA,
                      log_abund = NA,
                      col_labeling = NA,
                      circlelinetypes = NA,
                      circlelinecolor = NA,
                      maxcol_facet= NA,
                      include_name = NA,alpha=NA,
                      font=NA,otherfontsize = NA,
                      legendtitlesize =NA,
                      cohortsize = NA,
                      include_legend = NA,show_P=NA)
    
    #enables normalize box if table contains NormAbund datatype
    observe({
      req(!is.na(tb()[1,1]))
      if ("NormAbund" %in% tb()$datatype) {
        enable("norm")
      } else {
        disable("norm")
        updateCheckboxInput(inputId = "norm",value= F )
      }
    })
    
    observe({
      req(!is.na(tb()[1,1]))
      if ("Isotopologues" %in% tb()$datatype) {
        enable("P_isotopologues")
      } else {
        disable("P_isotopologues")
        updateCheckboxInput(inputId = "P_isotopologues",value= F )
      }
    })
    
    
    #enables min_lab_dist setting if FC_position set to slice
    observe({
      if (input$FC_position == "slice") {
        enable("min_lab_dist")
      } else {
        disable("min_lab_dist")
      }
    })
      
    
    ###Inputs for figure
    ###Give as input options possible compounds, factors and factor levels
    #compound choices
    choises_compound<-reactive({
      req(!is.na(tb()[1,1]))
      compnames<-colnames(tb())[(which(colnames(tb())=="datatype")+1):
                                  ncol(tb())]
      return(compnames)
    })
    observeEvent(choises_compound(), {
      updateSelectInput(inputId = "compound", choices = choises_compound()) 
    })
    
    #factor choices
    choises_factor<-reactive({
      req(!is.na(tb()[1,1]))
      colnames(tb())[2:(which(colnames(tb())=="datatype")-1)]
    })
    observeEvent(choises_factor(), {
      updateSelectInput(inputId = "factor", choices = choises_factor(),
                        selected = choises_factor()[1])
    })
    
    #factor levels: observeEvent to see when changes, still need req because the event can
    #activate before the factor input was updated with the new default choice
    observeEvent(input$factor,{
      req(!is.na(tb()[1,1]))
      req(input$factor)
      choises_fac_levels<-unique(pull(tb(),input$factor))
      updateSelectInput(inputId = "fac_levels", choices = choises_fac_levels,
                        selected=choises_fac_levels)
    })
    
    
    ###Show/hide Layout options
    observeEvent(input$layout_showhide, {
      if(input$layout_showhide %% 2 == 1){
        shinyjs::show(id = "Layout_options")
      }else{
        shinyjs::hide(id = "Layout_options")
      }
    })
    
    ###Show/hide font options
    observeEvent(input$font_showhide, {
      if(input$font_showhide %% 2 == 1){
        shinyjs::show(id = "font_options")
      }else{
        shinyjs::hide(id = "font_options")
      }
    })
    
    #make reactives vector colours and cirle linetypes for input to function
    #convert to numeric for functions compatibility
    col_labeling<-reactive(c(input$colour_unlab,input$colour_lab))
    circlelinetypes<-reactive(as.numeric(c(input$type_circle_25,
                                           input$type_circle_50,
                                           input$type_circle_75,
                                           input$type_circle_100)))
    
    ###Reset all to default
    observeEvent(input$reset, {
      #reset compound, factor and normalisation choices
      updateSelectInput(inputId = "compound", choices = choises_compound()) 
      updateSelectInput(inputId = "factor", choices = choises_factor())
      choises_fac_levels<-unique(pull(tb(),input$factor))
      updateSelectInput(inputId = "fac_levels", choices = choises_fac_levels,
                        selected=choises_fac_levels)
      if ("NormAbund" %in% tb()$datatype) {
        updateCheckboxInput(inputId = "norm",value= T )
      } else {
        updateCheckboxInput(inputId = "norm",value= F )
      }
      updateCheckboxInput(inputId = "P_isotopologues",value= F )

      #General layout
      updateNumericInput(inputId = "maxcol_facet",value = 4)
      updateRadioButtons(inputId = "FC_position",selected ="center")
      updateNumericInput(inputId = "min_lab_dist",value = 0.42)
      updateNumericInput(inputId = "label_decimals",value = 0)
      updateCheckboxInput(inputId = "percent_add",value= F)
      updateCheckboxInput(inputId = "log_abund",value= F)
      updateCheckboxInput(inputId = "include_name",value= T)
      updateCheckboxInput(inputId = "show_P",value= T)
      updateCheckboxInput(inputId = "include_legend",value= T)
      
      #Colours reset
      updateColourInput(inputId = "colour_lab", value = "#ffd966", 
                        session = getDefaultReactiveDomain() )
      updateColourInput(inputId = "colour_unlab", value = "#bfbfbf", 
                        session = getDefaultReactiveDomain() )
      updateColourInput(inputId = "colour_circle",value= "gray", 
                        session = getDefaultReactiveDomain() )
      updateSliderInput(inputId = "alpha",value = 0.7)
      
      #circle linetypes reset
      updateSelectInput(inputId = "type_circle_100",selected=1)
      updateSelectInput(inputId = "type_circle_75",selected=1)
      updateSelectInput(inputId = "type_circle_50",selected=1)
      updateSelectInput(inputId = "type_circle_25",selected=1)
      
      #fontoptions reset
      updateSelectInput(inputId = "font", selected = "sans") 
      updateNumericInput(inputId = "fontsize_cohort",value = 18)
      updateNumericInput(inputId = "fontsize_legend",value = 16)
      updateNumericInput(inputId = "fontsize_other",value = 16)
    })
    
    ###datatables: split calculation over different table so that if 
    #parameter that only affects later step is modified, not everything has
    #to be recalculated
    
    #Get table for selected compound
    comptb<- reactive({
      v_settings$ready<-F
      req(input$factor)
      req(input$fac_levels)
      req(all(input$fac_levels %in% unique(pull(tb(),input$factor))))
      
      comptb<-obtain_compounddata(tb(),compound = input$compound,
                          fact_name = input$factor,fact_order=input$fac_levels,
                          normalize = input$norm)
      
      #either remove isotopologues or parse them into one entry per isotopologue
      #then make sure the value column is numeric for further analysis
      if (!input$P_isotopologues) {
        comptb<-filter(comptb,!datatype=="Isotopologues") %>%
          mutate(across(!!input$compound,as.numeric))
      } else {
        comptb<-parse_isos_torow(comptb,valuecolumn = input$compound) %>%
          mutate(across(!!input$compound,as.numeric))
      }
      return(comptb)
    })

    #Get table for plot
    plottb<- reactive({
      v_settings$ready<-F
      #somehow comptb() is calculated before it has any data, so require
      #that the dataype column contains Abund which it always should
      req("Abund" %in% comptb()$datatype)
      
      #check validity of numeric inputs
      valid_num_input(input$label_decimals,min = 0,label="FC label decimals")
      valid_num_input(input$min_lab_dist,min = 0,max=1,
                      label="Slice label distance from center")
      
      tb<-prepare_slicedata(comptb(),fact_name = input$factor,
                            compound=input$compound,
                            label_decimals = input$label_decimals,
                            min_lab_dist = input$min_lab_dist,
                            percent_add = input$percent_add,
                            FC_position = input$FC_position,
                            P_isotopologues = input$P_isotopologues)
      return(tb)
    })
    
    #make pie chart plot
    pies<-reactive({
      v_settings$ready<-F
      req(plottb())
      
      #test what happens if fraction input maxcol facet
      #check validity of numeric inputs
      valid_num_input(input$maxcol_facet,min = 1,
                      label="Maximum charts per row")
      valid_num_input(input$fontsize_cohort,min = 0,
                      label="Fontsize cohort labels")
      valid_num_input(input$fontsize_legend,min = 0,
                      label="Fontsize legend title")
      valid_num_input(input$fontsize_other,min = 0,
                      label="Fontsize other text")
      
      
      pies<-make_piechart(plottb(),compound=input$compound,
                          fact_name = input$factor,
                          log_abund = input$log_abund,
                          circlelinecolor = input$colour_circle,
                          circlelinetypes = circlelinetypes(),
                          maxcol_facet= input$maxcol_facet,
                          include_name = input$include_name,
                          show_P=input$show_P,
                          col_labeling = col_labeling(), alpha=input$alpha,
                          font=input$font,otherfontsize = input$fontsize_other,
                          legendtitlesize =input$fontsize_legend,
                          cohortsize = input$fontsize_cohort,
                          include_legend = input$include_legend)
      
      return(pies)
    })
    
    #update input settings to those used for latest pies plot
    observeEvent(pies(),{
      #save settings and other information to output to parent molecule
      v_settings$ready<-T
      v_settings$norm<-input$norm
      v_settings$fact_name <- input$factor
      v_settings$fact_order<-input$fac_levels
      v_settings$normalize <- input$norm
      v_settings$P_isotopologues <- input$P_isotopologues
      v_settings$label_decimals <- input$label_decimals
      v_settings$min_lab_dist <- input$min_lab_dist
      v_settings$percent_add <- input$percent_add
      v_settings$FC_position <- input$FC_position
      v_settings$log_abund <- input$log_abund
      v_settings$col_labeling <- col_labeling()
      v_settings$circlelinetypes <- circlelinetypes()
      v_settings$circlelinecolor <- input$colour_circle
      v_settings$maxcol_facet<- input$maxcol_facet
      v_settings$include_name <- input$include_name
      v_settings$show_P<-input$show_P
      v_settings$alpha<-input$alpha
      v_settings$font<-input$font
      v_settings$otherfontsize <- input$fontsize_other
      v_settings$legendtitlesize <-input$fontsize_legend
      v_settings$cohortsize <- input$fontsize_cohort
      v_settings$include_legend <- input$include_legend


    })
    
    ###output piechart as image 
    output$pie_image <- renderImage({
      req(pies())
      # Height needs to adapt to the number of rows of the image
      height<-10*ceiling(length(unique(input$fac_levels))/input$maxcol_facet)
      

            
      tempimage<-tempfile(fileext = "png")
      ggsave(tempimage,pies(),width=24.6,height=height,
             units = "cm",device = "png")
      
      # Return a list containing the filename
      list(src = tempimage,
           contentType = 'image/png',
           width="100%",
           alt = paste0("Change input parameters to adapt this plot, use the ",
           "save options to save similar images for all compounds in the input"))
    }, deleteFile = TRUE)
    
    ###Show/hide caption
    observeEvent(input$caption_showhide, {
      if(input$caption_showhide %% 2 == 1){
        shinyjs::show(id = "caption")
      }else{
        shinyjs::hide(id = "caption")
      }
    })
    
    #render caption
    output$caption<-renderText({
      req(pies())
      
      caption<-create_caption(fact_order = v_settings$fact_order,
                     log_abund = v_settings$log_abund,
                     circlelinetypes = v_settings$circlelinetypes,
                     FC_position = v_settings$FC_position,
                     show_P = v_settings$show_P,
                     P_isotopologues = v_settings$P_isotopologues)
    })
    
    return(v_settings)
  })
}

#Wrap modules in caller function and test-------------------------------------------------------------------
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
# #example without normalized abundance
# # example_tb<-tibble(Sample = c("S1","S2","S3","S4","S1","S2","S3","S4",
# #                               "S1","S2","S3","S4"),
# #                    Cohort=c("coh1","coh1","coh2","coh3",
# #                             "coh1","coh1","coh2","coh3",
# #                             "coh1","coh1","coh2","coh3"),
# #                    Mugwort2=c("Mugw1","Mugw1","Mugw2","Mugw2",
# #                               "Mugw1","Mugw1","Mugw2","Mugw2",
# #                               "Mugw1","Mugw1","Mugw2","Mugw2"),
# #                    datatype=c("Abund","Abund","Abund","Abund",
# #                               "FracCont","FracCont","FracCont","FracCont",
# #                               "Abund","Abund","Abund","Abund"),
# #                    Result = c(40005,45858,7000000,5000000,
# #                               .5,.453,.32,1,
# #                               400005,458058,700000,500000),
# #                    Result2 = c(40005,4585824,7552123,50000,
# #                                .8,.75,.32,0.3,
# #                                400005,45805824,755123,5000))
# datatype_index<-which(colnames(example_tb)=="datatype")
# 
# example_tb<-mutate(example_tb,Sample= as.character(Sample),
#                   across(2:(datatype_index-1), as.factor))

# example_tb<-read_csv(
#   file = "~/GitHub/mec-shiny-apps/shiny-server/TraVisPies/Example_data/Standardized input/Input_Example_standardized w isotopologues.csv")

travis_vizApp<- function() {
  ui <- fluidPage(
    travis_viz_ui("viz")
  )
  
  server <- function(input, output, session) {
    travis_viz_server("viz",reactive(example_tb))
  }
  
  shinyApp(ui, server)  
}
travis_vizApp()

# test --------------------------------------------------------------------
# comptb<-example_tb[,c(2:4)]
# comptb %>% filter(!datatype=="Isotopologues") %>% 
#   mutate(Phosphoenolpyruvic_acid=as.numeric(Phosphoenolpyruvic_acid))
# 

