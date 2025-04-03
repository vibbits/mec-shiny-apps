# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)


###Summary: This code aims to produce the pie charts plots proposed in the 
# TraVis Pies: A Guide for Stable Isotope Metabolomics Interpretation Using an 
# Intuitive Visualization, based on 3 input files:metadata, raw abundance and fractional
# contribution. It is in addition possible to overlay them on a metabolic map 
# image if desired. In addition to the one-factor plots proposed in the article, 
# two-factor plotting is also available where the rows and columns of pie charts 
# in the final plot each correspond to a different factor.All of this can also
# be done for data including more than one tracer.
# The plots are generated using ggplot2 for each compound in the input data,
# several other tidyverse packages are used to get the data in the right format 
# to generate the required figures.

###Script sections:
# "Functions and libraries" loads the required libraries and defines the
# required functions
# "User input" to input file locations and select formatting options
# "Code" to execute the functions while producing intermediate tables and 
# warnings allowing to check if input data was correct or for troubleshooting 
# "Overlay pies on map" overlays the pie charts generated on a metabolic map 
# based on a separately provided .csv specifying the XY pixel coordinates for 
# overlaying

###Output
# Depending on the input options selected, the charts will either be generated
# to the IDE (tested in Rstudio) plot window, or saved in a subfolder Pie charts
# created in the input data folder. If pathway_charts requested, a second folder 
# "pie charts for pathway" will be made with more concise figures that are 
# easier to use to put on a metabolic map.

###Input file tips
#The samples considered will be based on the metadata file provided. Samples 
#in other files that are not in the metadatafile will not be used to generate
#results.

#Todo nonUI
#choose normalize or not when summarizing

# todo make it work for isotopologue pies, if FC pies fine adapt function isoslice
#to match FC slice

#make it work for 1cohort or 1replicate studies


#Todo UI
#Set variables = "None" as =NULL
#adapt UI functions to work with new setup

#todo general
#longest step with loads of data is summarize pies because of repeat over
#all groups in data, but likely not a big issue when few groups

# Functions and libraries ---------------------------------------------------------------
#libraries for functions used
library(here)    #to make r source from file location instead of magic stuff
library(vroom)   #for error messages on box
library(readxl)       #for reading excel files
library(forcats)      #for factor manipulation
library(dplyr)        #for faster.easier manipulation of data
library(tibble)       #for manipulating tibbles
library(readr)        #for writing .csv file of merged output
library(tidyr)        #for restructuring data tibbles
library(broom)        #for using regression models in dplyr pipes
library(RColorBrewer) #for generating colors
library(ggplot2)      #for generating the pie chart plots

#load functions to support the app 
source(here::here("Functions and modules/metabolomics_tibble_functions.R"))
source(here::here("Functions and modules/TraVis_Pies_functions.R"))

# User input -------------------------------------------------------------------
#Data file locations and specifications
#rawpath is windows copied folder with files, need this command to properly 
#read in R without having to modify the strings by hand
#example: rawpath<-r"(C:\User\Projects\Pie charts)"
#can also specify relative to the project folder using the here::here command

#set optional variables to NULL, functions designed to handle absence
isostring<-comparative_factor_column<-sampletype_column<-factor_column <- norm_column <- 
  tracer_column<-inputtype<-libfile<-lib_tb<-factor_levels_ordered<-col_labeling<-NULL

#set required variables to a default value that can be changed for specific 
#projects further in input
metastring<-"meta"
abundstring<-"abund"
labelstring<-"iso"
minfract_detected<-0
iso_charts<-F

# #test  excel with labeled data 1factor 1 tracer
isostring<-"_C13-0"
path<-here::here("Example_data/Experimental examples for nonUI app/Excel_1factor_1tracer")
savepath<-path
excelfile<-"example excel.xlsx"
inputpath<-paste(path,excelfile,sep = "/")
read_excel(inputpath,
           which(grepl("meta",tolower(excel_sheets(inputpath)))))
read_excel(inputpath,
           which(grepl("iso",tolower(excel_sheets(inputpath)))))%>%
  colnames()
sample_column <-"Sample"
factor_column <- "cohort"   #"None" if not present, or 1 or two element vector
factor_levels<-read_excel(inputpath,
           which(grepl("meta",tolower(excel_sheets(inputpath)))))%>%
  pull(factor_column)%>%
  unique()
factor_levels_ordered<-factor_levels[2:4]
factor_columns <- c(factor_column,comparative_factor_column)
norm_column <- "Normalisation"   #"None" if not present
sampletype_column<-"sample_type"
libfile<-"Lib_excelexample.csv"
lib_tb<-vroom::vroom(paste0(path,"/",libfile),delim = ",")

#todo test data 1-factor no replicates
# rawpath<-r"(D:\Documents\GitHub\mec-shiny-apps\shiny-server\TraVisPies\Example_data\Input noreplicate one cohort)"
# inputpath<-path<-gsub("\\\\", "/", rawpath)
# inputpath<-path<-here::here("Example_data/Input noreplicate one cohort")
# savepath<-path
# list.files(inputpath)
# metastring<-"metadata."
# abundstring<-"_RA."
# labelstring<-"_iso"
# isostring<-"_C13-label"
# loadfile_stringmatch(inputpath,metastring)%>%colnames(.)
# testload<-loadfile_stringmatch(inputpath,abundstring)
# testload<-loadfile_stringmatch(inputpath,labelstring)
# sample_column <-"Sample"
# factor_column <- "Cohort"   #"None" if not present, or 1 or two element vector


#todo test data 1-factor FC or iso input
# rawpath<-r"(C:\Users\u0134881\Documents\R\Create figures\Pie charts\Pie charts inputfiles\Pie charts 1factor)"
# inputpath<-path<-gsub("\\\\", "/", rawpath)
# inputpath<-path<-here::here("Example_data/Original input")
# savepath<-path
# mapcoordsfile<-"Pathway figure coords.csv"
# list.files(inputpath)
# metastring<-"metadata."
# abundstring<-"_RA."
# labelstring<-"_FC"
# # labelstring<-"_iso"
# isostring<-"C13-label"
# loadfile_stringmatch(inputpath,metastring)%>%colnames(.)
# testload<-loadfile_stringmatch(inputpath,abundstring)
# testload<-loadfile_stringmatch(inputpath,labelstring)
# sample_column <-"Sample"
# factor_column <- "Cohort"   #"None" if not present, or 1 or two element vector
# norm_column<-"Normalisation"

#test data 2-factor
# rawpath<-r"(D:\Documents\GitHub\mec-shiny-apps\shiny-server\TraVisPies\Example_data\Experimental examples for nonUI app\Pie charts 2factor)"
# inputpath<-path<-gsub("\\\\", "/", rawpath)
# inputpath<-path<-here::here("Example_data/Experimental examples for nonUI app/Pie charts 2factor")
# savepath<-path
# list.files(inputpath)
# metastring<-"metadata."
# abundstring<-"_RA."
# labelstring<-"_FC"
# loadfile_stringmatch(inputpath,metastring)%>%colnames(.)
# testload<-loadfile_stringmatch(inputpath,abundstring)
# testload<-loadfile_stringmatch(inputpath,labelstring)
# sample_column <-"Sample"
# factor_column <- "Time"   #"None" if not present
# comparative_factor_column <- "Condition"   #"None" if not present, factor on each level of which the first factor is compared
# norm_column <- "Normalisation"   #"None" if not present


#test data 1-factor different tracers
# rawpath<-r"(F:\Documents\Code\R\Create figures\TraVis Pies\Pie charts inputfiles\Pie charts 1 factor multitracer)"
# inputpath<-path<-gsub("\\\\", "/", rawpath)
# inputpath<-path<-here::here("Example_data/Experimental examples for nonUI app/Pie charts 1 factor multitracer")
# savepath<-path
# list.files(inputpath)
# metastring<-"metadata."
# abundstring<-"_RA."
# labelstring<-"_FC"
# loadfile_stringmatch(inputpath,metastring)%>%colnames(.)
# testload<-loadfile_stringmatch(inputpath,abundstring)
# testload<-loadfile_stringmatch(inputpath,labelstring)
# sample_column <-"Sample"
# factor_column <- "Condition"   #"None" if not present
# norm_column <- "Normalisation"   #"None" if not present
# tracer_column <-"Tracer"                #"None" if not present

#test data 2-factor different tracers
# rawpath<-r"(F:\Documents\Code\R\Create figures\TraVis Pies\Pie charts inputfiles\Pie charts 2factor multitracer)"
# inputpath<-path<-gsub("\\\\", "/", rawpath)
# inputpath<-path<-here::here("Example_data/Experimental examples for nonUI app/Pie charts 2factor multitracer")
# savepath<-path
# list.files(inputpath)
# metastring<-"2factor_multitrace_metadata.csv"
# abundstring<-"2factor_multitrace_RA.csv"
# labelstring<-"2factor_multitrace_FC.csv"
# loadfile_stringmatch(inputpath,metastring)%>%colnames(.)
# sample_column <-"Sample"
# factor_column <- "Condition"   #"None" if not present
# comparative_factor_column <- "Supplementation"   #"None" if not present, factor on each level of which the first factor is compared
# norm_column <- "Normalisation"   #"None" if not present
# tracer_column <-"Tracer"                #"None" if not present


#other data 2-factor different tracers
# rawpath<-r"(D:\Documents\Articles\Own\Sugar separation PGM1\Manuscript\Figures pies)"
# inputpath<-path<-gsub("\\\\", "/", rawpath)
# savepath<-path
# list.files(inputpath)
# metastring<-"PGM1gal_multitrace_metadata_SDC.csv"
# abundstring<-"_RA.csv"
# labelstring<-"_FC.csv"
# loadfile_stringmatch(inputpath,metastring)%>%colnames(.)
# sample_column <-"Sample"
# factor_column <- "Condition"   #"None" if not present
# comparative_factor_column <- "Supplementation"   #"None" if not present, factor on each level of which the first factor is compared
# norm_column <- "Normalisation"   #"None" if not present
# tracer_column <-"Tracer"                #"None" if not present
# col_labeling<-c("#63B2F3","#FFE699","#bfbfbf")   #colors for 2 tracers and unlabeled fraction

#excel multitracer multifactor
# isostring<-"_C13-0"
# path<-here::here("Example_data/Experimental examples for nonUI app/Excel_2factor_2tracer")
# savepath<-path
# excelfile<-"TraVis pies input PGM1.xlsx"
# inputpath<-paste(path,excelfile,sep = "/")
# read_excel(inputpath,
#            which(grepl("meta",tolower(excel_sheets(inputpath)))))
# read_excel(inputpath,
#            which(grepl("iso",tolower(excel_sheets(inputpath)))))%>%
#   colnames()
# sample_column <-"Sample"
# factor_column <- "Condition"   #"None" if not present, or 1 or two element vector
# comparative_factor_column <- "Supplementation"   #"None" if not present, factor on each level of which the first factor is compared
# norm_column <- "Normalisation"   #"None" if not present
# sampletype_column<-"sample_type"


#Miscellaneous
P_isotopologues<-T                    #leave at false, only input fraction contribution data 
log_abund<-F
detail_charts<-T                      #makes images with detail for solo use
pathway_charts<-F                       #also generate images fit for pathway
iso_charts<-T
normalize<-(length(norm_column)>0)                         #normalize abundances?
print_tables<-F                       #print generated tables to console?
compounds<-NULL                      #which compounds included; NULL => all
show_P<-T                              #show P values on pie plots

#figure appearance parameters if not specified in project
#any color input recognized by ggplot2::scale_fill_manual can be used
col_labeling<-c("#ffd966","#bfbfbf")   #colors for labeled and unlabeled fractions, NULL to use default color scheme 
maxcol_facet<-3                       #maximum amount of images horizontal
include_name<-T                        #include compound name on figure
include_legend<-T                      #include legend on figure

#axis names and fonts, load font library
#Font: set to "sans" to use standard font. For other available options run
#windowsFonts() after loading the extrafont library.If desired font not present,
#check link below on importing fonts:
#https://www.r-bloggers.com/2013/02/change-fonts-in-ggplot2-and-create-xkcd-style-graphs/
xAxLab<-""                            #xlabel, best "" for grids of pies
yAxLab<-""                            #ylabel, best "" for grids of pies
font<-"Calibri"                       #
# font<-"sans"

#figure size parameters in cm for detailed pie charts. A4 landscape is 
#recommended: width=24.6 and height 16
width<- 24.6                          
height<-16                            

#figure size parameters in cm for summary pie charts for pathway map.
#Recommended width=6.15 and height=2.81
#CURRENTLY NOT USED
# mapwidth<-6.15                          
# mapheight<-2.81    
# mapwidth<-5                           
# mapheight<-5 

#fontsizes on  detailed pie charts: Cohort names above chart, legend names above legend, all others separate)
#when width 24.6=height=16 and font= calibri
# for 2 cohorts 28 24 24 recommended
# for 3 cohorts 18 16 16 recommended
cohortsize<-18                        #text size of cohort names
legendtitlesize<-16                   #set legend font size
otherfontsize<-12                     #adapt text size of all but those above


#fontsizes on summary pie charts for pathway
#when width 6.15=height=2.81 and font= calibri
# for 2 cohorts 28 24 recommended
# for 3 cohorts 24 20 recommended
mapcohortsize<-16
mapotherfontsize<-18

#textlabel fractional contribution parameters
#FC_position sets where FC should be displayed. "center" to display in center,
#"slice" to display in labeled slice
#min_lab_dist sets minimal distance at which FC label is plotted. 0 is center
#1 is the outer circle. If distance would be smaller based on pie abundance,
#label is plotted at min_lab_dist  distance from the circle center
FC_position<-"center"  
min_lab_dist<-0.7                       
label_decimals<-1                     #amount of decimals in FC label
percent_add<-T                         #if true adds "%" to the FC label

#linetypes and colour of concentric circles
# circlelinetypes<-c(3,4,2,6)  #concentric lines from inside to out dotted, dot dash, dash, doubledash
# circlelinetypes<-c(3,3,3,3)  #all concentric lines dotted
# circlelinetypes<-c(2,2,2,2)  #all concentric lines dashed
circlelinetypes<-c(1,1,1,1)  #all concentric lines solid
# circlelinetypes<-c(0,0,0,0)  #no concentric circles

circlelinecolor<-"gray"

#Other settings
alpha<-0.7
format<-"png"


# Code merging data-------------------------------------------------------------------
debug(dolly_to_longtibble)
tb<-dolly_to_longtibble(path,inputpath,metastring = metastring,
                        abundstring = abundstring,labelstring = labelstring,
                        isostring = isostring,
                        sample_column = sample_column,
                        factor_column = factor_column,
                        comparative_factor_column =
                          comparative_factor_column,
                        factor_levels_ordered = factor_levels_ordered,
                        norm_column = norm_column,
                        tracer_column = tracer_column,
                        sampletype_column = sampletype_column,lib_tb = lib_tb,
                        savedata=T)

# print(paste(unique(tb$datatype)))
# Code pies ---------------------------------------------
#derive variables used later on
factor_symbols<-sym_or_null(factor_columns,returnlist = T)
nutrient_symbols<-NULL

#prepare tracer visualization parameters, correct when necessary
#make sure FC_position is set to slice when multiple tracer nutrients
if(length(tracer_column)>0){
  nutrient_symbols<-rlang::syms(unique(pull(tb,any_of(tracer_column))))
  tracernumber<-length(nutrient_symbols)
  if (tracernumber>1 & 
      FC_position =="center") {
    FC_position <- "slice"
    print(paste0("As multiple tracers are supplied, FC will be displayed in ",
                 "the slice."))
  }
  
  #checks if the right amount of colors is set, sets right amount of default 
  #distinctive colors (amount of tracers +1 for unlabeled fraction) if not
  if (!tracernumber == length(col_labeling)-1){
    print(paste0("Using default color scheme as for ",tracernumber," tracers ",
                 tracernumber+1," colors are required but ",
                 length(col_labeling), "were supplied."))
    if(tracernumber==1) {
      col_labeling<-c("#ffd966","#bfbfbf")
    } else {
      library(RColorBrewer)
      col_labeling<-brewer.pal(tracernumber+1,"Accent")
    }
  }
} else if("FracCont" %in% tb$datatype) {
  tracernumber<-1
  if (!tracernumber == length(col_labeling)-1){
    print(paste0("Using default color scheme as for ",tracernumber," tracers ",
                 tracernumber+1," colors are required but ",
                 length(col_labeling), "were supplied."))
    if(tracernumber==1) {
      col_labeling<-c("#ffd966","#bfbfbf")
    } else {
      library(RColorBrewer)
      col_labeling<-brewer.pal(tracernumber+1,"Accent")
    }
  }
}


#prepare summarized table with means and p values of differences
#of selected factor levels with desired factor order
sum_tb<-tb %>%
  # filter(compound=="Glucose-6-phosphate 1mox 4prop")%>%

  #Select only desired columns and filter only supported datatypes.
  #Extract data only for desired factor levels and set factor order
  prepare_piedata(factor_columns = factor_columns,
                  tracer_column = tracer_column,
                  factor_order = factor_levels_ordered)%>%
  
  #summarize data per combination of compounds, factors, tracer types and data types
  #for each comparison factor level, test differences of first factor using P value
  #if only one factor simply test differences of first factor once.
  summarise_piedata(prepare_tb,factor_column = factor_column,
                    comparative_factor_column = comparative_factor_column,
                    tracer_column = tracer_column,factor_order = 
                      factor_levels_ordered)

#obtain isotopologue slice tb for marking iso difference in FC plot and for 
#plotting iso plots if isotopologues provided
isos_calculated<-F
if(!any(grepl("iso",tolower(sum_tb$datatype)))) {
  if(exists("isoslice_tb")) rm("isoslice_tb")
} else if (length(nutrient_symbols)>1){
  if(exists("isoslice_tb")) rm("isoslice_tb")
  print(paste0("Isotopologues provided, but multiple tracer nutrients used. ",
               "This is not supported currently, isotopologue data will be ignored"))
} else {
  isoslice_tb<- sum_tb%>%
    separate_wider_delim(datatype,"_",names=c("datatype","Isotopologue"),
                         too_few = "align_start")%>%
    make_labelslices(label_type = "isotopologues",factor_columns=factor_columns,
                     label_column = "Isotopologue",
                     normalize=normalize,add_unlab = F,p_string = "P_iso",
                     label_decimals=label_decimals,
                     percent_add=percent_add,
                     FC_position=FC_position,min_lab_dist=min_lab_dist)
  
  #label factor name if any isotopologue has a significant difference, 
  #regardless of comparative factor level
  signi_iso_tb<-isoslice_tb%>%
    mutate(iso_sign_label=if_else(P_iso>=0.05|is.na(P_iso),
                                  "",
                                  "*"))%>%
    select(compound,!!!factor_symbols,iso_sign_label)%>%
    group_by(compound,!!!factor_symbols[1])%>%
    summarise(iso_sign_label=if_else(any(iso_sign_label=="*"),
                                     "*",
                                     ""))
  
  isos_calculated<-T
}

#obtain fraccont slice tb for plotting, if desired add star to cohort name if any
#isotopologues significant if isotopologues were calculated
FCslice_tb<- sum_tb%>%
  make_labelslices(label_type = "frac_con",factor_columns=factor_columns,
                   label_column = tracer_column,labelstring = "frac",
                   normalize=normalize,add_unlab = T,p_string = "P_FC",
                   label_decimals=label_decimals,
                   percent_add=percent_add,
                   FC_position=FC_position,min_lab_dist=min_lab_dist)%>%
  {
    if(P_isotopologues & isos_calculated) {
      left_join(.,signi_iso_tb) %>%
        mutate(iso_sign_label=if_else(is.na(iso_sign_label),
                                      "",
                                      iso_sign_label),
               !!factor_symbols[[1]]:=paste0(!!factor_symbols[[1]],iso_sign_label))%>%
        select(-iso_sign_label)
    } else {
      .
    }
  }

make_piechart(FCslice_tb,
              factor_columns = factor_columns,
              tracer_column = tracer_column,
              log_abund=log_abund,
              circlelinecolor = circlelinecolor,
              selected_compound= unique(FCslice_tb$compound)[2],
              circlelinetypes = circlelinetypes,
              maxcol_facet = maxcol_facet,
              include_name = include_name,col_labeling = col_labeling,
              alpha=alpha,font=font,otherfontsize = otherfontsize,
              legendtitlesize =legendtitlesize,
              cohortsize = cohortsize,include_legend = include_legend,
              show_P=show_P)

#todo test color selection, now only too few
if (exists("isoslice_tb") & iso_charts){
  
  make_piechart(isoslice_tb,
                factor_columns = factor_columns,
                tracer_column = "Isotopologue",
                log_abund=log_abund,
                circlelinecolor = circlelinecolor,
                selected_compound= "L-Isoleucine",
                circlelinetypes = circlelinetypes,
                maxcol_facet = maxcol_facet,
                include_name = include_name,col_labeling = col_labeling,
                alpha=alpha,font=font,otherfontsize = otherfontsize,
                legendtitlesize =legendtitlesize,
                cohortsize = cohortsize,include_legend = include_legend,
                show_P=show_P)
}

#loop over each compound in input tibble to create 
if(detail_charts|pathway_charts) {
  # debug(generate_pies)
  generate_pies(FCslice_tb,
                compound_col="compound",detail_charts=detail_charts,
                pathway_charts=pathway_charts,savepath=savepath,
                normalize=normalize,
                factor_columns = factor_columns,
                tracer_column = tracer_column,
                log_abund=log_abund,
                circlelinecolor = circlelinecolor,
                circlelinetypes = circlelinetypes,
                maxcol_facet = maxcol_facet,
                include_name = include_name,col_labeling = col_labeling,
                alpha=alpha,font=font,otherfontsize = otherfontsize,
                legendtitlesize =legendtitlesize,
                cohortsize = cohortsize,include_legend = include_legend,
                show_P=show_P,format=format, width=width,height=height)
}

if (iso_charts) {
  generate_pies(isoslice_tb,
                compound_col="compound",detail_charts=detail_charts,
                pathway_charts=F,savepath=savepath,subfolder="pies isos",
                normalize,
                factor_columns = factor_columns,
                tracer_column = "Isotopologue",
                log_abund=log_abund,
                circlelinecolor = circlelinecolor,
                circlelinetypes = circlelinetypes,
                maxcol_facet = maxcol_facet,
                include_name = include_name,col_labeling = col_labeling,
                alpha=alpha,font=font,otherfontsize = otherfontsize,
                legendtitlesize =legendtitlesize,
                cohortsize = cohortsize,include_legend = include_legend,
                show_P=show_P,format=format,width=width,height=height)
}

#save caption as text file
fileConn<-file(paste0(savepath,"/caption.txt"))
writeLines(
  create_caption(factor_order = factor_levels_ordered,log_abund = log_abund,
                 circlelinetypes = circlelinetypes,FC_position = FC_position,
                 show_P = show_P,P_isotopologues = P_isotopologues),
  fileConn)
close(fileConn)

# Function pies ---------------------------------------------
plot_fromtibble(tb,charttype="FC",compound_col="compound",
                selected_compound=unique(tb$compound)[6],
                factor_columns=factor_columns,
                tracer_column=tracer_column,
                log_abund=log_abund,
                circlelinecolor=circlelinecolor,
                circlelinetypes=circlelinetypes,
                maxcol_facet=maxcol_facet,
                include_name=include_name,
                col_labeling=col_labeling,
                alpha=alpha,
                otherfontsize=otherfontsize,
                font=font,
                legendtitlesize=legendtitlesize,
                cohortsize=cohortsize,
                include_legend=include_legend,
                show_P=show_P)

# debug(plot_fromtibble)
plot_fromtibble(tb,charttype="iso",compound_col="compound",
                selected_compound=unique(tb$compound)[6],
                factor_columns=factor_columns,
                tracer_column=tracer_column,
                log_abund=log_abund,
                circlelinecolor=circlelinecolor,
                circlelinetypes=circlelinetypes,
                maxcol_facet=maxcol_facet,
                include_name=include_name,
                col_labeling=col_labeling,
                alpha=alpha,
                otherfontsize=otherfontsize,
                font=font,
                legendtitlesize=legendtitlesize,
                cohortsize=cohortsize,
                include_legend=include_legend,
                show_P=show_P)


    
mergedtibble_to_pies(tb,compound_col="compound",
                     detail_charts=detail_charts,
                     pathway_charts=pathway_charts,
                     iso_charts=iso_charts,
                     savepath=savepath,
                     normalize=normalize,
                     factor_columns=factor_columns,
                     tracer_column=tracer_column,
                     log_abund=log_abund,
                     circlelinecolor=circlelinecolor,
                     circlelinetypes=circlelinetypes,
                     maxcol_facet=maxcol_facet,
                     include_name=include_name,
                     col_labeling=col_labeling,
                     alpha=alpha,
                     otherfontsize=otherfontsize,
                     font=font,
                     legendtitlesize=legendtitlesize,
                     cohortsize=cohortsize,
                     include_legend=include_legend,
                     format=format,
                     show_P=show_P,
                     width=width,height=height) 

#todo add code to generate caption
# create_caption<-function(factor_order,log_abund,circlelinetypes,FC_position,show_P,
#                          P_isotopologues) {
# Overlay pies on map --------------------------------------------------------
#how to assign coordinates: get bitmap format empty map, eg. import empty map template in r then export as png, use this as base empty map
#if using powerpoint have to make image with bounds larger than what will be needed, then save as bitmap in this step!!!
#reason: if any background layer is transparent, you get problems. Possibly can be solved in R as well
#import empty map png in Inkscape, in document properties set drawing to fit borders to png and set units to px 
#then import image generated by generatie pies code above and select default import resolution
#then note coordinates from top left corner of image in a .csv file linking them to the metabolite name used in the analysis

library(magick)

#read in empty path and figure to path mapping in subfolder
plotfilepath<-paste0(path,"/Pie charts pathway/")         #file with pie charts meant for pathway
plotfilepath<-paste0(path,"/Pie charts pathway gluc/")         #file with pie charts meant for pathway
mapfile<-""
figurepath<-paste0(path,"/pathway/")
pathway.img <- image_read(paste0(figurepath,mapfile))
fig.coords<-read.csv(paste0(figurepath,mapcoordsfile)) %>%
  filter(Compound %in% compounds)                         #drops entries on map that don't have their name among the pies (or that don't have a compound name, aka they shouldn't be plotted)

#plot figure per figure
i<-1
for (i in 1:nrow(fig.coords)) {
  #obtain right file and read image
  if (normalize) {
    plotfilename<-paste0(plotfilepath,"pies normalized ",fig.coords$Compound[i],".png")
  } else {
    plotfilename<-paste0(plotfilepath,"pies ",fig.coords$Compound[i],".png")
  }
  pie.img<-image_read(plotfilename)
  #plot file over 
  pathway.img<-image_composite(pathway.img,pie.img,offset=
                                 geometry_point(fig.coords$Xoffset[i],fig.coords$Yoffset[i]))
  
}

# print(pie.img)  #check generated image in r studio, very slow!
# 

#write image
image_write(pathway.img, paste0(figurepath,"output pathway.png"), format = "png")
gc()              #needed to reproducibly release image objects, otherwise they will not correctly generate output if this function is run again in quick succession

