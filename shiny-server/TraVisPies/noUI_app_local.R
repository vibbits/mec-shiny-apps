# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

#note: this non-shiny code does not yet support input with isotopologue data
#please use fraction contribution input instead

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

# Functions and libraries ---------------------------------------------------------------
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

# User input -------------------------------------------------------------------
#Data file locations and specifications
#rawpath is windows copied folder with files, need this command to properly 
#read in R without having to modify the strings by hand
#example: rawpath<-r"(C:\User\Projects\Pie charts)"
#can also specify relative to the project folder using the here::here command

#test data 1-factor
# rawpath<-r"(C:\Users\u0134881\Documents\R\Create figures\Pie charts\Pie charts inputfiles\Pie charts 1factor)"
# path<-gsub("\\\\", "/", rawpath)
# savepath<-path
# metadatafile<-"Input_Example_metadata.csv"
# abundancefile<-"Input_Example_RA.csv"
# fracconfile<-"Input_Example_FC.csv"
# mapcoordsfile<-"Pathway figure coords.csv"
# read_csv_clean(file=paste(path,metadatafile,sep = "/"),
#                remove_empty = T)%>%colnames()
# sample_column <-"Sample"
# factor_column <- "Cohort"   #"None" if not present, or 1 or two element vector
# norm_column <- "None"   #"None" if not present
# tracer_column <-"None"                #"None" if not present

#test data 1-factor iso input
# rawpath<-r"(F:\Documents\Code\Github\mec-shiny-apps\shiny-server\TraVisPies\Example_data\Crashes)"
# path<-gsub("\\\\", "/", rawpath)
path<-here::here("Example_data/Crashes")
savepath<-path
metadatafile<-"Input_MCF1711_metadata.csv"
abundancefile<-"Input_MCF1711_RA.csv"
tracerfile<-"Input_MCF1711_isotopologues.csv"
norm_column<-"None"
(meta_tb<-read_csv_clean(paste0(path,"/",metadatafile),remove_empty = T,
                         remove_rowempty = T))
sample_column<-colnames(meta_tb)[1]
factor_column<-colnames(meta_tb)[2]
tracer_column <-"None"                #"None" if not present


# meta_formatted_tb<-format_metadata(meta_tb = meta_tb,
#                                    sample_column = sample_column,
#                                    factor_column = cohort_column,
#                                    norm_column = norm_column)
# 
# abund_tb<-read_csv_clean(paste0(path,"/",abundancefile),remove_empty = T,
#                          remove_rowempty = T)
# compounds<-colnames(abund_tb[-which(colnames(abund_tb)==sample_column)])


#test data 2-factor
# rawpath<-r"(C:\Users\u0134881\Documents\R\Create figures\Pie charts\Pie charts inputfiles\Pie charts 2factor)"
# path<-gsub("\\\\", "/", rawpath)
# path<-here::here("Example_data/Other examples for nonUI app/Pie charts 2factor")
# savepath<-path
# metadatafile<-"2factorpies_metadata.csv"
# abundancefile<-"2factorpies_RA.csv"
# fracconfile<-"2factorpies_FC.csv"
# read_csv_clean(file=paste(path,metadatafile,sep = "/"),
#                remove_empty = T)%>%colnames()
# sample_column <-"Sample"
# factor_column <- c("Time","Condition")   #"None" if not present, or 1 or two element vector
# norm_column <- "Normalisation"   #"None" if not present
# tracer_column <-"None"                #"None" if not present


#test data 1-factor different tracers
# rawpath<-r"(F:\Documents\Code\R\Create figures\TraVis Pies\Pie charts inputfiles\Pie charts 1 factor multitracer)"
# path<-gsub("\\\\", "/", rawpath)
# path<-here::here("Example_data/Other examples for nonUI app/Pie charts 1 factor multitracer")
# savepath<-path
# metadatafile<-"MCF001748,74_multitrace_metadata.csv"
# abundancefile<-"MCF001748,74_multitrace_RA.csv"
# fracconfile<-"MCF001748,74_multitrace_FC.csv"
# read_csv_clean(file=paste(path,metadatafile,sep = "/"),
#                remove_empty = T)%>%colnames()
# read_csv_clean(file=paste(path,fracconfile<-"MCF001748,74_multitrace_FC.csv"
# ,sep = "/"),
#                remove_empty = T)%>%colnames()
# sample_column <-"Sample"
# factor_column <- c("Condition")   #"None" if not present, or 1 or two element vector
# norm_column <- "Normalisation"   #"None" if not present
# tracer_column <-"Tracer"                #"None" if not present

#test data 2-factor different tracers
# rawpath<-r"(F:\Documents\Code\R\Create figures\TraVis Pies\Pie charts inputfiles\Pie charts 2factor multitracer)"
# path<-gsub("\\\\", "/", rawpath)
# path<-here::here("Example_data/Other examples for nonUI app/Pie charts 2factor multitracer")
# savepath<-path
# metadatafile<-"2factor_multitrace_metadata.csv"
# abundancefile<-"2factor_multitrace_RA.csv"
# fracconfile<-"2factor_multitrace_FC.csv"
# read_csv_clean(file=paste(path,metadatafile,sep = "/"),
#                remove_empty = T)%>%colnames()
# read_csv_clean(file=paste(path,fracconfile,sep = "/"),
#                remove_empty = T)%>%colnames()
# sample_column <-"Sample"
# factor_column <- c("Condition","Supplementation")   #"None" if not present, or 1 or two element vector
# norm_column <- "Normalisation"   #"None" if not present
# tracer_column <-"Tracer"                #"None" if not present

#Miscellaneous
P_isotopologues<-F                     #leave at false, only input fraction contribution data 
log_abund<-F
detail_charts<-T                      #makes images with detail for solo use
pathway_charts<-F                        #also generate images fit for pathway 
save_chart<-T                         #save chart images? If false plots in ID
normalize<-(!norm_column=="None")                         #normalize abundances?
print_tables<-F                       #print generated tables to console?
compounds<-NULL                      #which compounds included; NULL => all
show_P<-T                              #show P values on pie plots
compounds<-"Glucose 1,6-bisphosphate"  #TODO:remove

#factors and factor level order
#fact.invert False makes first factor X factor and second Y for grid, 
#True is opposite
#factX.levels and factY.levels can be used to set the factor level order
#in the pie chart if the order in the input file is not satisfactory.
#important as the first level will be used as a reference

fact.invert<-F 
factX.levels<-NULL #only used in two factor analysis. 
# factX.levels<-c("NT","25mM 2DG", "10uM AMA")                    #character vector with order of factor levels in X factor, set NULL to take order in file
factY.levels<-NULL                    #character vector with order of factor levels in Y factor, set NULL to take order in file
# factY.levels<-c("Reference","Treatment1","Treatment2")                    #character vector with order of factor levels in Y factor, set NULL to take order in file

#figure appearance parameters
#any color input recognized by ggplot2::scale_fill_manual can be used
colLabeling<-c("#bfbfbf","#ffd966")   #colors for labeled and unlabeled fraction 
# colLabeling<-c("#bfbfbf","#ffd966","lightblue")   #colors for 2 tracers and unlabeled fraction
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
otherfontsize<-16                     #adapt text size of all but those above


#fontsizes on summary pie charts for pathway
#when width 6.15=height=2.81 and font= calibri
# for 2 cohorts 28 24 recommended
# for 3 cohorts 24 20 recommended
mapcohortsize<-16
mapotherfontsize<-18

#textlabel fractional contribution parameters
#FCposition sets where FC should be displayed. "center" to display in center,
#"slice" to display in labeled slice
#minLabDist sets minimal distance at which FC label is plotted. 0 is center
#1 is the outer circle. If distance would be smaller based on pie abundance,
#label is plotted at minLabDist  distance from the circle center
FCposition<-"center"  
minLabDist<-0.7                       
labelDecimals<-1                     #amount of decimals in FC label
percentAdd<-T                         #if true adds "%" to the FC label

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


# Code-------------------------------------------------------------------
#derive settings from input
if(length(factor_column)==1){
  twofactor=F
} else if(length(factor_column)==2) {
  twofactor=T
} else {
  stop("Factor should be set to 'None' if not present, or be a 1 or 2 element vector")
}

#If no tracer column, make dummy column and set tracer column name to Tracer
if(tracer_column=="None") {
  tracer_column<-"Labeling"
}

#input metadata and abundance data, put in right format for following functions
meta_formatted_tb<-read_csv_clean(file=paste(path,metadatafile,sep = "/"),
                                  remove_empty = T) %>%
  format_metadata(sample_column = sample_column,
                  factor_column = factor_column,
                  norm_column = norm_column,
                  tracer_column=tracer_column)

abund_tb<-read_csv_clean(paste(path,abundancefile,sep = "/"),remove_empty = T)

#read isotopologue or fractional contribution data. Set empty isotopologue tibble
#if no isotopologue data supplied
if(grepl("iso",tracerfile)) {
  iso_tb<-extract_col_isotopologues(
    read_csv_clean(paste0(path,"/",tracerfile),remove_empty = T,
                   remove_rowempty = T),
    iso_suffix_sep = "_")
  frac_tb<-calculate_FC(iso_tb)
} else {
  frac_tb<-read_csv_clean(paste0(path,"/",tracerfile),remove_empty = T,
                          remove_rowempty = T)
  iso_tb<-NULL
  merged<-merge_input(meta_tb = meta_formatted_tb,
                      abund_tb = abund_tb,
                      frac_tb = frac_tb,
                      sample_col = sample_column,
                      compounds = compounds)
}

#Per compound adapt FC's below 0 (artefacts due to natural abundance
#correction) to be positive to avoid problems with the visualisations
#later on.
# for (i in (2:ncol(frac_tb))) {
#   if (any(frac_tb[,i]<0)) {
#     FCs<-pull(frac_tb[,i])
#     FCs[which(FCs<0)]<-FCs[which(FCs<0)]-min(FCs[which(FCs<0)])
#     frac_tb[,i]<-FCs
#   }
# }

#make sure FCposition is set to slice when multiple tracer nutrients
if (length(unique(pull(meta_formatted_tb[,tracer_column])))>1 & 
    FCposition =="center") {
  FCposition <- "slice"
}

#checks if the right amount of colors is set, sets right amount of default 
#distinctive colors (amount of tracers +1 for unlabeled fraction) if not
if (!length(unique(pull(meta_formatted_tb[,tracer_column]))) == length(colLabeling)-1){
  if(length(unique(pull(meta_formatted_tb[,tracer_column])))==1) {
    colLabeling<-c("#bfbfbf","#ffd966")
  } else {
    library(RColorBrewer)
    colLabeling<-brewer.pal(length(unique(pull(meta_formatted_tb[,tracer_column])))+1,"Accent")
  }
}


#generate error or warning messages if any
check_output<-check_samples_compounds(
  meta_tb = meta_formatted_tb,
  abund_tb = abund_tb,
  frac_tb = frac_tb,
  sample_column = sample_column,
  norm_column = norm_column)

if (check_output$error) {
  validate(check_output$message)
} else {
  outputtext<-check_output$message
}

#Getcompound names from input
if (length(compounds)<length(colnames(abund_tb)[2:ncol(abund_tb)])) {
  compounds<-colnames(abund_tb)[2:ncol(abund_tb)]
}

#merge all input into one table, get compounds and factor orders
tb<-merge_input(meta_tb = meta_formatted_tb,abund_tb = abund_tb,frac_tb = frac_tb,
                compounds=compounds, sample_col = sample_column,
                iso_tb=iso_tb)

compounds_updated<-colnames(tb)[which(!colnames(tb)%in%
                                        c(colnames(meta_formatted_tb),
                                          "datatype"))]

compounds_updated<-"Sorbitol"  #for checking specific compound only

#make list of factor orders per factor to allow multiple factors
fact_order<-list(NULL)
for (i in 1:length(factor_column)) {
  fact_order[[i]]<-unique(pull(tb,!!factor_column[i]))
}

#generate figures per compound with specified settings and save if requested
generate_multiple_pies(tb,compounds=compounds_updated,
                       detail_charts=detail_charts,
                       pathway_charts=pathway_charts,
                       savepath=savepath,
                       normalize=normalize,
                       fact_name=factor_column,
                       tracer_column=tracer_column,
                       fact_order=fact_order, 
                       P_isotopologues=P_isotopologues,
                       log_abund=log_abund,
                       label_decimals=labelDecimals,
                       percent_add=percentAdd,
                       FC_position=FCposition,
                       min_lab_dist=minLabDist,
                       circlelinecolor=circlelinecolor,
                       circlelinetypes=circlelinetypes,
                       maxcol_facet=maxcol_facet,
                       include_name=include_name,
                       show_P=show_P,
                       col_labeling=colLabeling,
                       alpha=alpha,
                       otherfontsize=otherfontsize,
                       font=font,
                       legendtitlesize=legendtitlesize,
                       cohortsize=cohortsize,
                       include_legend=include_legend,
                       format=format,
                       mapotherfontsize=mapotherfontsize,
                       mapcohortsize=mapcohortsize)

#Todo nonshiny code check FC p value gal6phos = 0.22 without fraction??

debug(corFC_addUnlab)
undebug(corFC_addUnlab)

debug(prepare_slicedata)
undebug(prepare_slicedata)

debug(add_FClabels)
undebug(add_FClabels)

debug(summarize_compounddata)
undebug(summarize_compounddata)

debug(summarize_addP)
undebug(summarize_addP)

debug(generate_pie)
undebug(generate_pie)

debug(make_piechart)
undebug(make_piechart)

debug(obtain_compounddata)
undebug(obtain_compounddata)


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


#make barcharts-------------------------------------------------------------------
save_chart<-T                        #save chart images? If false plots in ID

# convert raw path to R usable path, input files and check input
path<-gsub("\\\\", "/", rawpath)
meta_tb<-read_metacsv_clean(file=paste(path,metadatafile,sep = "/"),
                            normalize = normalize,twofactor = twofactor,
                            fact.invert = fact.invert,
                            factX.levels = factX.levels,
                            factY.levels = factY.levels)
abund_tb<-read_csv_clean(paste(path,abundancefile,sep = "/"),remove_empty = T)
fraccon_tb<-read_csv_clean(paste(path,fracconfile,sep = "/"),remove_empty = T)

check_input(meta_tb,abund_tb,fraccon_tb,FCposition = FCposition,
            colLabeling = colLabeling)


#Get factor and compound names from input
fact.names<-get_factornames(tb=meta_tb,twofactor = twofactor)
if (length(compounds)<length(colnames(abund_tb)[2:ncol(abund_tb)])) {
  compounds<-colnames(abund_tb)[2:ncol(abund_tb)]
  
}

#merge all input into one table
tb<-merge_input(meta_tb = meta_tb,abund_tb = abund_tb,fraccon_tb = fraccon_tb,
                compounds=compounds, normalize = normalize)
compounds<-colnames(tb)[(ncol(meta_tb)+1):(ncol(tb)-1)]
# compounds<-c("Hexose")

#initialize messages variable, convert raw input filepath to R path
messages<-NULL

#check input
if (any(!compounds %in% colnames(tb))) {
  stop(paste0("Some requested metabolite names in 'compounds' are not ",
              "among the abundance column names. Make sure all requested",
              " compounds appear with the same name in the input files"))
}


#loop over each compound in input tibble
for (compound in compounds) {
  print(paste0("Processing ",compound))
  
  #prepare filename if saving required
  if (save_chart) {
    if (normalize) {
      plotfilename<-paste0("pies normalized ",compound,".png")
    } else {
      plotfilename<-paste0("pies ",compound,".png")
    }
  }
  
  #get table with only measured compound data, then a table summarizing
  #derived means and p values per cohort for abundance and one for fractional
  #contribution, then put together table with inputformat for pie function
  print(paste0(compound,": Extracting compound data"))
  compound_tb<-obtain_compounddata(tb,compound,fact.names)
  if (print_tables) print(compound_tb)
  
  #rename P variable for fusing with FC table that also has P column, and 
  #compound variable to Abund as all values are abundances
  print(paste0(compound,": Summarizing abundance data"))
  sum_tb_ab<-summarize_compounddata(filter(compound_tb,datatype=="Abund"),
                                    compound,fact.names)%>%
    rename(Abund=compound,P.RA=P)
  if (print_tables) print(sum_tb_ab)
  
  #rename P variable for fusing with abundance table that also has P column, 
  # and compound variable to FracCont as all values are fractional
  #contributions. Also adds entry for unlabeled fraction
  print(paste0(compound,": Summarizing fractional contribution data"))
  sum_tb_FC<-summarize_compounddata(filter(compound_tb,datatype=="FracCont"),
                                    compound,fact.names)%>%
    add_unlabeled_sum(compound,fact.names = fact.names)%>%
    rename(FracCont=compound,P.FC=P)
  if (print_tables) print(sum_tb_FC)
  
  #gather abundance and fraccont data together in one input table with one 
  #entry per pie slice (per combination cohort and labeling origin) with all 
  #other required info including labels for plotting function
  print(paste0(compound,": combining abundance and fractional contribution",
               "data and adding info required for plotting pie slices and",
               "labels"))
  slice_tb<-prepare_slicedata(compound_tb,sum_tb_FC,fact.names = fact.names,
                              compound=compound,labelDecimals = labelDecimals,
                              minLabDist = minLabDist,percentAdd = percentAdd,
                              FCposition = FCposition)
  
  if (print_tables) print(slice_tb)
  
  savepath<-path
  #plot based on information in slice table
  print(paste0(compound,": Building chart"))
  
  
  #abundance in area through width, and fraction in heigth
  # widths<-sumgluctb$NormAbund*max(sumgluctb$Time)/length(unique(sumgluctb$Time)) #maximal width for proportional area of bars
  widths<-slice_tb$Abund
  
  barcohortsize<-10
  #plot detailed pie chart based on information in slice table
  print(paste0(compound,": Building detailed pie chart"))
  chart<-slice_tb %>% ggplot(aes(x = Cohort, y = Fraction, fill = Labeling,width = Abund)) + 
    geom_bar(stat = "identity", position = "fill", color = "black",size=0.5) + #make basic rectangle plot, fill causes height to be standardized so ratio can be inspected easily
    geom_hline(yintercept=seq(0,1,0.2),color="white",alpha=0.5)+  #add shape of biggest bar area as dottet rectangle
    geom_hline(yintercept=c(0,1),linetype="dotted",size=1)+  #add shape of biggest bar area as dottet rectangle
    geom_vline(xintercept=c((1:length(unique(slice_tb$Cohort)))-max(widths)/2,length(unique(slice_tb$Cohort))+max(widths)/2),linetype="dotted",size=1)+ #add gridlines to later mark size fraction
    scale_x_discrete(breaks = unique(slice_tb$Cohort)) +  #set time indications in center of each bar and on correct time axis
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    # facet_wrap(~Cohort,ncol=1) + #put in facet plots below each other
    theme_bw(base_size = barcohortsize) +#Change plots to black on white
    theme(panel.grid= element_blank(),  #remove minor grid lines of all axes
          strip.background = element_rect(fill = NA, colour = NA),
          axis.ticks = element_blank(),
          plot.title = element_text(size = barcohortsize, face = "bold"),
          legend.title = element_text(size = legendtitlesize),
          strip.text = element_text(size = barcohortsize))+
    scale_fill_manual(values=colLabeling[length(colLabeling):1],guide=guide_legend(reverse=T))+
    labs(x=NULL, y=NULL)
  
  #removes legend if desired
  chart <-chart + theme(legend.position = "none")  
  
  
  #save detailed chart if required or print to rstudio plot
  print(paste0(compound,": Plotting or saving pie chart"))
  if (save_chart) {
    #set folder path to save pie charts if saving requested
    plotfilefolder<-paste0(savepath,"/Pie charts pathway/")
    plotfilepath<-paste0(plotfilefolder,plotfilename)
    if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder))
    ggsave(plotfilepath,chart,width=mapwidth,height=mapheight,units = "cm",
           device = "png")
  } else {
    print(chart)
  }
  
  
  print(c(messages,"Finished"))
}

#make XYplots-------------------------------------------------------------------
save_chart<-T                        #save chart images? If false plots in ID
mapcoordsfile<-"MCF968 glycoshort  coords - scatter.csv"


# convert raw path to R usable path, input files and check input
path<-gsub("\\\\", "/", rawpath)
meta_tb<-read_metacsv_clean(file=paste(path,metadatafile,sep = "/"),
                            normalize = normalize,twofactor = twofactor,
                            fact.invert = fact.invert,
                            factX.levels = factX.levels,
                            factY.levels = factY.levels)
abund_tb<-read_csv_clean(paste(path,abundancefile,sep = "/"),remove_empty = T)
fraccon_tb<-read_csv_clean(paste(path,fracconfile,sep = "/"),remove_empty = T)

check_input(meta_tb,abund_tb,fraccon_tb,FCposition = FCposition,
            colLabeling = colLabeling)


#Get factor and compound names from input
fact.names<-get_factornames(tb=meta_tb,twofactor = twofactor)
if (length(compounds)<length(colnames(abund_tb)[2:ncol(abund_tb)])) {
  compounds<-colnames(abund_tb)[2:ncol(abund_tb)]
  
}

#merge all input into one table
tb<-merge_input(meta_tb = meta_tb,abund_tb = abund_tb,fraccon_tb = fraccon_tb,
                compounds=compounds, normalize = normalize)
compounds<-colnames(tb)[(ncol(meta_tb)+1):(ncol(tb)-1)]
# compounds<-c("Hexose")

#initialize messages variable, convert raw input filepath to R path
messages<-NULL

#check input
if (any(!compounds %in% colnames(tb))) {
  stop(paste0("Some requested metabolite names in 'compounds' are not ",
              "among the abundance column names. Make sure all requested",
              " compounds appear with the same name in the input files"))
}


#loop over each compound in input tibble
for (compound in compounds) {
  print(paste0("Processing ",compound))
  
  #prepare filename if saving required
  if (save_chart) {
    if (normalize) {
      plotfilename<-paste0("pies normalized ",compound,".png")
    } else {
      plotfilename<-paste0("pies ",compound,".png")
    }
  }
  
  #get table with only measured compound data, then a table summarizing
  #derived means and p values per cohort for abundance and one for fractional
  #contribution, then put together table with inputformat for pie function
  print(paste0(compound,": Extracting compound data"))
  compound_tb<-obtain_compounddata(tb,compound,fact.names)
  if (print_tables) print(compound_tb)
  
  #rename P variable for fusing with FC table that also has P column, and 
  #compound variable to Abund as all values are abundances
  print(paste0(compound,": Summarizing abundance data"))
  sum_tb_ab<-summarize_compounddata(filter(compound_tb,datatype=="Abund"),
                                    compound,fact.names)%>%
    rename(Abund=compound,P.RA=P)
  if (print_tables) print(sum_tb_ab)
  
  #rename P variable for fusing with abundance table that also has P column, 
  # and compound variable to FracCont as all values are fractional
  #contributions. Also adds entry for unlabeled fraction
  print(paste0(compound,": Summarizing fractional contribution data"))
  sum_tb_FC<-summarize_compounddata(filter(compound_tb,datatype=="FracCont"),
                                    compound,fact.names)%>%
    add_unlabeled_sum(compound,fact.names = fact.names)%>%
    rename(FracCont=compound,P.FC=P)
  if (print_tables) print(sum_tb_FC)
  
  #gather abundance and fraccont data together in one input table with one 
  #entry per pie slice (per combination cohort and labeling origin) with all 
  #other required info including labels for plotting function
  print(paste0(compound,": combining abundance and fractional contribution",
               "data and adding info required for plotting pie slices and",
               "labels"))
  slice_tb<-prepare_slicedata(compound_tb,sum_tb_FC,fact.names = fact.names,
                              compound=compound,labelDecimals = labelDecimals,
                              minLabDist = minLabDist,percentAdd = percentAdd,
                              FCposition = FCposition)
  
  if (print_tables) print(slice_tb)
  
  savepath<-path
  #plot based on information in slice table
  print(paste0(compound,": Building chart"))
  
  
  #abundance in area through width, and fraction in heigth
  # widths<-sumgluctb$NormAbund*max(sumgluctb$Time)/length(unique(sumgluctb$Time)) #maximal width for proportional area of bars
  widths<-slice_tb$Abund
  
  #plot detailed pie chart based on information in slice table
  print(paste0(compound,"Building detailed pie chart"))
  chart<- slice_tb %>% filter(Labeling=="Labeled") %>%
    ggplot(aes(x = FracCont, y = Abund, shape=Cohort,color=Cohort)) + 
    geom_hline(yintercept=seq(0,1,0.2),color="gray",alpha=0.5,size=0.5)+
    geom_vline(xintercept=seq(0,1,0.2),color="gray",alpha=0.5,size=0.5)+
    geom_hline(yintercept=c(0,1),color="black",alpha=0.5,size=0.5)+
    geom_vline(xintercept=c(0,1),color="black",alpha=0.5,size=0.5)+
    geom_point(size=2)+
    # lims(x=c(-0.2,1.2),y=c(-0.2,1.2))+
    scale_x_continuous(breaks = seq(0,1,0.2),limits = c(-0.05,1.05)) +
    scale_y_continuous(breaks = seq(0,1,0.2),limits = c(-0.05,1.05)) +
    theme_bw(base_size = mapotherfontsize*1.5) +#Change plots to black on white
    theme(panel.grid= element_blank(),  #remove minor grid lines of all axes
          strip.background = element_rect(fill = NA, colour = NA),
          axis.ticks = element_blank(),
          plot.title = element_text(size = mapcohortsize*1.5, face = "bold"),
          legend.title = element_text(size = legendtitlesize),
          strip.text = element_text(size = mapcohortsize))+
    labs(x="Fractional contribution", y="Abundance")
  #removes legend if desired
  chart <-chart + theme(legend.position = "none") 
  
  mapwidth<-2.81*1.75                          
  mapheight<-2.81*1.25
  
  #save detailed chart if required or print to rstudio plot
  print(paste0(compound,": Plotting or saving pie chart"))
  if (save_chart) {
    #set folder path to save pie charts if saving requested
    plotfilefolder<-paste0(savepath,"/Pie charts pathway/")
    plotfilepath<-paste0(plotfilefolder,plotfilename)
    if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder))
    ggsave(plotfilepath,chart,width=mapwidth,height=mapheight,units = "cm",
           device = "png")
  } else {
    print(chart)
  }
  
  
  print(c(messages,"Finished"))
}

#make basic pie charts-------------------------------------------------------------------
save_chart<-T                        #save chart images? If false plots in ID


# convert raw path to R usable path, input files and check input
path<-gsub("\\\\", "/", rawpath)
meta_tb<-read_metacsv_clean(file=paste(path,metadatafile,sep = "/"),
                            normalize = normalize,twofactor = twofactor,
                            fact.invert = fact.invert,
                            factX.levels = factX.levels,
                            factY.levels = factY.levels)
abund_tb<-read_csv_clean(paste(path,abundancefile,sep = "/"),remove_empty = T)
fraccon_tb<-read_csv_clean(paste(path,fracconfile,sep = "/"),remove_empty = T)

check_input(meta_tb,abund_tb,fraccon_tb,FCposition = FCposition,
            colLabeling = colLabeling)


#Get factor and compound names from input
fact.names<-get_factornames(tb=meta_tb,twofactor = twofactor)
if (length(compounds)<length(colnames(abund_tb)[2:ncol(abund_tb)])) {
  compounds<-colnames(abund_tb)[2:ncol(abund_tb)]
  
}

#merge all input into one table
tb<-merge_input(meta_tb = meta_tb,abund_tb = abund_tb,fraccon_tb = fraccon_tb,
                compounds=compounds, normalize = normalize)
compounds<-colnames(tb)[(ncol(meta_tb)+1):(ncol(tb)-1)]
# compounds<-c("Hexose")

#initialize messages variable, convert raw input filepath to R path
messages<-NULL

#check input
if (any(!compounds %in% colnames(tb))) {
  stop(paste0("Some requested metabolite names in 'compounds' are not ",
              "among the abundance column names. Make sure all requested",
              " compounds appear with the same name in the input files"))
}


#loop over each compound in input tibble
for (compound in compounds) {
  print(paste0("Processing ",compound))
  
  #prepare filename if saving required
  if (save_chart) {
    if (normalize) {
      plotfilename<-paste0("pies normalized ",compound,".png")
    } else {
      plotfilename<-paste0("pies ",compound,".png")
    }
  }
  
  #get table with only measured compound data, then a table summarizing
  #derived means and p values per cohort for abundance and one for fractional
  #contribution, then put together table with inputformat for pie function
  print(paste0(compound,": Extracting compound data"))
  compound_tb<-obtain_compounddata(tb,compound,fact.names)
  if (print_tables) print(compound_tb)
  
  #rename P variable for fusing with FC table that also has P column, and 
  #compound variable to Abund as all values are abundances
  print(paste0(compound,": Summarizing abundance data"))
  sum_tb_ab<-summarize_compounddata(filter(compound_tb,datatype=="Abund"),
                                    compound,fact.names)%>%
    rename(Abund=compound,P.RA=P)
  if (print_tables) print(sum_tb_ab)
  
  #rename P variable for fusing with abundance table that also has P column, 
  # and compound variable to FracCont as all values are fractional
  #contributions. Also adds entry for unlabeled fraction
  print(paste0(compound,": Summarizing fractional contribution data"))
  sum_tb_FC<-summarize_compounddata(filter(compound_tb,datatype=="FracCont"),
                                    compound,fact.names)%>%
    add_unlabeled_sum(compound,fact.names = fact.names)%>%
    rename(FracCont=compound,P.FC=P)
  if (print_tables) print(sum_tb_FC)
  
  #gather abundance and fraccont data together in one input table with one 
  #entry per pie slice (per combination cohort and labeling origin) with all 
  #other required info including labels for plotting function
  print(paste0(compound,": combining abundance and fractional contribution",
               "data and adding info required for plotting pie slices and",
               "labels"))
  slice_tb<-prepare_slicedata(compound_tb,sum_tb_FC,fact.names = fact.names,
                              compound=compound,labelDecimals = labelDecimals,
                              minLabDist = minLabDist,percentAdd = percentAdd,
                              FCposition = FCposition)
  
  if (print_tables) print(slice_tb)
  
  savepath<-path
  #plot based on information in slice table
  print(paste0(compound,": Building chart"))
  
  #plot detailed pie chart based on information in slice table
  print(paste0(compound,"Building detailed pie chart"))
  
  #invert labeling for pies
  lablevels_inv<-
    levels(slice_tb$Labeling)[length(levels(slice_tb$Labeling)):1]
  slice_tb$Labeling<- factor(slice_tb$Labeling,levels=lablevels_inv)
  
  #create starting barplot. X= halved abundances required, adds gridlines that
  #will become reference circles.  
  plotrect<-slice_tb %>% ggplot(aes(x = Abund/2, y = Fraction, fill = Labeling, 
                                    width = Abund)) + 
    geom_vline(xintercept=c(0.25),colour=circlelinecolor,
               linetype=circlelinetypes[1])+ 
    geom_vline(xintercept=c(0.5),colour=circlelinecolor,
               linetype=circlelinetypes[2])+ 
    geom_vline(xintercept=c(0.75),colour=circlelinecolor,
               linetype=circlelinetypes[3])+ 
    geom_vline(xintercept=c(1),colour=circlelinecolor,
               linetype=circlelinetypes[4])+ 
    geom_bar(stat = "identity", position = "fill")
  
  
  #add name of compound if desired, and the assign colors and thier legend order
  # if (include_name) plotrect<-plotrect+ggtitle(compound)
  plotrect<-plotrect  +
    scale_fill_manual(values=colLabeling,guide=guide_legend(reverse=T))
  
  #transform bar to pie chart and plot pies on grid, depending on amount of 
  #factors.
  if (twofactor) {
    gridformula<-as.formula(paste0(fact.names[2],"~",fact.names[1]))
    #switch="both" to set labels to same side as axis titles
    piebasic<-plotrect+
      facet_grid(gridformula,switch="both") +   
      coord_polar("y", start = 0, direction = 1) 
  } else {
    piebasic<-plotrect+
      facet_wrap(vars(!!rlang::sym(fact.names[1])),ncol=maxcol_facet) +   
      coord_polar("y", start = 0, direction = 1)
  }
  
  #apply final formatting to pie plots. Removes x and y labels entirely, 
  #including the space reserved for them on the plot
  #sets relative abundance p values in upper right corner of pie plots
  pies<-piebasic +
    labs(x=NULL, y=NULL)+                           
    #Change plots to black on white, remove text axes (fraction) that interfere
    #with circles, axis ticks, fraction grid lines. set legend title size,
    #remove rectangles and background around factor levels, set factor levels
    #to right text size
    theme_bw(base_size = mapotherfontsize) +                      
    theme(axis.text = element_blank(),              
          axis.ticks = element_blank(),             
          panel.grid = element_blank(),            
          # legend.title = element_text(size = legendtitlesize),
          strip.background = element_rect(fill = NA, colour = NA), 
          strip.text = element_text(size = mapcohortsize))
  #removes legend if desired
  chart <-pies + theme(legend.position = "none")  
  
  mapwidth<-6.15                          
  mapheight<-2.81  
  
  #save detailed chart if required or print to rstudio plot
  print(paste0(compound,": Plotting or saving pie chart"))
  if (save_chart) {
    #set folder path to save pie charts if saving requested
    plotfilefolder<-paste0(savepath,"/Pie charts pathway/")
    plotfilepath<-paste0(plotfilefolder,plotfilename)
    if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder))
    ggsave(plotfilepath,chart,width=mapwidth,height=mapheight,units = "cm",
           device = "png")
  } else {
    print(chart)
  }
  
  
  print(c(messages,"Finished"))
}
#prep table outside function-------------------------------------------------------------------
# convert raw path to R usable path, input files and check input
path<-gsub("\\\\", "/", rawpath)
meta_tb<-read_metacsv_clean(file=paste(path,metadatafile,sep = "/"),
                            normalize = normalize,twofactor = twofactor,
                            fact.invert = fact.invert,
                            factX.levels = factX.levels,
                            factY.levels = factY.levels)
abund_tb<-read_csv_clean(paste(path,abundancefile,sep = "/"),remove_empty = T)
fraccon_tb<-read_csv_clean(paste(path,fracconfile,sep = "/"),remove_empty = T)

check_input(meta_tb,abund_tb,fraccon_tb,FCposition = FCposition,
            colLabeling = colLabeling)


#Get factor and compound names from input
fact.names<-get_factornames(tb=meta_tb,twofactor = twofactor)
if (length(compounds)==0) {
  compounds<-colnames(abund_tb)[2:ncol(abund_tb)]
  
}

#merge all input into one table
tb<-merge_input(meta_tb = meta_tb,abund_tb = abund_tb,fraccon_tb = fraccon_tb,
                compounds=compounds, normalize = normalize)
compounds<-colnames(tb)[(ncol(meta_tb)+1):(ncol(tb)-1)]
# compounds<-c("Hexose")

#initialize messages variable, convert raw input filepath to R path
messages<-NULL

#check input
if (any(!compounds %in% colnames(tb))) {
  stop(paste0("Some requested metabolite names in 'compounds' are not ",
              "among the abundance column names. Make sure all requested",
              " compounds appear with the same name in the input files"))
}


#loop over each compound in input tibble
for (compound in compounds) {
  print(paste0("Processing ",compound))
  
  #prepare filename if saving required
  if (save_chart) {
    if (normalize) {
      plotfilename<-paste0("pies normalized ",compound,".png")
    } else {
      plotfilename<-paste0("pies ",compound,".png")
    }
  }
  
  #get table with only measured compound data, then a table summarizing
  #derived means and p values per cohort for abundance and one for fractional
  #contribution, then put together table with inputformat for pie function
  print(paste0(compound,": Extracting compound data"))
  compound_tb<-obtain_compounddata(tb,compound,fact.names)
  if (print_tables) print(compound_tb)
  
  #rename P variable for fusing with FC table that also has P column, and 
  #compound variable to Abund as all values are abundances
  print(paste0(compound,": Summarizing abundance data"))
  sum_tb_ab<-summarize_compounddata(filter(compound_tb,datatype=="Abund"),
                                    compound,fact.names)%>%
    rename(Abund=compound,P.RA=P)
  if (print_tables) print(sum_tb_ab)
  
  #rename P variable for fusing with abundance table that also has P column, 
  # and compound variable to FracCont as all values are fractional
  #contributions. Also adds entry for unlabeled fraction
  print(paste0(compound,": Summarizing fractional contribution data"))
  sum_tb_FC<-summarize_compounddata(filter(compound_tb,datatype=="FracCont"),
                                    compound,fact.names)%>%
    add_unlabeled_sum(compound,fact.names = fact.names)%>%
    rename(FracCont=compound,P.FC=P)
  if (print_tables) print(sum_tb_FC)
  
  #gather abundance and fraccont data together in one input table with one 
  #entry per pie slice (per combination cohort and labeling origin) with all 
  #other required info including labels for plotting function
  print(paste0(compound,": combining abundance and fractional contribution",
               "data and adding info required for plotting pie slices and",
               "labels"))
  
  slice_tb<-prepare_slicedata(compound_tb,sum_tb_FC,fact.names = fact.names,
                              compound=compound,labelDecimals = labelDecimals,
                              minLabDist = minLabDist,percentAdd = percentAdd,
                              FCposition = FCposition)
  
  if (print_tables) print(slice_tb)
}

# End ---------------------------------------------------------------------
cbind(pull(frac_tb[,c(1)]),frac_tb[,c(4)])
refvalues<-c(1.01,1.01,1.01,1.01)
tgtvalues<-c(1.01,1.01,0,0)
kruskal.test(c(refvalues,tgtvalues),
             c(rep("Reference",length(refvalues)),
               rep("Target",length(tgtvalues))))$p.value

# Donotuse old Functions and libraries ---------------------------------------------------------------
#load font library. For windows only it loads these fonts for bitmap output
# as well, not required for other operating systems
library(extrafont) 
if (Sys.info()[['sysname']]=="Windows") loadfonts(device="win") 
#for reading input data as tibbles fully compatible with dplyr and ggplot2 functions
library(readr)        
#for manipulating data as ggplot2 compatible tibbles
library(dplyr)        
#for intuitive conversion of tibble columns to factors by "as_factor"
library(forcats)      
#for restructuring data tibbles to allow different calculations
library(tidyr)         
#for generating the pie chart plots
library(ggplot2)    
#to deal with overlaying labels
library(ggrepel)
#for overlaying pie chart plots on a metabolic map
library(magick)       


# function for checking if any column cell contains non-NA data
has_data <- function(x) { sum(!is.na(x)) > 0 } 

# function for checking if any column cell is different from 0
has_nonzero <- function(x) { any(x != 0)}         

# function for loading and cleaning abundance and FC files
read_csv_clean<- function(file,remove_empty=FALSE){
  input_tb<-read_csv(file = file,show_col_types = FALSE)
  if (remove_empty) {
    input_tb<-select_if(input_tb,has_data)          #drop empty columns
  }
  
  if (! "Sample" %in% colnames(input_tb))  {
    stop(paste0("No column names 'Sample' found in input, please put sample ",
                "names in a column named 'Sample'"))
  }
  
  return(input_tb)
}
# function for loading and cleaning metadata files
read_metacsv_clean<- function(file,normalize=T,twofactor=F,fact.invert=F,
                              factX.levels=NULL,factY.levels=NULL){
  #initialize messages variable
  messages<-NULL
  
  #load and clean metadatafile
  input_tb<-read_csv_clean(file = file,remove_empty = T)
  
  #check if normalisation column present with name "Normalisation"
  #rename if different capitulisation
  #if normalisation requested but no column give warning
  #add placeholder normalisation column with only 1's if not present
  if ("normalisation" %in% tolower(colnames(input_tb))) {
    colnames(input_tb)[which(tolower(colnames(input_tb))=="normalisation")]<-
      "Normalisation"
  } else {
    if (normalize==T) {
      messages<-c(messages,paste0("Normalisation requested but not applied as",
                                  " normalisation column empty or not provided",
                                  " in metadata"))
    }
    input_tb$Normalisation<-1
  }
  
  #check if tracer column present with name "Tracer"
  #rename if different capitulisation
  #add placeholder tracer column with only name "Labeled" if not present
  #are in the input file, as center labeling is not possible
  if ("tracer" %in% tolower(colnames(input_tb))) {
    colnames(input_tb)[which(tolower(colnames(input_tb))=="tracer")]<-"Tracer"
  } else {
    input_tb$Tracer<-"Labeled"
  }
  
  input_tb<-select(input_tb,Sample,Tracer,Normalisation,everything())
  
  #checks regarding cohort factors given
  #first check if secondary factor requested while twofactor is not enabled
  #then if at least two factors are present if twofactor requested
  #then make placeholder factor with 1 level if no factor provided
  if (twofactor==F & length(factY.levels)>0) {
    messages<-c(messages,paste0("One factor analysis selected, but levels for",
                                " 2nd factor specified. Proceeding with ",
                                "one factor"))
  }
  if (twofactor==T & ncol(input_tb)<5) {
    messages<-c(messages,paste0("Two factor analysis selected, but less than",
                                " two cohortfactors columns in metadata. ",
                                "Proceeding with one factor analysis."))
    twofactor<-F
  }
  if (twofactor==F & ncol(input_tb)<4) {
    messages<-c(messages,paste0("No factor column is present, making dummy ",
                                "factor with one level, will result in only",
                                "one cohort based on all samples in figures"))
    input_tb$Cohort<-"Sample mean"
  }
  
  #Give message with factors used
  if (twofactor) {
    messages<-c(messages,paste0("Factors used for two-factor analysis: ",
                                paste0(colnames(input_tb)[c(4,5)],
                                       collapse = ", ")))
    if (ncol(input_tb)>5) {
      messages<-c(messages,paste0("Unused metadata columns: ",
                                  paste0(colnames(input_tb)[6:ncol(input_tb)],
                                         collapse = ", ")))
    }
  } else {
    messages<-c(messages,paste0("Factors used for one-factor analysis: ",
                                colnames(input_tb)[4]))
    if (ncol(input_tb)>4) {
      messages<-c(messages,paste0("Unused metadata columns: ",
                                  paste0(colnames(input_tb)[5:ncol(input_tb)],
                                         collapse = ", ")))
    }
  }
  
  #make sure metadata columns except normalisation are of right type
  if (twofactor) {
    input_tb<-mutate(input_tb,across(c(1,2),as.character))
    input_tb<-mutate(input_tb,across(c(4,5),as.factor))
  } else {
    input_tb<-mutate(input_tb,across(c(1,2),as.character))
    input_tb<-mutate(input_tb,across(c(4),as.factor))
  }
  
  if (twofactor) {
    #invert factor order if requested
    if (fact.invert) input_tb<-input_tb %>% relocate(5,4,.after=3)                       
    
    #set factor level order, if none provided take order in 
    #input file
    fact.names<-colnames(input_tb)[c(4,5)]
    if (length(factX.levels)>0) {                                     
      input_tb[,fact.names[1]]<-fct_relevel(pull(input_tb[,fact.names[1]]),
                                            factX.levels)
    } else {
      factorder_input<-unique(as.character(pull(input_tb[,fact.names[1]])))
      input_tb[,fact.names[1]]<-fct_relevel(pull(input_tb[,fact.names[1]]),
                                            factorder_input)
    }
    if (length(factY.levels)>0) {                                     
      input_tb[,fact.names[2]]<-fct_relevel(pull(input_tb[,fact.names[2]]),
                                            factY.levels)
    } else {
      factorder_input<-unique(as.character(pull(input_tb[,fact.names[2]])))
      input_tb[,fact.names[2]]<-fct_relevel(pull(input_tb[,fact.names[2]]),
                                            factorder_input)
    }
  } else {
    if (fact.invert) {
      messages<-c(messages,(paste0("Factor inversion requested is only meaningful for ",
                                   "twofactor analysis. Ignored since performing one-factor",
                                   "analysis")))
    }
    fact.names<-colnames(input_tb)[c(4)]
    if (length(factX.levels)>0) {
      input_tb[,fact.names[1]]<-fct_relevel(pull(input_tb[,fact.names[1]]),
                                            factX.levels)
    } else {
      factorder_input<-unique(as.character(pull(input_tb[,fact.names[1]])))
      input_tb[,fact.names[1]]<-fct_relevel(pull(input_tb[,fact.names[1]]),
                                            factorder_input)
    }
  }
  print(messages)
  return(input_tb)
}

#get factor names from table
get_factornames<-function(tb,twofactor=F) {
  if (twofactor) {
    fact.names<-colnames(meta_tb)[c(4,5)]
  } else {
    fact.names<-colnames(meta_tb)[c(4)]
  }
}

#checks inputdata  with requested analysis parameters and generates warnings 
#when incompatible. Uses <<- to set variables outside function environment
check_input<-function(meta_tb,abund_tb,fraccon_tb,FCposition,colLabeling){
  #initialize messages variable
  messages<-NULL
  
  #make sure FCposition is set to slice when multiple tracer nutrients
  if (length(unique(meta_tb$Tracer))>1 & FCposition =="center"){
    FCposition <<- "slice"
    messages<-c(messages,paste0("FC label was requested to be in center, but",
                                " as multiple tracer nutrients were used, ",
                                "several labels will exist per pie. Putting",
                                " label in slice instead."))
  }
  
  #checks if the right amount of colors is set, sets right amount of default 
  #distinctive colors (amount of tracers +1 for unlabeled fraction) if not
  if (!length(unique(meta_tb$Tracer)) == length(colLabeling)-1){
    if(length(unique(meta_tb$Tracer))==1) {
      colLabeling<<-c("#bfbfbf","#ffd966")
    } else {
      library(RColorBrewer)
      colLabeling<<-brewer.pal(length(unique(meta_tb$Tracer))+1,"Accent")
    }
    messages<-c(messages,paste0("Amount of label colors (colors for pie",
                                "slices) is not equal",
                                " to the amount of tracer + 1 for unlabeled ",
                                "fraction. Default colours used instead"))
  }
  
  if (any(!c(abund_tb$Sample,fraccon_tb$Sample) %in% meta_tb$Sample )) {
    messages<-c(messages,paste0("There are samples in the abundance and/or ",
                                "fractional contribution file that are not in ",
                                "the metadata file. Only samples in the ",
                                "metadata file will be taken into  account,",
                                "make sure that other samples can be safely ",
                                "ignored."))
  }
  
  if (any(!meta_tb$Sample %in% abund_tb$Sample)) {
    stop(paste0("Not all samples requested in the metadata file are present",
                "in the abundance file. Make sure all requested",
                " samples appear with the same name in the sample column of",
                " the abundance file"))
  }
  if (any(!meta_tb$Sample %in% fraccon_tb$Sample)) {
    stop(paste0("Not all samples requested in the metadata file are present",
                "in the fractional contribution file. Make sure all requested",
                " samples appear with the same name in the sample column of",
                " the fractional contribution file"))
  }
  
  #print messages
  print(messages)
}


# function to merge different input files into a tibble with all info needed
# to generate pies for all compounds
merge_input<-function(meta_tb,abund_tb,fraccon_tb,compounds,
                      normalize=T) {
  #initialize messages variable
  messages<-NULL
  
  #check input
  if (any(!compounds %in% colnames(abund_tb)[2:ncol(abund_tb)])) {
    stop(paste0("Some requested metabolite names in 'compounds' are not ",
                "among the abundance column names. Make sure all requested",
                " compounds appear with the same name in the input files"))
  }
  #check if compounds in fractional contribution table are absent from abundance
  #table
  if (any(!colnames(fraccon_tb) %in% colnames(abund_tb))) {
    messages<-c(messages,paste0("Following metabolites in the fractional ",
                                "contribution file were not in the ",
                                "abundance file. They are ",
                                "dropped from the analysis:"))
    messages<-c(messages,
                paste0(colnames(fraccon_tb)[(which(!colnames(fraccon_tb) %in% 
                                                     colnames(abund_tb)))],
                       collapse = ", "))
    
    #remove metabolites without missing in abundance data
    #from fractional contribution data
    fraccon_tb<-fraccon_tb[,which(colnames(fraccon_tb) %in% colnames(abund_tb))]
  }
  
  #add metadata to abundance and fractional contribution data respectively
  #retaining only selected samples, and drop metabolites with 0 abundance
  #in every sample to avoid errors
  abund_tb<-left_join(meta_tb,abund_tb,by="Sample") %>%
    select(1:ncol(meta_tb),any_of(compounds)) %>%
    select_if(has_nonzero)
  
  fraccon_tb<-left_join(meta_tb,fraccon_tb,by="Sample") %>%
    select(1:ncol(meta_tb),any_of(compounds))
  
  if (any(!colnames(fraccon_tb) %in% colnames(abund_tb))) {
    messages<-c(messages,paste0("Following metabolites in the fractional ",
                                "contribution file had 0 abundance in every ",
                                "selected sample. They are ",
                                "dropped from the analysis:"))
    messages<-c(messages,
                paste0(colnames(fraccon_tb)[(which(!colnames(fraccon_tb) %in% 
                                                     colnames(abund_tb)))],
                       collapse = ", "))
    
    #remove metabolites without or with 0 abundance (filtered out on input)
    #from fractional contribution data
    fraccon_tb<-fraccon_tb[,which(colnames(fraccon_tb) %in% colnames(abund_tb))]
  }
  
  #add fractional contribution equal to 100% unlabeled to compounds in abundance
  #but not fraction labeling table
  if (any(!colnames(abund_tb) %in% colnames(fraccon_tb))) {
    messages<-c(messages,paste0("Following metabolites in the abundance ",
                                "file were not in the ", 
                                "fractional contribution file.",
                                " Their fractional contribution is ",
                                "considered to be 100% unlabeled in all",
                                "samples: "))
    
    
    nolabnames<-colnames(abund_tb)[which(! colnames(abund_tb) %in%
                                           colnames(fraccon_tb))]
    for (i in nolabnames) {
      fraccon_tb$new<-0
      colnames(fraccon_tb)[ncol(fraccon_tb)]<-i
    }
    messages<-c(messages,nolabnames)
  }
  
  #normalize abundances if requested and possible,
  #print message noting whether normalisation was applied
  if (normalize) {
    abund_tb[,(ncol(meta_tb)+1):ncol(abund_tb)] <- 
      abund_tb[,(ncol(meta_tb)+1):ncol(abund_tb)]  / meta_tb$Normalisation
    messages<-c(messages,"Normalisation applied")
  } else {
    messages<-c(messages,"No normalisation applied")
    
  }
  
  #Per compound adapt FC's below 0 (artefacts due to natural abundance 
  #correction) to be positive to avoid problems with the visualisations
  #later on
  for (i in (ncol(meta_tb)+1):ncol(fraccon_tb)) {
    if (any(fraccon_tb[,i]<0)) {
      FCs<-pull(fraccon_tb[,i])
      FCs[which(FCs<0)]<-FCs[which(FCs<0)]-min(FCs[which(FCs<0)])     
      fraccon_tb[,i]<-FCs
    }
  }
  
  #prepare tables for joining and join
  fraccon_tb$datatype<-"FracCont"
  abund_tb$datatype<-"Abund"
  tb<-full_join(fraccon_tb,abund_tb,by=colnames(abund_tb))  #join separate tb's
  
  print(messages)
  
  return(tb)
}


#Extract data for one compound in merged input
#add explicitly the unlabeled fraction
#need to use !! for dynamic variable names in tidyverse selection
#see https://stackoverflow.com/questions/50537164/summarizing-by-dynamic-column-name-in-dplyr 
obtain_compounddata<-function(tb,compound,fact.names){
  compound_tb<-tb %>% select(Tracer,!!fact.names,datatype,
                             !!compound)
}

summarize_addP<-function(tb,cohortcolumn,valuecolumn,
                         Dtype=c("checkColumn","Abundance","FracCont")){
  #if datatype is provided in column, sort tb per datatype to make sure order is
  # ok for rest of function. Otherwise, check if datatype provided as variable,
  # and add column with only that type.
  #If so, set to that datatype, if not, error.
  if (Dtype=="checkColumn"){
    if ("datatype" %in% colnames(tb)) {
      tb<-tb[order(tb$datatype),]
    } else {
      print(paste0("summarize_P function requested to check for datatype in ",
                   "tibble column called 'datatype' (default option), but no such column ",
                   "provided. Either provide column name or specify datatype in function",
                   "call"))
    }
  } else {
    tb$datatype<-Dtype
  }
  
  
  #initialize tibble for output with one entry per factor level each for abund
  #and fraccont, with initialized column for p values, and an index noting
  #the last row in the P column that received data
  tb_out<-unique(tb[,-which(colnames(tb)==valuecolumn)])
  tb_out$P<-NA
  index<-0
  
  #loop over datatypes supplied
  for (j in unique(tb$datatype)){
    #create a separate tibble for each datatype to extract values
    datatype_selected<-j
    tb_type<-filter(tb,datatype==datatype_selected)
    
    #get cohorts names, extract first cohort as reference cohort, 
    #and obtain values of this cohort
    cohorts<-unique(pull(tb_type[,cohortcolumn]))
    refcohort<-cohorts[1]
    refvalues<-pull(tb_type[which(pull(tb_type[,cohortcolumn])==refcohort),
                            valuecolumn])
    
    #loop over target (non-reference) cohorts 
    for (i in 2:length(cohorts)) {  
      #extract values for current cohort
      tgtcohort<-cohorts[i]
      tgtvalues<-pull(tb_type[which(pull(tb_type[,cohortcolumn])==tgtcohort),
                              valuecolumn])
      
      
      #make P resultstring. If only one entry in cohort, show that no P could be
      #calculated by setting value to 99.
      #Otherwise perform appropriate test depending on datatype.
      #t.test for abundance data and kruskal wallis for fraccont
      #Set P=1 if all values are the same(likely 0) resulting in NaN. Make 
      #string depending on datatype
      if (length(tgtvalues)==1) {
        tb_out$P[index+i]<-99
      } else {
        if (datatype_selected=="Abund"){
          p<-t.test(refvalues,tgtvalues,)$p.value
          if (is.nan(p)) p<-1                 
          tb_out$P[index+i]<-p
        } else if (datatype_selected =="FracCont"){
          p<-kruskal.test(c(refvalues,tgtvalues),
                          c(rep("Reference",length(refvalues)),
                            rep("Target",length(tgtvalues))))$p.value             
          if (is.nan(p)) p<-1                 
          tb_out$P[index+i]<-p
        } else {
          stop(paste0("Datatype "),datatype_selected,
               " is not supported for P calculations P calculations only for Abund",
               " or FracCont")
        }
      }
    }
    #raise index by amount of cohorts in last set
    index<-index+i
  }
  
  return(tb_out)
}

#Make table with averages of datatype per cohort
#Calculates p values of significance tests of both relative abundance, and
#fractional contribution for each tracer, for printing on pie charts
#Group the table by cohort and calculate the mean per cohort and datatype
#then adds the P values calculated of each tracer per 
#combination of tracer and cohort factors
#drop grouping structure afterwards to avoid unexpected issues in the future
summarize_compounddata<-function(compound_tb,compound,fact.names){
  #factors and compounds need to be symbolized to use in 
  #tidyverse grouping function
  fact_symbols<-rlang::syms(fact.names) #list of symbols if multiple names
  comp_symbol<- rlang::sym(compound) #one symbol
  
  #get mean abundance and fractional contribution of each tracer per 
  #combination of tracer and cohort factors
  #need to use !! for dynamic variable names from one symbol and to 
  #use !!!  to symbolize list of symbols for group/summarise strings
  #see https://stackoverflow.com/questions/50537164/summarizing-by-dynamic-column-name-in-dplyr
  sum_tb<-group_by(compound_tb,Tracer,!!! fact_symbols,datatype)%>%
    summarise(!!compound := mean(!! comp_symbol),.groups = "drop")
  
  #Calculates p values of significance tests of both relative abundance, and
  #fractional contribution for each tracer per combination of tracer and cohort 
  #factors. Then joins to means and move P column to end
  tb_withP<-compound_tb %>% select(Tracer,!!fact.names,datatype,
                                   !!compound)
  if (length(fact.names)==2){
    tb_withP<-group_by(tb_withP,Tracer,!!!rlang::syms(fact.names[2]))
  } else if (length(fact.names)==1){
    tb_withP<-group_by(tb_withP,Tracer)
  }
  tb_withP<-group_modify(tb_withP,~summarize_addP(.x,cohortcolumn = fact.names[1],
                                                  valuecolumn = compound,Dtype = "checkColumn"))%>%
    ungroup()%>%
    right_join(sum_tb)%>%
    relocate(P, .after = last_col())
}

#add average unlabeled FC to summarized table with labeled FC's
add_unlabeled_sum<-function(sum_tb_FC,compound,fact.names){
  #tidyverse grouping function
  tracer_symbol<-rlang::syms(unique(sum_tb_FC$Tracer))
  
  #Calculate the unlabeled fraction for each sample. Then put back
  #in right format by joining to required info and entering missing info
  FC_tb<-sum_tb_FC%>%
    select(!P)%>%
    pivot_wider(names_from=c(Tracer),values_from=compound) %>%
    rowwise()%>%   #require to make sum function on next line work per row
    mutate(Unlabeled = 1-sum(!!!tracer_symbol)) %>%
    ungroup()%>%        #undo rowwise grouping
    pivot_longer(c(!!!tracer_symbol,Unlabeled),names_to = "Tracer",
                 values_to = compound)%>%
    left_join(select(sum_tb_FC,!c(compound,datatype)),
              by=c(fact.names,"Tracer"))%>%
    mutate(datatype=if_else(is.na(datatype),"FracCont",datatype))%>%
    relocate(P, .after = last_col())
  
  #set labeling as factor
  FC_tb$Tracer<-as.factor(FC_tb$Tracer)
  
  return(FC_tb)
}


#add fractional contribution labels and positions to pie table with requested
#formatting. 
add_FClabels<-function(slice_tb,labelDecimals,percentAdd,fact.names,FCposition,
                       minLabDist){
  slice_tb<-rowwise(slice_tb) %>%    #needed to apply some functions per row  
    #Get label, set to ND if not detected in any sample in group. Set label
    #next to empty if labeling is requested in center
    mutate(FracCont=round(FracCont,labelDecimals+2),
           labFC=if_else(percentAdd,paste0(FracCont*100,"%"),
                         as.character(FracCont*100)),
           labFC=if_else(Abund==0,"ND",labFC),
           labFC=if_else(FCposition=="center"& Tracer=="Unlabeled","",
                         labFC))%>%
    group_by(!!!rlang::syms(fact.names)) %>%
    #get labeling positions on FC and abundance axes. Depends if in
    #center or in slice. Center if not detected (label ND)
    #If slice, set posFC as sum of current and all 
    #preceding FC's-half the current FC. Set posAb in slice at minLabDist radius  
    #if abundance smaller than twice minLabDist. 
    mutate(FClab_posAngle=if_else(FCposition=="center"|Abund==0,0,
                                  cumsum(FracCont)-FracCont/2),
           FClab_posDist=if_else(FCposition=="center"|Abund==0,0,
                                 if_else(Abund<minLabDist*2,minLabDist,
                                         Abund/2)))%>%
    ungroup()                   #undo grouping
}

#make table with summarized data in the right format for pie creation,
#per slice. The average abundance normalized to the largest average abundance 
#is the pie radius. The fractions of the above parameter multiplied with the 
#labeled and unlabeled fraction correspond to the desired slices of a pie with
#this radius 
prepare_slicedata<-function(compound_tb,sum_tb_FC,compound,fact.names,
                            labelDecimals,percentAdd,FCposition,minLabDist){
  
  #factors and compounds need to be symbolized to use in 
  #tidyverse grouping function
  fact_symbols<-rlang::syms(fact.names) #list of symbols if multiple names
  comp_symbol<- rlang::sym(compound) #one symbol
  
  #Get abundance per sample and drop tracer column as we want to sum 
  #disregarding tracer,and P if present as it will be recalculated
  sum_tb_ab<-filter(compound_tb,datatype=="Abund")%>%
    select(!!! fact_symbols, !!comp_symbol,Abund=compound)
  
  #get average abundance per combination of cohort factors
  #meaning averaging over tracers as abundance should not be not tracer
  #dependent.
  slice_ab_tb<-sum_tb_ab%>%
    group_by(!!! fact_symbols) %>%
    summarise(Abund := mean(Abund),.groups = "drop")
  
  #Get abundance p values taken together independent of tracer
  #after dropping existing P column, join.
  if (length(fact.names)==2){
    sum_tb_ab<-group_by(sum_tb_ab,!!!rlang::syms(fact.names[2]))
  }
  slice_ab_tb<-group_modify(sum_tb_ab,~summarize_addP(.x,
                                                      cohortcolumn = fact.names[1],
                                                      valuecolumn = "Abund",Dtype = "Abund"))%>%
    ungroup()%>%
    right_join(slice_ab_tb,by=fact.names)%>%
    select(P.RA=P,!c(datatype,P))%>%
    relocate(P.RA, .after = last_col())
  
  
  #join abundance table to all tracer to add results to unlabeled as well,
  #then join with all and calculate abundance normalized to biggest abundance and
  #fraction of abundance per type of label
  slice_tb<-full_join(slice_ab_tb,sum_tb_FC,by=fact.names)%>%
    mutate(Abund=Abund/max(Abund),Fraction=FracCont*Abund)
  
  #Add FC labels and their positions, make P label and add P label radius 
  #positions, set informative table names and clean up unneccesary columns
  slice_tb<-add_FClabels(slice_tb,labelDecimals=labelDecimals,
                         percentAdd=percentAdd,fact.names = fact.names,
                         FCposition=FCposition,minLabDist=minLabDist)%>%
    rowwise()%>%
    mutate(P.FClab=case_when(
      is.na(P.FC) ~ "",
      P.FC==99 ~ "N=1,P=NA",
      length(unique(slice_tb$Tracer))>2 & P.FC<0.05 ~ "*",
      length(unique(slice_tb$Tracer))>2 & P.FC<0.1 ~ "`",
      length(unique(slice_tb$Tracer))>2 & P.FC>=0.05 ~ "",
      length(unique(slice_tb$Tracer))<=2 & P.FC<0.05 ~paste0("pFC=",
                                                             round(P.FC,2),
                                                             "*"),
      length(unique(slice_tb$Tracer))<=2 & P.FC>=0.05 ~ paste0("pFC=",
                                                               round(P.FC,2))),
      P.RAlab=case_when(
        is.na(P.RA) ~ "",
        P.RA==99 ~ "N=1,P=NA",
        P.RA<0.05 ~ paste0("pRA=",round(P.RA,2),"*"),
        P.RA>=0.05 ~ paste0("pRA=",round(P.RA,2))),
      Labeling=Tracer)%>%
    select(!Tracer)%>%
    ungroup()
}



#makes pie chart based on table with required data per pie slice
make_piechart<-function(slice_tb,twofactor=twofactor,compound,
                        fact.names=fact.names,circlelinecolor="gray",
                        maxcol_facet=4,
                        circlelinetypes=c(1,1,1,1),yAxLab="",xAxLab="",
                        include_name=F,colLabeling=colLabeling,
                        otherfontsize=10,font="sans",legendtitlesize=10,
                        cohortsize=12,include_legend=T,invert_FCangle=T){
  
  #required actions when inverting FC and FC label angle if plotted opposite 
  # way automatically
  if (invert_FCangle) {
    lablevels_inv<-
      levels(slice_tb$Labeling)[length(levels(slice_tb$Labeling)):1]
    slice_tb$Labeling<- factor(slice_tb$Labeling,levels=lablevels_inv)
    RAlab_y=7/8
  } else {
    RAlab_y=1/8
  }
  
  #create starting barplot. X= halved abundances required, adds gridlines that
  #will become reference circles.  
  plotrect<-slice_tb %>% ggplot(aes(x = Abund/2, y = Fraction, fill = Labeling, 
                                    width = Abund)) + 
    geom_vline(xintercept=c(0.25),colour=circlelinecolor,
               linetype=circlelinetypes[1])+ 
    geom_vline(xintercept=c(0.5),colour=circlelinecolor,
               linetype=circlelinetypes[2])+ 
    geom_vline(xintercept=c(0.75),colour=circlelinecolor,
               linetype=circlelinetypes[3])+ 
    geom_vline(xintercept=c(1),colour=circlelinecolor,
               linetype=circlelinetypes[4])+ 
    geom_bar(stat = "identity", position = "fill") + #makes basic rectangle plot
    #due to the pie manipulations the X and Y axes get inverted
    #so labels are assigned inverted too
    labs(x=yAxLab,y=xAxLab) 
  
  #add name of compound if desired, and the assign colors and thier legend order
  if (include_name) plotrect<-plotrect+ggtitle(compound)
  plotrect<-plotrect  +
    scale_fill_manual(values=colLabeling,guide=guide_legend(reverse=T))
  
  #positions of text at specified locations. GGrepel used when multiple tracer
  # to avoid labels overlapping. Fontsize needs to be adjusted for reasons: 
  #https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  if (length(unique(slice_tb$Labeling))>2) {
    plotrect<-plotrect  +
      geom_text_repel(aes(label=paste0(labFC,P.FClab)),
                      x = slice_tb$FClab_posDist,y=slice_tb$FClab_posAngle,
                      size=otherfontsize*5/14, family=font,
                      point.size=NA,direction = "x")
  } else {
    plotrect<-plotrect  +
      geom_text(aes(label=labFC),x = slice_tb$FClab_posDist,
                y=slice_tb$FClab_posAngle,size=otherfontsize*5/14, family=font)+  
      geom_text(aes(label=P.FClab),x=1.6,y=5/8,size=otherfontsize*5/14,
                hjust="inward",vjust="inward",family=font)      
  }
  
  #transform bar to pie chart and plot pies on grid, depending on amount of 
  #factors.
  if (twofactor) {
    gridformula<-as.formula(paste0(fact.names[2],"~",fact.names[1]))
    #switch="both" to set labels to same side as axis titles
    piebasic<-plotrect+
      facet_grid(gridformula,switch="both") +   
      coord_polar("y", start = 0, direction = 1) 
  } else {
    piebasic<-plotrect+
      facet_wrap(vars(!!rlang::sym(fact.names[1])),ncol=maxcol_facet) +   
      coord_polar("y", start = 0, direction = 1)
  }
  
  #apply final formatting to pie plots. Removes x and y labels entirely, 
  #including the space reserved for them on the plot
  #sets relative abundance p values in upper right corner of pie plots
  pies<-piebasic +
    labs(x=NULL, y=NULL)+                           
    geom_text(aes(label=P.RAlab),x=1.6,y=RAlab_y,size=otherfontsize*5/14,
              hjust="inward",vjust="inward",family=font) + 
    #Change plots to black on white, remove text axes (fraction) that interfere
    #with circles, axis ticks, fraction grid lines. set legend title size,
    #remove rectangles and background around factor levels, set factor levels
    #to right text size
    theme_bw(base_size = otherfontsize) +                      
    theme(axis.text = element_blank(),              
          axis.ticks = element_blank(),             
          panel.grid = element_blank(),            
          plot.title = element_text(size = cohortsize, face = "bold"),
          legend.title = element_text(size = legendtitlesize),
          strip.background = element_rect(fill = NA, colour = NA), 
          strip.text = element_text(size = cohortsize))
  
  #removes legend if desired
  if (!include_legend) pies <-pies + theme(legend.position = "none")    
  
  return(pies)
}


# Generate pie chart plot for each compound and save if requested
generate_pies<-function(tb,compounds,pathway_charts,save_chart,savepath,normalize=T,
                        fact.names,labelDecimals,percentAdd,FCposition,
                        minLabDist,twofactor,circlelinecolor,circlelinetypes,
                        maxcol_facet=maxcol_facet,
                        yAxLab,xAxLab,include_name,colLabeling,otherfontsize,font,
                        legendtitlesize,cohortsize,include_legend,
                        invert_FCangle,
                        mapotherotherfontsize,mapcohortsize,width,height,
                        mapwidth,mapheight) {
  #initialize messages variable, convert raw input filepath to R path
  messages<-NULL
  
  #check input
  if (any(!compounds %in% colnames(tb))) {
    stop(paste0("Some requested metabolite names in 'compounds' are not ",
                "among the abundance column names. Make sure all requested",
                " compounds appear with the same name in the input files"))
  }
  
  
  #loop over each compound in input tibble
  for (compound in compounds) {
    print(paste0("Processing ",compound))
    
    #prepare filename if saving required
    if (save_chart) {
      if (normalize) {
        plotfilename<-paste0("pies normalized ",compound,".png")
      } else {
        plotfilename<-paste0("pies ",compound,".png")
      }
    }
    
    #get table with only measured compound data, then a table summarizing
    #derived means and p values per cohort for abundance and one for fractional
    #contribution, then put together table with inputformat for pie function
    print(paste0(compound,": Extracting compound data"))
    compound_tb<-obtain_compounddata(tb,compound,fact.names)
    if (print_tables) print(compound_tb)
    
    #rename P variable for fusing with FC table that also has P column, and 
    #compound variable to Abund as all values are abundances
    print(paste0(compound,": Summarizing abundance data"))
    sum_tb_ab<-summarize_compounddata(filter(compound_tb,datatype=="Abund"),
                                      compound,fact.names)%>%
      rename(Abund=compound,P.RA=P)
    if (print_tables) print(sum_tb_ab)
    
    #rename P variable for fusing with abundance table that also has P column, 
    # and compound variable to FracCont as all values are fractional
    #contributions. Also adds entry for unlabeled fraction
    print(paste0(compound,": Summarizing fractional contribution data"))
    sum_tb_FC<-summarize_compounddata(filter(compound_tb,datatype=="FracCont"),
                                      compound,fact.names)%>%
      add_unlabeled_sum(compound,fact.names = fact.names)%>%
      rename(FracCont=compound,P.FC=P)
    if (print_tables) print(sum_tb_FC)
    
    #gather abundance and fraccont data together in one input table with one 
    #entry per pie slice (per combination cohort and labeling origin) with all 
    #other required info including labels for plotting function
    print(paste0(compound,": combining abundance and fractional contribution",
                 "data and adding info required for plotting pie slices and",
                 "labels"))
    slice_tb<-prepare_slicedata(compound_tb,sum_tb_FC,fact.names = fact.names,
                                compound=compound,labelDecimals = labelDecimals,
                                minLabDist = minLabDist,percentAdd = percentAdd,
                                FCposition = FCposition)
    
    if (print_tables) print(slice_tb)
    
    #plot detailed pie chart based on information in slice table
    print(paste0(compound,"Building detailed pie chart"))
    pies<-make_piechart(slice_tb,twofactor = twofactor,fact.names = fact.names,
                        circlelinecolor = circlelinecolor,compound=compound,
                        circlelinetypes = circlelinetypes,
                        include_name = include_name,colLabeling = colLabeling,
                        font=font,otherfontsize = otherfontsize,
                        legendtitlesize =legendtitlesize,
                        cohortsize = cohortsize,include_legend = include_legend,
                        invert_FCangle = invert_FCangle)
    
    #save detailed chart if required or print to rstudio plot
    print(paste0(compound,": Plotting or saving pie chart"))
    if (save_chart) {
      #set folder path to save pie charts if saving requested
      plotfilefolder<-paste0(savepath,"/Pie charts/")
      plotfilepath<-paste0(plotfilefolder,plotfilename)
      if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder))
      ggsave(plotfilepath,pies,width=width,height=height,units = "cm",
             device = "png")
    } else {
      print(pies)
    }
    
    #plot summary pie chart for pathway based on information in slice table
    if (pathway_charts&!twofactor) {
      print(paste0(compound,"Building summary pie chart"))
      pies<-make_piechart(slice_tb,twofactor = twofactor,fact.names = fact.names,
                          circlelinecolor = circlelinecolor,compound=compound,
                          circlelinetypes = circlelinetypes,
                          maxcol_facet=maxcol_facet,
                          include_name = F,colLabeling = colLabeling,
                          font=font,otherfontsize = mapotherfontsize,
                          cohortsize = mapcohortsize,include_legend = F,
                          invert_FCangle = invert_FCangle)
      
      #save summary pie chart for metabolites if required or print to rstudio plot
      print(paste0(compound,": Plotting or saving pie chart for pathway"))
      if (save_chart) {
        plotfilefolder<-paste0(savepath,"/Pie charts pathway/")
        plotfilepath<-paste0(plotfilefolder,plotfilename)
        if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder))
        ggsave(plotfilepath,pies,width=mapwidth,height=mapheight,units = "cm",
               device = "png")
      } else {
        print(pies)
      }
    }
    
    
  }
  print(c(messages,"Finished"))
}






