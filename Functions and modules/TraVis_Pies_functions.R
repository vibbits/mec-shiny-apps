# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: Functions used by TraVis pies to convert raw input data to pie 
# charts. Called in modules of TraVis pies, but not inherently linked to R shiny
# functionality and could be used and useful outside shiny framework.


# Functions and libraries ---------------------------------------------------------------
#libraries for UI
library(dplyr)        #for faster.easier manipulation of data as tibbles
library(vroom)        #for easier file loading
library(forcats)      #for factor manipulation
library(readr)        #for writing .csv file of merged output
library(tidyr)        #for restructuring data tibbles
library(ggplot2)      #for generating the pie chart plots

#Need to install additional software for fonts using remotes package. 
# Checks if fonts are imported fonts so they can be used in TraVis Pies
# Imports fonts if not
check_install_fonts<-function(import_dir=NULL) {
  #Installs remotes if not yet installed 
  if (!require(remotes)) {
    install.packages("remotes")
    library(remotes)
  }
  
  #If rttf2pt not yet installed, uses remote to install the last known version of  
  #rttf2pt1 compatible with extrafont .
  if (!require(Rttf2pt1)) {
    remotes::install_version("Rttf2pt1", version = "1.3.8")
    library(Rttf2pt1)
  }
  
  #If rttf2pt installed but not the last known version of compatible with 
  #extrafont, uninstalls then uses remote to install the right version
  if (!packageVersion("Rttf2pt1")=="1.3.8") {
    detach("package:extrafont",unload = TRUE)
    detach("package:Rttf2pt1",unload = TRUE)
    remove.packages("Rttf2pt1")
    remotes::install_version("Rttf2pt1", version = "1.3.8")
    library(Rttf2pt1)
    library(extrafont)
  }
  
  #Installs extrafont if not yet installed 
  if (!require(extrafont)) {
    install.packages("extrafont")
    library(extrafont)
    
  } 
  
  if (length(fonts())>1) {
    print("Fonts already imported")
    return()
  }
  
  
  #load font library. For windows only it loads these fonts for bitmap output
  # as well, not required for other operating systems. Don't load import library 
  # unless never done before on this system or unless 
  print(paste0("No font import detected, importing fonts (can take a few ",
  "minutes). Should only run once ever on a platform"))
  if (length(import_dir)==0) {
    font_import(prompt=F)
  } else {
    font_import(import_dir,prompt=F)
  }
  loadfonts()
  
}


# function for checking if any column cell contains non-NA data
has_data <- function(x) { sum(!is.na(x)) > 0 } 

# function for checking if any column cell is different from 0
has_nonzero <- function(x) { any(x != 0)}         

# function for loading and cleaning abundance and FC files
read_csv_clean<- function(file,remove_empty=FALSE,perc_to_num=T){
  input_tb<-vroom::vroom(file = file, delim = ",",show_col_types = FALSE)
    
  #drop empty columns if desired
  if (remove_empty) {
    input_tb<-select_if(input_tb,has_data)          
  }
  
  #set percentage strings to fractions if desired
  if (perc_to_num){
    percolumns<-grep("%",input_tb)


    input_tb<-mutate(input_tb,across(all_of(percolumns),function(x) 
      as.numeric(sub(pattern="%", replacement = "",x,fixed = T))/100))
  }
  
  return(input_tb)
}

#function to prepare metadata to uniform format
format_metadata<-function(meta_tb,sample_column,factor_column,norm_column) {
  sample_symbol<-rlang::sym(sample_column)
  fact_symbol<-rlang::sym(factor_column)
  
  #check normalisation, drop Normalisation column if exists but not selected
  #add dummy if no normalisation required, otherwise rename correct column
  if ("Normalisation" %in% colnames(meta_tb) & norm_column != "Normalisation") {
    meta_tb<-select(meta_tb,-Normalisation)
  }
  if (norm_column == "None") {
    meta_tb$Normalisation<-1
  } else {
    meta_tb<-rename(meta_tb,Normalisation=all_of(norm_column))
  }
  
  #check factor column, add dummy if no factors given
  if (factor_column=="None" ) {
    meta_tb$Cohort<-"SingleCohort"
    factor_column<-"Cohort"
    fact_symbol<-rlang::sym(factor_column)
  }
  
  #Order columns, drop all unrequired columns and set type
  #factor set to single type if None, or use pull for as.factor 
  #best use := to use !! demasking environmental variable as name 
  #(as_factor might also work but I had issues and dropped it)
  #drop normalisation column if dummy
  meta_tb<-transmute(meta_tb,
                     !!sample_symbol := as.character(pull(meta_tb,
                                                          sample_column)),
                     !!fact_symbol := as.character(pull(meta_tb,factor_column)),
                     Normalisation=as.numeric(Normalisation))
  if (norm_column=="None") meta_tb<-select(meta_tb,-Normalisation)
  
  return(meta_tb)
}

# function to merge different input files into a tibble with all info needed
# to generate pies for all compounds
merge_input<-function(meta_tb,abund_tb,fraccon_tb,sample_col="Sample",
                      compounds) {
  
  #rename sample column in all inputs
  meta_tb<-rename(meta_tb,Sample=all_of(sample_col))
  abund_tb<-rename(abund_tb,Sample=all_of(sample_col))
  fraccon_tb<-rename(fraccon_tb,Sample=all_of(sample_col))
  
  #add metadata to abundance and fractional contribution data respectively
  #retaining only selected samples, and drop metabolites with 0 abundance
  #in every sample to avoid errors
  abund_tb<-left_join(meta_tb,abund_tb,by="Sample") %>%
    select(1:ncol(meta_tb),any_of(compounds)) %>%
    select_if(has_nonzero)
  
  fraccon_tb<-left_join(meta_tb,fraccon_tb,by="Sample") %>%
    select(1:ncol(meta_tb),any_of(colnames(abund_tb))) 
  
  #add fractional contribution equal to 100% unlabeled to compounds in abundance
  #but not fraction labeling table
  if (any(!colnames(abund_tb) %in% colnames(fraccon_tb))) {
    nolabnames<-colnames(abund_tb)[which(! colnames(abund_tb) %in%
                                           colnames(fraccon_tb))]
    for (i in nolabnames) {
      fraccon_tb$new<-0
      colnames(fraccon_tb)[ncol(fraccon_tb)]<-i
    }
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
  
  #prepare tables for joining and join, then order
  fraccon_tb$datatype<-"FracCont"
  abund_tb$datatype<-"Abund"
  tb<-full_join(fraccon_tb,abund_tb,by=colnames(abund_tb)) %>%
    select(colnames(meta_tb),datatype,everything())
  
  
  #add normalized abundances if normalization column provided
  #then remove normalisation factor
  if ("Normalisation" %in% colnames(meta_tb)) {
    normabund_tb<-tb %>% 
      filter(datatype=="Abund") %>%
      mutate(across((ncol(meta_tb)+2):ncol(tb),
                    function(x) x/Normalisation))
    normabund_tb$datatype<-"NormAbund"            
    
    tb<-full_join(tb,normabund_tb,by=colnames(tb)) %>%
      select(-Normalisation)
  }
  
  return(tb)
}

#Extract data for one compound in merged input, only for desired factor
#levels and set factor order
#need to use !! for dynamic variable names in tidyverse selection
#see https://stackoverflow.com/questions/50537164/summarizing-by-dynamic-column-name-in-dplyr 
obtain_compounddata<-function(tb,compound,fact_name,
                              fact_order=unique(pull(tb,fact_name)),
                              normalize=F){
  #prepare factor name symbol to use as target column name for mutate
  fact_symbol<-rlang::sym(fact_name)
  
  #select only one compound, filter to include normalized or non normalized
  #abundances, then select only given factor levels, then drops unused levels
  compound_tb<-tb %>% select(!!fact_name,datatype,!!compound) %>%
    filter(datatype %in% c("FracCont",
                           if_else(normalize,"NormAbund","Abund"))) %>%
    filter(!!fact_symbol %in% fact_order) %>%
    droplevels() %>%
    
    #Change factor variable from text into actual factor for visualisation and
    #significance testing. Set datatype to Abund if normalized 
    #abundances were used, then arrange data order to match the factor levels
    mutate(!!fact_symbol:=factor(!!fact_symbol,levels = fact_order),
           datatype=if_else(datatype=="NormAbund","Abund",datatype)) %>%
    arrange(!!fact_symbol)
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
  
  #if only 1 cohort is provided, set P to 1 for further checking
  if (!length(unique(pull(tb[,cohortcolumn])))>1) {
    tb_out$P<-99
    return(tb_out)
  }
  
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
      if (length(tgtvalues)==1|length(refvalues)==1) {
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

#add fractional contribution labels and positions to pie table with requested
#formatting. 
add_FClabels<-function(slice_tb,label_decimals,percent_add,fact_name,FC_position,
                       min_lab_dist){
  slice_tb<-rowwise(slice_tb) %>%    #to apply following functions per row  
    #Get label, set to ND if not detected in any sample in group. Set label
    #of unlabeled fraction to empty if labeling is requested in center
    mutate(FracCont=round(FracCont,label_decimals+2),
           labFC=if_else(percent_add,paste0(FracCont*100,"%"),
                         as.character(FracCont*100)),
           labFC=if_else(Abund==0,"ND",labFC),
           labFC=if_else(any(FC_position=="center"& Labeling=="Unlabeled",
                             FracCont==0 & Labeling=="Unlabeled"),"",labFC))%>%
    group_by(!!rlang::sym(fact_name)) %>%
    
    #get labeling positions on FC and abundance axes. Depends if in
    #center or in slice. Center if not detected (label ND)
    #If slice, set posFC as sum of current and all 
    #preceding FC's-half the current FC. Set posAb in slice at min_lab_dist radius 
    #if abundance smaller than twice min_lab_dist. 
    mutate(FClab_posAngle=if_else(FC_position=="center"|Abund==0,0,
                                  cumsum(FracCont)-FracCont/2),
           FClab_posDist=if_else(FC_position=="center"|Abund==0|FracCont==1,
                                 as.double(0),
                                 if_else(Abund<min_lab_dist*2,
                                         as.double(min_lab_dist),
                                         as.double(Abund/2))))%>%
    ungroup()                   #undo grouping
}

#make table with summarized data in the right format for pie creation,
#per slice. The average abundance normalized to the largest average abundance 
#is the pie radius. The fractions of the above parameter multiplied with the 
#labeled and unlabeled fraction correspond to the desired slices of a pie with
#this radius. P values for relative abundance and fractional contribution
#are calculated and a label for these on the pie chart is generated
prepare_slicedata<-function(compound_tb,compound,fact_name,
                            label_decimals,percent_add,FC_position,min_lab_dist){
  #factor and compound need to be symbolized to use in 
  #tidyverse grouping function
  fact_symbol<-rlang::sym(fact_name)
  comp_symbol<- rlang::sym(compound) 
  
  #Calculates p values of significance tests of both relative abundance, and
  #fractional contribution per cohort factor level.
  P_tb<-compound_tb %>% select(!!fact_name,datatype,
                               !!compound)%>%
    summarize_addP(cohortcolumn = fact_name,valuecolumn = compound,
                   Dtype = "checkColumn")%>%
    mutate(datatype=if_else(datatype=="Abund","P.RA","P.FC"))%>%
    pivot_wider(names_from=datatype,values_from=P) 
  
  #get mean abundance and fractional contribution per cohort factor level
  #then join with P data from above.
  #Then format as table with single entry per cohort with abundance and FC as
  #variables. Use to calculate abundance normalized to biggest abundance and
  #fraction of abundance labeled and unlabeled, and finally the fractional 
  #contribution of the unlabeled part. Format as table with two entries
  #factor level, one for the labeled part and one for the unlabeled part
  sum_tb<-group_by(compound_tb,!! fact_symbol,datatype)%>%
    summarise(!!compound := mean(!!comp_symbol),.groups = "drop")%>%
    left_join(P_tb)%>%
    pivot_wider(names_from=datatype,values_from=!!compound)%>% 
    mutate(Abund=Abund/max(Abund),Labeled=FracCont*Abund,
           Unlabeled=(1-FracCont)*Abund) %>%  
    pivot_longer(Labeled:Unlabeled,names_to="Labeling",values_to="Fraction") %>%
    mutate(FracCont=if_else(Labeling=="Unlabeled",1-FracCont,FracCont))
  
  
  #sets Labeling column factor order to unlabeled then labeled, makes 
  #make_piechart plotting function result more intuitive
  sum_tb$Labeling <- factor(sum_tb$Labeling,levels=c("Unlabeled","Labeled"))  
  
  #Add FC labels and their positions, make P label and add P label radius 
  #positions, set informative P label names
  slice_tb<-add_FClabels(sum_tb,label_decimals=label_decimals,
                         percent_add=percent_add,fact_name = fact_name,
                         FC_position=FC_position,min_lab_dist=min_lab_dist)%>%
    rowwise()%>%
    mutate(P.FClab=case_when(
        is.na(P.FC)       ~ "",
        P.FC==99          ~ "N=1,P=NA",
        P.FC<0.05         ~ paste0("pFC=",round(P.FC,2),"*"),
        P.FC>=0.05        ~ paste0("pFC=",round(P.FC,2))),
      P.RAlab=case_when(
        is.na(P.RA)       ~ "",
        P.RA==99          ~ "N=1,P=NA",
        P.RA<0.05         ~ paste0("pRA=",round(P.RA,2),"*"),
        P.RA>=0.05        ~ paste0("pRA=",round(P.RA,2))))%>%
    ungroup()
}


#makes pie chart based on table with required data per pie slice
make_piechart<-function(slice_tb,compound,
                        fact_name=fact_name,circlelinecolor="gray",
                        maxcol_facet=4,
                        circlelinetypes=c(1,1,1,1),
                        include_name=F,col_labeling=c("#bfbfbf","#ffd966"),
                        alpha=0.7,
                        otherfontsize=10,font="sans",legendtitlesize=10,
                        cohortsize=12,include_legend=T,show_P=T){
  
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
    geom_bar(stat = "identity", position = "fill",alpha=alpha) 
  
  #add name of compound if desired, and the assign colors and their legend order
  if (include_name) plotrect<-plotrect+ggtitle(compound)
  plotrect<-plotrect  +
    scale_fill_manual(values=col_labeling,guide=guide_legend(reverse=T))
  
  #positions of text at specified locations, if desired. 
  #Fontsize needs to be adjusted for reasons: 
  #https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  plotrect<-plotrect  +
    geom_text(aes(label=labFC),x = slice_tb$FClab_posDist,
              y=slice_tb$FClab_posAngle,size=otherfontsize*5/14)
  if (show_P) {
    plotrect<-plotrect  +  
      geom_text(aes(label=P.RAlab),x=1.6,y=7/8,size=otherfontsize*5/14,
                hjust="inward",vjust="inward") +
      geom_text(aes(label=P.FClab),x=1.6,y=5/8,size=otherfontsize*5/14,
                hjust="inward",vjust="inward") 
  }
  
  #transform bar to pie chart and plot pies on grid.
  piebasic<-plotrect+
    facet_wrap(vars(!!rlang::sym(fact_name)),ncol=maxcol_facet) +   
    coord_polar("y", start = 0, direction = 1)
  
  
  #apply final formatting to pie plots. Removes x and y labels entirely, 
  #including the space reserved for them on the plot
  #sets relative abundance p values in upper right corner of pie plots
  pies<-piebasic +
    labs(x=NULL, y=NULL)+                           
    #Change plots to black on white, remove text axes (fraction) that interfere
    #with circles, axis ticks, fraction grid lines. Set text font,
    #set legend title size,
    #remove rectangles and background around factor levels, set factor levels
    #to right text size
    theme_bw(base_size = otherfontsize) +                      
    theme(axis.text = element_blank(),              
          axis.ticks = element_blank(),             
          panel.grid = element_blank(),
          text=element_text(family = font),
          plot.title = element_text(size = cohortsize, face = "bold"),
          legend.title = element_text(size = legendtitlesize),
          strip.background = element_rect(fill = NA, colour = NA), 
          strip.text = element_text(size = cohortsize))
  
  #removes legend if desired
  if (!include_legend) pies <-pies + theme(legend.position = "none")    
  
  return(pies)
}

generate_pie<-function(tb,compound,detail_charts,pathway_charts,savepath,
                       normalize=T,fact_name,fact_order,label_decimals,percent_add,
                       FC_position,min_lab_dist,circlelinecolor,circlelinetypes,
                       maxcol_facet,include_name,col_labeling,
                       alpha,otherfontsize,
                       font,legendtitlesize,cohortsize,include_legend,
                       mapotherfontsize=16,mapcohortsize=18,format="png",
                       show_P=T) {
  
  print(paste0("Processing ",compound))
  
  #prepare filename
  if (normalize) {
    plotfilename<-paste0("pies normalized ",compound,".",format)
  } else {
    plotfilename<-paste0("pies ",compound,".",format)
  }
  
  #get table with only measured compound data, then a table summarizing
  #derived means and p values per cohort for abundance and one for fractional
  #contribution, then put together table with inputformat for pie function
  print(paste0("extracting compounddata"))
  
  compound_tb<-obtain_compounddata(tb,compound,fact_name,
                                   fact_order = fact_order,
                                   normalize = normalize)
  
  #make table with summarized data in the right format for pie creation
  #each entry containing the needed info for one slice of one of the pie 
  #charts.The average abundance normalized to the largest average abundance 
  #is the pie radius. The fractions of the above parameter multiplied with the
  #labeled and unlabeled fraction correspond to the desired slices of a pie 
  #with this radius 
  print(paste0("preparing slice data"))
  
  slice_tb<-prepare_slicedata(compound_tb,fact_name = fact_name,
                              compound=compound,label_decimals = label_decimals,
                              min_lab_dist = min_lab_dist,
                              percent_add = percent_add,
                              FC_position = FC_position)
  
  if (detail_charts) {
    #plot detailed chart based on information in slice table
    print(paste0("saving detailed chart"))
    
    pies<-make_piechart(slice_tb,fact_name = fact_name,
                        circlelinecolor = circlelinecolor,compound=compound,
                        circlelinetypes = circlelinetypes,
                        maxcol_facet = maxcol_facet,
                        include_name = include_name,col_labeling = col_labeling,
                        alpha=alpha,font=font,otherfontsize = otherfontsize,
                        legendtitlesize =legendtitlesize,
                        cohortsize = cohortsize,include_legend = include_legend,
                        show_P=show_P)
    
    #save detailed pie chart for pathway if required
    plotfilefolder<-paste0(savepath,"/Pie charts/")
    plotfilepath<-paste0(plotfilefolder,plotfilename)
    if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder))
    ggsave(plotfilepath,pies,width=24.6,height=16,units = "cm",
           device = format)
  }
  
  if (pathway_charts) {
    print(paste0("saving pathway chart"))
    
    #plot summary pie chart for pathway based on information in slice table
    pies<-make_piechart(slice_tb,fact_name = fact_name,
                        circlelinecolor = circlelinecolor,compound=compound,
                        circlelinetypes = circlelinetypes,include_name = F,
                        maxcol_facet = maxcol_facet,
                        col_labeling = col_labeling,font=font,alpha=alpha,
                        otherfontsize = mapotherfontsize,
                        cohortsize = mapcohortsize,include_legend = F,
                        show_P=show_P)

    #save summary pie chart for pathway if required
    plotfilefolder<-paste0(savepath,"/Pie charts pathway/")
    plotfilepath<-paste0(plotfilefolder,plotfilename)
    if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder))
    ggsave(plotfilepath,pies,width=24.6,height=16,units = "cm",
           device = format)
  }
  
}

# Generate pie chart plot for each compound and save if requested
generate_multiple_pies<-function(tb,compounds,detail_charts,pathway_charts,savepath,
                        normalize=T,fact_name,fact_order,label_decimals,percent_add,
                        FC_position,min_lab_dist,circlelinecolor,circlelinetypes,
                        maxcol_facet,include_name,col_labeling,
                        alpha,otherfontsize,
                        font,legendtitlesize,cohortsize,include_legend,
                        mapotherfontsize=16,mapcohortsize=18,format="png",
                        show_P=T) {
  
  #loop over each compound in input tibble
  for (compound in compounds) {
    print(paste0("Compound ",which(compounds==compound),
                 " of ",length(compounds),")"))
    generate_pie(tb=tb,compound=compound,detail_charts=detail_charts,
                 pathway_charts=pathway_charts,savepath=savepath,
                 normalize=normalize,fact_name=fact_name,fact_order=fact_order,
                 label_decimals=label_decimals,percent_add=percent_add,
                 FC_position=FC_position,min_lab_dist=min_lab_dist,
                 circlelinecolor=circlelinecolor,
                 circlelinetypes=circlelinetypes,maxcol_facet=maxcol_facet,
                 include_name=include_name,col_labeling=col_labeling,
                 alpha=alpha,otherfontsize=otherfontsize,
                 font=font,legendtitlesize=legendtitlesize,
                 cohortsize=cohortsize,include_legend=include_legend,
                 mapotherfontsize=mapotherfontsize,mapcohortsize=mapcohortsize,
                 format=format,show_P=show_P)
    
  }  
  print("Finished")
}

