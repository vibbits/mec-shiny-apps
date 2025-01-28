# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: Functions used by TraVis pies to convert raw input data to pie 
# charts. Called in modules of TraVis pies, but not inherently linked to R shiny
# functionality and could be used and useful outside shiny framework.

#Libraries ---------------------------------------------------------------
#libraries for UI
library(dplyr)        #for faster.easier manipulation of data
library(tibble)       #for manipulating tibbles
library(vroom)        #for easier file loading
library(forcats)      #for factor manipulation
library(readr)        #for writing .csv file of merged output
library(tidyr)        #for restructuring data tibbles
library(stringr)       #for padding leading zeros to isotopologue strings
library(ggplot2)      #for generating the pie chart plots
library(ggrepel)      #for avoiding overlapping FC labels in multitracer plots

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


# Functions ---------------------------------------------

#Select only desired columns and filter only supported datatypes.
#Extract data only for desired factor levels and set factor order
#need to use !! for dynamic variable names in tidyverse selection
#see https://stackoverflow.com/questions/50537164/summarizing-by-dynamic-column-name-in-dplyr 
prepare_piedata<-function(tb,sample_column="Sample",factor_columns,
                          tracer_column,
                          factor_order=NULL){
  #prepare factor name symbol to use as target column name for mutate
  #select only one compound, filter to include normalized or non normalized
  # abundances
  sample_symbol<-sym_or_null(sample_column,allownull = F)
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  tracer_symbol<-sym_or_null(tracer_column)
  
  #if factor order is supplied as a vector, assume it is only the order
  #of the first factor
  if(length(factor_order)>0 & !is.list(factor_order)){
    factor_order<-list(factor_order)
  }

  #select only desired and present columns, keep only supported datatypes
  compound_tb<-tb %>% select(!!sample_symbol,!!!factor_symbols,!!tracer_symbol,
                             compound,datatype,value) %>%
    filter(datatype %in% c("FracCont","NormAbund","Abund") |
             grepl("iso",tolower(datatype),fixed = T))
  
  if(length(factor_symbols)>0) {
    for(i in 1:length(factor_symbols)) {
      #if factor order is not given, just assume levels and order already
      #in data
      if(length(factor_order)<i){
        factor_order[[i]]<-unique(pull(tb,!!factor_symbols[[i]]))
      }
      
      # select only given factor levels, then drops unused levels
      compound_tb<-compound_tb %>%
        clean_order_factors(factor_columns = factor_columns[i],
                            factor_order[i]) %>%
        
        #Just in case factor is loaded as character,
        #Change factor variable from text into actual factor for visualisation
        #and significance testing, then arrange data order to match the factor
        #levels
        mutate(!!factor_symbols[[i]]:=factor(!!factor_symbols[[i]],
                                          levels = factor_order[[i]]))%>%
        arrange(!!factor_symbols[[i]])
    }
  }
  
  return(compound_tb)
}

#function for generating krusal results in tibble format per cohort
#compared to reference cohort which is the first in the factor
kruskal_piedata<-function(data,test_formula,factor_column,factor_order){
  #make symbol for selection
  factor_symbol<-rlang::sym(factor_column)
  
  test_formula<-as.formula(paste0("value ~ ",factor_column))
  ref_cohort<-factor_order[[1]][1]
  
  #prepare tibble, set NA as default value in case calculation is not possible
  kwtest_results<-tibble(!!factor_symbol:=factor_order[[1]][-1],
                         p.value=NA)
  for(i in 2:length(factor_order[[1]])){
    tgt_cohort<-factor_order[[1]][i]
    partdata<-data %>% filter(!!factor_symbol %in% c(ref_cohort,tgt_cohort))
    
    #tests to avoid errors and assign other values. If there is only one
    #factor level present, leave P at NA. If all FC are exactly the same (e.g.
    #all 100% or all 0%) set P to 1
    if(length(unique(na.omit(pull(partdata,factor_column))))<2) next
    if(length(unique(na.omit(pull(partdata,value))))==1){
      kwtest_results<-kwtest_results %>%
        mutate(p.value=if_else(!!factor_symbol==tgt_cohort,
                               1,
                               p.value))
      next
    }
    
    #do Kruskal wallis test to determine P value
    kwtest_results<-kwtest_results %>%
      mutate(p.value=if_else(!!factor_symbol==tgt_cohort,
                             kruskal.test(formula=test_formula,data=partdata)$p.value,
                             p.value))
  }
  
  return(kwtest_results)
}

#Make table with averages of datatype per cohort
#Calculates p values of significance tests of both relative abundance,
#fractional contribution and isotopologues for each tracer
#todo make sure it supports no factor and no tracer as well
summarise_piedata<-function(prepare_tb,abund_string="abun",factor_column,
                            comparative_factor_column=NULL,tracer_column,
                            factor_order)
  {
  #test
  if(length(factor_columns)>2 ) stop("More than two factors were supplied")
  
  #prepare factor and tracer symbols
  factor_columns <- c(factor_column,comparative_factor_column)
  sample_symbol<-sym_or_null(sample_column,allownull = F)
  factor_symbol<-sym_or_null(factor_column)
  compar_factor_symbol<-sym_or_null(comparative_factor_column)
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  tracer_symbol<-sym_or_null(tracer_column)
  
  #if factor order is supplied as a vector, assume it is only the order
  #of the first factor
  if(length(factor_order)>0 & !is.list(factor_order)){
    factor_order<-list(factor_order)
  }
  #if a factor is given but factor order is not given, just assume levels and 
  #order already in data
  if(length(factor_columns)>0 & length(factor_order)==0){
    factor_order[[1]]<-unique(pull(prepare_tb,!!factor_symbols[[1]]))
  }
  
  ab_sum_tb<-prepare_tb%>%
    filter(grepl("abun",tolower(datatype))) %>%
    {
      if (length(comparative_factor_column)>0) {
        group_by(.,compound,!!compar_factor_symbol,datatype)
      } else {
        group_by(.,compound,datatype)
      }
    } %>%
    summarise(data = list(pick(everything())),.groups="keep") %>%
    rowwise() %>%
    mutate(mod=list(lm(
      as.formula(paste0("value ~ ",factor_column)),
      data=data))) %>%
    reframe(tidy(mod))%>%
    select(-estimate,-std.error,-statistic) %>%
    filter(!term=="(Intercept)")%>%
    mutate(term=if_else(term=="(Intercept)",
                        factor_order[[1]][1],
                        gsub(factor_column,"",term)),
           term=factor(term,levels=factor_order[[1]]))%>%
    rename(!!factor_symbol:=term)%>%
    right_join(
      prepare_tb %>%
        filter(grepl("abun",tolower(datatype))) %>%
        group_by(.,compound,!!!factor_symbols,datatype)%>%
        summarise(average = mean(value))
    )%>%
    mutate(across(any_of(tracer_column),~ ""))
  
  
  #this code can serve for testing problems in the list mutate per row
  # test<- prepare_tb %>%
  #   filter(grepl("frac",tolower(datatype))|
  #            grepl("iso",tolower(datatype))) %>%
  #   {
  #     if (length(comparative_factor_column)>0) {
  #       group_by(.,compound,!!compar_factor_symbol,!!tracer_symbol,datatype)
  #     } else {
  #       group_by(.,compound,!!tracer_symbol,datatype)
  #     }
  #   } %>%
  #   summarise(data = list(pick(everything())),.groups="keep") 
  # 
  # for (i in 1:nrow(test)){
  #   print(paste0(test$compound[i]," ",test$datatype[i]))
  #   parttb<-test%>% ungroup() %>% slice(i) %>% rowwise()
  #   # debug(kruskal_piedata)
  #   
  #   mutate(parttb,mod=list(kruskal_piedata(
  #     data=data,
  #     test_formula= as.formula(paste0("value ~ ",factor_column)),
  #     factor_column = factor_column,
  #     factor_order = factor_order)))
  # }
  
  
  #filter entries for compounds that are not detected in the sample 
  #(abundance not higher than 0). Otherwise the calculated FC would be wrongly
  #lowered
  prepare_FC_tb<- prepare_tb %>%
  filter(datatype=="Abund")%>%
  filter(value>0)%>%
  select(-datatype,-value)%>%
  left_join(prepare_tb)
  
  sum_tb<- prepare_FC_tb %>%
    
    #only check labeling data and group for all compound-tracer-comparative  
    #(2nd) factor combinations and for isotopologues in addition per 
    #isotopologue, then make one row per group with the remaining data in a 
    #list in the last column
    filter(grepl("frac",tolower(datatype))|
             grepl("iso",tolower(datatype))) %>%
    {
      if (length(comparative_factor_column)>0) {
        group_by(.,compound,!!compar_factor_symbol,!!tracer_symbol,datatype)
      } else {
        group_by(.,compound,!!tracer_symbol,datatype)
      }
    } %>%
    summarise(data = list(pick(everything())),.groups="keep") %>%
    
    #per row group, calculate p value
    rowwise()%>%
    mutate(mod=list(kruskal_piedata(
      data=data,
      test_formula= as.formula(paste0("value ~ ",factor_column)),
      factor_column = factor_column,
      factor_order = factor_order)))%>%
    
    #put result list back into tibble format without lists
    reframe(mod)%>%
    mutate(p.value=if_else(is.nan(p.value),
                           NA,
                           p.value))%>%
    
    #add P to fractional contributions 
    right_join(
      prepare_FC_tb %>%
        filter(grepl("frac",tolower(datatype))|
                 grepl("iso",tolower(datatype))) %>%
        group_by(.,compound,!!!factor_symbols,!!tracer_symbol,datatype)%>%
        summarise(average = mean(value))
    )%>%
    
    #add one entry per missing compound-factor combination, with NA as value
    right_join(
      ab_sum_tb %>%
        select(-datatype,-p.value,-average)%>%
        unique()
    )%>%
    mutate( datatype=if_else(is.na(datatype),"FracCont",datatype),
    )%>%
    {
      if(length(tracer_column)>0) {
        mutate(.,
               !!tracer_symbol:=if_else(is.na(!!tracer_symbol),
                                        pull(prepare_tb,tracer_column)[1],
                                        !!tracer_symbol))
      } else .
    }%>%
    
    #add abundance rows, rename P value column for later
    full_join(ab_sum_tb)%>%
    rename(P=p.value)
  
  return(sum_tb)
}



#todo remove useless compound references everywhere as column is now "value"
# also comp_symbol
#add average unlabeled FC to summarized table with labeled FC's
add_UnlabFC<-function(sum_tb,factor_columns,tracer_column){
  #get symbols tracer; and each tracer nutrient used as a symbol vector
  tracer_symbol<-rlang::sym(tracer_column)
  nutrient_symbols<-sum_tb%>%
    filter(grepl("frac",tolower(datatype),fixed = T)) %>%
    pull(tracer_column)%>%
    unique()%>%
    rlang::syms() 
  
  #First make sure fractions sum to 100%, if not divide each fraction by sum of 
  # fractions. Calculate the unlabeled fraction for each sample. Then put back
  #in right format by joining to required info and entering missing info
  sum_tb %>%
    #add unlabeled fractions as fractional contribution unlabeled
    filter(grepl("frac",tolower(datatype),fixed = T))%>%
    select(!P)%>%
    pivot_wider(names_from=!!tracer_symbol,values_from=average,
                values_fill = 0) %>%
    rowwise()%>%   #require to make sum function on next line work per row
    mutate(across(c(!!!nutrient_symbols),
                  .fns = ~ if_else(sum(!!!nutrient_symbols)>1,
                                   .x/sum(!!!nutrient_symbols),.x)),
           Unlabeled = 1-sum(!!!nutrient_symbols)) %>%
    ungroup()%>%        #undo rowwise grouping
    pivot_longer(c(!!!nutrient_symbols,Unlabeled),names_to = tracer_column,
                 values_to = "average")%>%
    left_join(select(sum_tb %>%
                       filter(grepl("frac",tolower(datatype),fixed = T)),
                     !c(average,datatype)),
              by=c("compound",factor_columns,tracer_column))%>%
    full_join(sum_tb)
}

#add fractional contribution labels and positions to pie table with requested
#formatting. 
add_FClabels<-function(sum_tb,fraction_column,label_decimals,percent_add,factor_columns,
                       tracer_column,FC_position,min_lab_dist){
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  tracer_symbol<-sym_or_null(tracer_column)
  if(length(tracer_symbol)>0) {
    nutrient_symbols<-sum_tb%>%
      filter(grepl("frac",tolower(datatype),fixed = T)) %>%
      pull(tracer_column)%>%
      unique()%>%
      rlang::syms()
  } else nutrient_symbols<-NULL
  
  fraction_symbol <- sym_or_null(fraction_column)
  tracer_symbol<-sym_or_null(tracer_column)
  if(length(tracer_column)>0) {
    tracers<-sum_tb%>%
      filter(grepl("frac",tolower(datatype),fixed = T)) %>%
      pull(tracer_column)%>%
      unique()
  } else tracers<-NULL
  
  
  if(length(tracers)>2){
    FC_position="slice"
  }
  
  test<-rowwise(sum_tb) %>%    #to apply following functions per row  
    #Get label, set to ND if not detected in any sample in group. Set label
    #of unlabeled fraction to empty if labeling is requested in center
    mutate(!!fraction_symbol:=round(!!fraction_symbol,label_decimals+2),
           !!fraction_symbol:=if_else(is.na(!!fraction_symbol),
                                      0,
                                      !!fraction_symbol),
           labFC=if_else(FC_position=="slice" & !!fraction_symbol==0,
                         paste0("<",10^-label_decimals/2),
                         as.character(!!fraction_symbol*100)),
           labFC=if_else(percent_add,paste0(labFC,"%"),
                         labFC),
           labFC=if_else(FC_position=="center"&
                           !!tracer_symbol=="Unlabeled","",labFC),
           labFC=if_else(Abund==0 | is.na(Abund),"ND",labFC)
    )%>%
    group_by(!!!rlang::syms(factor_columns)) %>%
    
    #get labeling positions on FC and abundance axes. Depends if in
    #center or in slice. Center if not detected (label ND)
    #If slice, set posFC as sum of current and all 
    #preceding FC's-half the current FC. Set posAb in slice at min_lab_dist radius 
    #if abundance smaller than twice min_lab_dist. 
    mutate(FClab_posAngle=if_else(FC_position=="center"|Abund==0|is.na(Abund),
                                  0,
                                  -cumsum(!!fraction_symbol)+!!fraction_symbol/2),
  FClab_posDist=if_else(FC_position=="center"|Abund==0|is.na(Abund)|
                          !!fraction_symbol==1,
                                 as.double(0),
                                 if_else(Abund<min_lab_dist*2,
                                         as.double(min_lab_dist),
                                         as.double(Abund/2))))%>%
    ungroup()                   #undo grouping
}

#todo update from FC functions
#transform in a shape useful for the pie chart plotting function isotopologue
#contribution. Keep only abundance and isotopologue data, keep only normalized
#or unnormalized abundances 
make_iso_slices<-function(sum_tb,factor_columns,tracer_column){
  #get symbols factor and tracer; and each tracer nutrient used
  factor_symbols<-rlang::syms(factor_columns)
  
  # only keep fractions larger than 0
  sum_tb%>%
    filter(grepl("iso",tolower(datatype),fixed = T),
           average>0) %>%
    rename(IsoCont=average,P_FC=P,Fraction=!!tracer_column)%>%
    left_join(
      sum_tb%>%
        {
          if (normalize) {
            filter(.,grepl("norm",tolower(datatype),fixed = T) & 
                     grepl("abund",tolower(datatype),fixed = T))
          } else {
            filter(.,!grepl("norm",tolower(datatype),fixed = T) & 
                     grepl("abund",tolower(datatype),fixed = T))
          }
        } %>%
        select(compound,!!!factor_symbols,P,average)%>%
        rename(Abund=average,P_RA=P),
      by = join_by(compound, !!!factor_symbols)
    )%>%
    group_by(compound)%>%
    mutate(Abund=Abund/max(Abund),
           Fraction=IsoCont*Abund)%>%
    select(-P_FC,-IsoCont)%>%
    mutate(IsoCont=Fraction/Abund)%>%
    na.omit()%>%
    left_join(
      sum_tb%>%
        rename(P_FC=P) %>%
        select(P_FC,compound,datatype,!!!factor_symbols)
    )%>%
    add_FClabels(fraction_column="IsoCont",label_decimals=label_decimals,
                 percent_add=percent_add,factor_columns = factor_columns,
                 tracer_column = "datatype",
                 FC_position="slice",min_lab_dist=min_lab_dist)%>%
    rowwise()%>%
    mutate(
      P_FClab=case_when(
        is.na(P_FC)       ~ "",
        P_FC==99          ~ "N=1,P=NA",
        P_FC<0.05         ~ paste0("pFC=",round(P_FC,2),"*"),
        P_FC>=0.05        ~ paste0("pFC=",round(P_FC,2))),
      P_RAlab=case_when(
        is.na(P_RA)       ~ "",
        P_RA==99          ~ "N=1,P=NA",
        P_RA<0.05         ~ paste0("pRA=",round(P_RA,2),"*"),
        P_RA>=0.05        ~ paste0("pRA=",round(P_RA,2)))
    ) %>%
    ungroup()%>%
    mutate(Isotopologue=factor(datatype))%>%
    select(-datatype)
}


#transform summarized data  in a shape useful for the pie chart plotting function 
#fractional contribution. Keep only abundance and fractional data, keep only normalized
#or unnormalized abundances  
make_FC_slices<-function(sum_tb,factor_columns,tracer_column,
                         abund_norm_comparlevel=F){
  
  #limit dataframe to useful data to make sure that no compounds remain that 
  #are never above 0 (negative values often overcorrections)
  compounds_detected<-sum_tb%>%
    filter(grepl("abund",tolower(datatype),fixed = T))%>%
    group_by(compound,datatype)%>%
    summarise(average=sum(average))%>%
    filter(average>0)%>%
    pull(compound)%>%
    unique()

  
  #get symbols factor and tracer; and each tracer nutrient used
  #prepare factor and tracer symbols. Need a placeholder at the start if no
  #tracer information was provided
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  
  if(length(tracer_column)==0){
    tracer_column<-"Labeling"
    sum_mod_tb<-sum_tb %>%
      mutate(Labeling="Labeled")
  } else sum_mod_tb<-sum_tb
  
  tracer_symbol<-sym_or_null(tracer_column)
  nutrients<-sum_mod_tb%>%
    filter(grepl("frac",tolower(datatype),fixed = T)) %>%
    pull(tracer_column)%>%
    unique()
  nutrient_symbols<-nutrients%>%
    rlang::syms()

  
  #take fractional contribution data only, keep only slices above 0
  # join abundances or if desired
  #normalized abundances to it, calculate the fractions as abundance*fraccon
  #for each tracer, and assume a single unnamed tracer when no
  #tracer column is provided. Rejoin the fraccon p values and add labels for
  #for the plots
  slice_tb<-sum_mod_tb%>%
    filter(grepl("frac",tolower(datatype),fixed = T)) %>%
    rename(FracCont=average,P_FC=P)%>%
    
    #if total fraccont is over 100% already, set total of fraccont to 100% 
    #by normalizing to sum of all contributions
    #removing small subtraction and rounding errors
    group_by(compound,!!!factor_symbols)%>%
    mutate(FracCont=if_else(FracCont<0,0,FracCont),
           FracCont=FracCont/(max(1,sum(FracCont,na.rm = T)))
    )%>%
    
    #join desired abundance data
    left_join(
      sum_mod_tb%>%
        {
          if (normalize) {
            filter(.,grepl("norm",tolower(datatype),fixed = T) & 
                     grepl("abund",tolower(datatype),fixed = T))
          } else {
            filter(.,!grepl("norm",tolower(datatype),fixed = T) & 
                     grepl("abund",tolower(datatype),fixed = T))
          }
        } %>%
        select(compound,!!!factor_symbols,P,average)%>%
        rename(Abund=average,P_RA=P),
      by = join_by(compound, !!!factor_symbols)
    ) %>%
    
    #per compound, and if desired per comparative factor level, mutate abundance 
    #to be at least 0, and then normalize to the maximal abundance present
    #in the group, then calculate fraction (slice) of the 
    #abundance (pie) for each entry
    {
      if(length(factor_columns>1)&abund_norm_comparlevel){
        group_by(.,compound,!!factor_symbols[[2]])
      } else {
        group_by(.,compound)
      }
    }%>%
    mutate(Abund=if_else(Abund<0,0,Abund),
           Abund=Abund/max(Abund),
           Abund=if_else(is.nan(Abund),0,Abund),
           Fraction=FracCont*Abund)%>%
    select(-P_FC,-FracCont)%>%
    {
      if (length(tracer_symbol)>0) {
        pivot_wider(.,names_from = !!tracer_symbol,values_from = Fraction)%>%
          rowwise()%>%
          mutate(Unlabeled=Abund-sum(!!!nutrient_symbols,na.rm = T))%>%
          pivot_longer(c(!!!nutrient_symbols,Unlabeled),names_to=tracer_column,
                       values_to="Fraction")
      } else {
        rename(.,Labeled=Fraction) %>%
          rowwise()%>%
          mutate(Unlabeled=Abund-sum(!!!nutrient_symbols,na.rm = T))%>%
          pivot_longer(c(Labeled,Unlabeled),names_to=Labeling,
                       values_to="Fraction")
      }
    }%>%
    
    #Obtain back fractional contribution from data, only keep compounds that
    #are detected in at least one group, and for groups that are not detected
    #only keep one entry (the Unlabeled entry will always be present)
    mutate(FracCont=Fraction/Abund)%>%
    filter(compound %in% compounds_detected,
           Abund>0|!!tracer_symbol=="Unlabeled")%>%
    left_join(
      sum_mod_tb%>%
        rename(P_FC=P) %>%
        select(P_FC,compound,datatype,!!tracer_symbol,!!!factor_symbols)
    )%>%
    
    
    #add FC labels and modify P labels accordingly
    add_FClabels(fraction_column="FracCont",label_decimals=label_decimals,
                 percent_add=percent_add,factor_columns = factor_columns,
                 tracer_column = tracer_column,
                 FC_position=FC_position,min_lab_dist=min_lab_dist)%>%
    rowwise()%>%
    mutate(
      P_FClab=case_when(
        is.na(P_FC)       ~ "",
        P_FC==99          ~ "N=1,P=NA",
        P_FC<0.05         ~ paste0("pFC=",round(P_FC,2),"*"),
        P_FC>=0.05        ~ paste0("pFC=",round(P_FC,2))),
      P_RAlab=case_when(
        is.na(P_RA)       ~ "",
        P_RA==99          ~ "N=1,P=NA",
        P_RA<0.05         ~ paste0("pRA=",round(P_RA,2),"*"),
        P_RA>=0.05        ~ paste0("pRA=",round(P_RA,2)))
    ) %>%
    ungroup()%>%
    mutate(
      !!tracer_symbol:=factor(!!tracer_symbol,
                              levels=c(nutrients,"Unlabeled"))
    )

  return(slice_tb)
}

#makes pie chart based on table with required data per pie slice
make_piechart<-function(slice_tb,selected_compound,tracer_column=tracer_column,
                        factor_columns=factor_columns,log_abund=F,
                        circlelinecolor="gray",maxcol_facet=4,
                        circlelinetypes=c(1,1,1,1),
                        include_name=F,col_labeling=c("#bfbfbf","#ffd966"),
                        alpha=0.7,
                        otherfontsize=10,font="sans",legendtitlesize=10,
                        cohortsize=12,include_legend=T,show_P=T){
    
  if(length(tracer_column)==0){
    tracer_column<-"Labeling"
  }
  tracer_symbol<-sym_or_null(tracer_column)
  
  #turn tracer column into factor if it isn't one yet,
  #save original levels in order of appearance
  #extract data of selected compound only, keep only fractions above 0
  #keep only tracer colors linked to existing levels
  if(is.factor(pull(slice_tb,tracer_column))) {
    slice_mod_tb <- slice_tb
  } else {
    slice_mod_tb <- slice_tb %>%
      mutate(!!tracer_symbol:=factor(!!tracer_symbol,
                                     levels=unique(!!tracer_symbol)))
  }
  orig_levels<-levels(pull(slice_mod_tb,tracer_column))
    
  slice_mod_tb <- slice_mod_tb %>%
    filter(compound==selected_compound,
           Fraction>0|Abund==0)
  col_labeling<-col_labeling[which(
    orig_levels %in% as.character(unique(pull(slice_mod_tb,tracer_column))))]
  #create starting barplot. X= halved abundances required, take log if requested
  #Adds gridlines that will become reference circles at 0.25 0.5 0.75 and 1 on 
  #normal scale or 0.001 0.01 0.1 and 1 on log scale. 
  if (log_abund) {
    #set minimal value to include on log axis, changing not recommended
    #and calculate minimal position distance on new scale
    minvalue<-0.0001
    
    #width can only be symmetric, so modify abundance to the value on normal 
    #scale corresponding to average of log scale minimal limit 
    #and logscale abundance, and abundance width to the difference of the log 
    #scale abundance and log scale minimal limit
    slice_mod_tb <- slice_mod_tb %>% 
      rowwise() %>%
      mutate(modAbund=10^((log10(Abund)+log10(minvalue))/2),
             modAbund_width=-(log10(minvalue)-log10(Abund)),
             #reset label distance position to work on logscale, either to middle
             #or to intended distance depending on whether one was given.
             FClab_posDist=if_else(FClab_posDist>0,
                                   log10(modAbund),
                                   log10(minvalue)))
    
    plotrect<-slice_mod_tb %>% ggplot(aes(x = modAbund, y = Fraction,  
                                      fill = !!tracer_symbol, 
                                      width = modAbund_width)) + 
      scale_x_log10(limits= c(minvalue, 1)) +
      geom_vline(xintercept=c(0.001),colour=circlelinecolor,
                 linetype=circlelinetypes[1])+ 
      geom_vline(xintercept=c(0.01),colour=circlelinecolor,
                 linetype=circlelinetypes[2])+ 
      geom_vline(xintercept=c(0.1),colour=circlelinecolor,
                 linetype=circlelinetypes[3])+ 
      geom_vline(xintercept=c(1),colour=circlelinecolor,
                 linetype=circlelinetypes[4])+ 
      geom_bar(stat = "identity", position = "fill",alpha=alpha) 
    
  } else {
    plotrect<-slice_mod_tb %>% ggplot(aes(x = Abund/2, y = Fraction,
                                      fill = !!tracer_symbol,  
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
  }
  
  
  #add name of compound if desired, and the assign colors and their legend order
  if (include_name) plotrect<-plotrect+ggtitle(selected_compound)
  if(length(col_labeling)>0) {
    plotrect<-plotrect  +
      scale_fill_manual(values=col_labeling,guide=guide_legend(reverse=F))
  } else {
    plotrect<-plotrect  +
      scale_fill_discrete()
    # scale_fill_manual(values=col_labeling,guide=guide_legend(reverse=T))
  }
  
  
  #positions of text at specified locations. GGrepel used when multiple tracer
  # to avoid labels overlapping. Fontsize needs to be adjusted for reasons:
  #https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  # plotrect<-plotrect  +
  #   geom_text(aes(label=labFC),x = slice_mod_tb$FClab_posDist,
  #             y=slice_mod_tb$FClab_posAngle,size=otherfontsize*5/14)
  if (show_P) {
    plotrect<-plotrect  +
      geom_text(aes(label=P_RAlab),x=1.6,y=7/8,size=otherfontsize*5/14,
                hjust="inward",vjust="inward")
    
    if (length(unique(pull(slice_tb[,tracer_column])))>2) {
      #modify label to include P label if exists
      slice_mod_tb<-slice_mod_tb%>%
        rowwise()%>%
        mutate(full_label=if_else(nchar(P_FClab)>1,
                                  paste0(labFC,"\n",P_FClab),
                                  labFC)
        )%>%
        ungroup()
      plotrect<-plotrect  +
        geom_text_repel(data=slice_mod_tb,
                        aes(label=full_label),
                        x = slice_mod_tb$FClab_posDist,y=slice_mod_tb$FClab_posAngle,
                        size=otherfontsize*5/14, family=font,
                        point.size=NA,direction = "both",
                        arrow = arrow())

    } else {
      plotrect<-plotrect  +
        geom_text(aes(label=labFC),x = slice_mod_tb$FClab_posDist,
                  y=slice_mod_tb$FClab_posAngle,size=otherfontsize*5/14, family=font)+  
        geom_text(aes(label=P_FClab),x=1.6,y=5/8,size=otherfontsize*5/14,
                  hjust="inward",vjust="inward",family=font)      
    }
  } else {
    if (length(unique(pull(slice_mod_tb[,tracer_column])))>2) {
      plotrect<-plotrect  +
        geom_text_repel(data=slice_mod_tb,
                        aes(label=paste0(labFC)),
                        x = slice_mod_tb$FClab_posDist,y=slice_mod_tb$FClab_posAngle,
                        size=otherfontsize*5/14, family=font,
                        point.size=NA,direction = "x",
                        arrow = arrow())
    } else {
      plotrect<-plotrect  +
        geom_text(aes(label=labFC),x = slice_mod_tb$FClab_posDist,
                  y=slice_mod_tb$FClab_posAngle,size=otherfontsize*5/14, family=font)+  
        geom_text(aes(label=P_FClab),x=1.6,y=5/8,size=otherfontsize*5/14,
                  hjust="inward",vjust="inward",family=font)      
    }
  }
  
  #transform bar to pie chart and plot pies on grid, depending on amount of 
  #factors.
  if(length(factor_columns)<2) {
    piebasic<-plotrect+
      facet_wrap(vars(!!rlang::sym(factor_columns)),ncol=maxcol_facet) +
      coord_polar("y", start = 0, direction = -1)
  } else {
    gridformula<-as.formula(paste0(factor_columns[2],"~",factor_columns[1]))
    #switch="both" to set labels to same side as axis titles
    piebasic<-plotrect+
      facet_grid(gridformula,switch="both") +   
      coord_polar("y", start = 0, direction = -1)
  }
  
  
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

#prepare tibble for input to function for generating pie charts
make_slicetibble<-
  function(tb,normalize=T,factor_columns,tracer_column,factor_order,label_decimals,
           percent_add,FC_position,min_lab_dist,P_isotopologues,show_P=T,
           abund_norm_comparlevel=F) {
    #derive variables used later on
    factor_columns <- c(factor_column,comparative_factor_column)
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
      # filter(compound=="Gal-6-P")%>%
      
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
    
    #obtain isotopologue slice tb for plotting if isotopologues provided
    isos_calculated<-F
    if(!any(grepl("iso",tolower(sum_tb$datatype)))) {
      if(exists("isoslice_tb")) rm("isoslice_tb")
    } else if (length(nutrient_symbols)>1){
      if(exists("isoslice_tb")) rm("isoslice_tb")
      print(paste0("Isotopologues provided, but multiple tracer nutrients used. ",
                   "This is not supported currently, isotopologue data will be ignored"))
    } else {
      isoslice_tb<- sum_tb%>%
        
        make_iso_slices(factor_columns=factor_columns, tracer_column = tracer_column)
      
      #label factor name if any isotopologue has a significant difference, 
      #regardless of comparative factor level
      signi_iso_tb<-isoslice_tb%>%
        mutate(iso_sign_label=if_else(P_FC>=0.05|is.na(P_FC),
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
    # debug(make_FC_slices)
    # undebug(add_FClabels)
    
    FCslice_tb<- sum_tb%>%
      make_FC_slices(factor_columns=factor_columns, 
                     tracer_column = tracer_column,
                     abund_norm_comparlevel=abund_norm_comparlevel)%>%
      {
        if(P_isotopologues & isos_calculated) {
          left_join(.,signi_iso_tb) %>%
            mutate(iso_sign_label=if_else(is.na(iso_sign_label),
                                          "",
                                          iso_sign_label),
                   !!factor_symbol:=paste0(!!factor_symbol,iso_sign_label))%>%
            select(-iso_sign_label)
        } else {
          .
        }
      }
    
    return(FCslice_tb)
  }

generate_pies<-
  function(slice_tb,detail_charts,pathway_charts,savepath,
           normalize=T,factor_columns,tracer_column,factor_order,label_decimals,
           percent_add,
           FC_position,min_lab_dist,P_isotopologues,log_abund,circlelinecolor,
           circlelinetypes,maxcol_facet,include_name,col_labeling,
           alpha,otherfontsize,
           font,legendtitlesize,cohortsize,include_legend,
           mapotherfontsize=16,mapcohortsize=18,format="png",
           show_P=T,width=24.6,height=16) {
    
    #loop over each compound in input tibble
    compounds<-unique(slice_tb$compound)
    for (compound in compounds) {
      print(paste0("Processing compound ",which(compounds==compound),
                   " of ",length(compounds)))
      
      #prepare filename, remove problematic characters
      if (normalize) {
        plotfilename<-paste0("pies normalized ",compound,".",format)
      } else {
        plotfilename<-paste0("pies ",compound,".",format)
      }
      plotfilename<-gsub("/","-",plotfilename)
      
      if (detail_charts) {
        #plot detailed chart based on information in slice table
        print(paste0("saving detailed chart"))
        
        pies<-make_piechart(slice_tb,
                            factor_columns = factor_columns,
                            tracer_column = tracer_column,
                            log_abund=log_abund,
                            circlelinecolor = circlelinecolor,
                            selected_compound=compound,
                            circlelinetypes = circlelinetypes,
                            maxcol_facet = maxcol_facet,
                            include_name = include_name,col_labeling = col_labeling,
                            alpha=alpha,font=font,otherfontsize = otherfontsize,
                            legendtitlesize =legendtitlesize,
                            cohortsize = cohortsize,include_legend = include_legend,
                            show_P=show_P)
        
        #save detailed pie chart if required
        plotfilefolder<-paste0(savepath,"/Pie charts/")
        plotfilepath<-paste0(plotfilefolder,plotfilename)
        if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder),
                                                    recursive = T)
        ggsave(plotfilepath,pies,width=width,height=height,units = "cm",
               device = format)
      }
      
      if (pathway_charts) {
        print(paste0("saving pathway chart"))
        
        #plot summary pie chart for pathway based on information in slice table
        pies<-make_piechart(slice_tb,
                            factor_columns = factor_columns,
                            tracer_column = tracer_column,
                            log_abund=log_abund,
                            circlelinecolor = circlelinecolor,
                            selected_compound=compound,
                            circlelinetypes = circlelinetypes,
                            maxcol_facet = maxcol_facet,
                            include_name = F,col_labeling = col_labeling,
                            alpha=alpha,font=font,otherfontsize = mapotherfontsize,
                            legendtitlesize =mapcohortsize,
                            cohortsize = cohortsize,include_legend = F,
                            show_P=show_P)
        
        #save summary pie chart for pathway if required
        plotfilefolder<-paste0(savepath,"/Pie charts pathway/")
        plotfilepath<-paste0(plotfilefolder,plotfilename)
        if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder),
                                                    recursive = T)
        ggsave(plotfilepath,pies,width=24.6,height=16,units = "cm",
               device = format)
      }
      print("Finished")
    }
  }
  
#replaces last occurence of a specified character pattern by a replacement
#pattern that can differ in length
replace_lastchar<-function(rawstring,pattern,replacement) {
  #check if and where there is a match
  patternpos<-gregexpr(pattern,rawstring,fixed = T)[[1]]
  
  #if no match patternpos== -1 and original string can be returned
  if (patternpos[1]==-1) return(rawstring)
  finalstring<-paste0(
    substr(rawstring,1,max(patternpos)-1),
    replacement,
    substr(rawstring,max(patternpos)+1,nchar(rawstring))
  )
}

#Generate a figure caption based on settings
create_caption<-function(factor_order,log_abund,circlelinetypes,FC_position,show_P,
                         P_isotopologues) {
  #start and add factor level order
  caption<-
    paste0("Pie chart visualizations by Travis Pies applied to ",
    length(factor_order),
    " cohorts: ",
           replace_lastchar(
             paste0(factor_order,collapse = ", "),
             pattern = ",",
             replacement = " and"
           ),
           ". For each metabolite, the pie radii correspond to the relative ",
           "abundance which can be compared between the cohorts of this ",
           "metabolite. ")
  
  #Add abundance info depending on log scale being used and which concentric
  #circles are shown
  ncircles<-length(which(!circlelinetypes==0))
  if(log_abund){
    #add note using log scale
    caption<-paste0(
      caption,
      "The radii are plotted on a base 10 log scale. ")

    if(ncircles>0) {
      circlevals<-c(0.001,0.01,0.1,1)
      caption<-paste0(
        caption,
        if_else(ncircles>1,
                "The concentric circles correspond from center outwards to ",
                "The concentric circle corresponds to "),
        replace_lastchar(
          paste0(circlevals[!circlelinetypes==0],collapse = ", "),
          pattern = ",",
          replacement = " and"
          
        ),
        " times the largest abundance. "
      )
    }
    
  } else {
    if(ncircles>0) {
      circlevals<-c(0.25,0.5,0.75,1)
      caption<-paste0(
        caption,
        if_else(ncircles>1,
                "The concentric circles correspond from center outwards to ",
                "The concentric circle corresponds to "),
        replace_lastchar(
          paste0(circlevals[!circlelinetypes==0],collapse = ", "),
          pattern = ",",
          replacement = " and"
          
        ),
        " times the largest abundance. "
      )
    }
    
  }
  
  #About FC depending on label
  if(FC_position == "center") {
    caption<-paste0(
      caption,"Both the labeled surface fraction of the pie and the ",
      "percentage displayed in the middle of each pie reflect the fractional ",
      "contribution. "
    )
  } else if (FC_position == "slice") {
    caption<-paste0(
      caption,"Both the area of a slice and the ",
      "percentage displayed in it reflect the fractional ",
      "contribution of the corresponding source to the tracer element. "
    )
  }
  
  #about P values if displayed
  if(show_P) {
    caption<-paste0(
      caption,"pRA and pFC indicate the significance of the difference in ",
      "respectively the relative abundance (t-test) or fractional contribution ",
      "(Kruskal-Wallis) with ",
      factor_order[1],
      " (* indicates a p value <0.05). "
    )
  }
  
  #about isotopologue info if displayed
  if(P_isotopologues) {
    caption<-paste0(
      caption,"Cohorts with an * next to their name have at least one ",
      "significantly different isotopologue (Kruskal-Wallis p value <0.05, ",
      "isotopologues not shown in figure)"
    )
  }
  
  return(caption)
}