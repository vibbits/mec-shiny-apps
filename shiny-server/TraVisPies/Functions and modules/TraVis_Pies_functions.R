# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: Functions used by TraVis pies to convert raw input data to pie 
# charts. Called in modules of TraVis pies, but not inherently linked to R shiny
# functionality and could be used and useful outside shiny framework.

# Functions and libraries ---------------------------------------------------------------
#libraries for UI
library(dplyr)        #for faster.easier manipulation of data
library(tibble)       #for manipulating tibbles
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

#function to check escher-trace like corrected isotopologue input data
#returns "OK" if all checks are passed, error message otherwise
# tb<-read_csv_clean("~/GitHub/mec-shiny-apps/shiny-server/TraVisPies/Example_data/Input_Example_Isotopologues.csv",
#                         remove_empty = F,perc_to_num = F)
# head(tb)
# check_iso_input(tb)
# modtb<-tb
# modtb[2,2]<-"Abundance"
# head(modtb)
# check_iso_input(modtb)
# any(i(filter(modtb,Fragment=="Abundance")$Metabolite)==0)# debug(check_iso_input)

check_iso_input<-function(tb){
  if (colnames(tb)[1]!= "Metabolite") return(
    paste0("The first column in the isotopologue input should be named ",
           "Metabolite"))
  
  if (colnames(tb)[2]!= "Fragment") return(
    paste0("The second column in the isotopologue input should be named ",
           "Fragment"))
  if (any(filter(tb,nchar(Metabolite)>0)$Fragment!="Abundance")) return(
    paste0("At least one Metabolite entry is not in a row with abundance data.",
           " Only include these entries in rows containing your fragment ",
           "abundance marked by setting the Fragment entry to 'Abundance'"))
  if (any(is.na(filter(tb,Fragment=="Abundance")$Metabolite))) return(
    paste0("At least one abundance data row does not contain a Metabolite ",
           "entry. Always include these entries in rows containing your ",
           "fragment abundance"))
  return("OK")
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

#transpose tibble, setting colnames to first column and first column to colnames
t_tibble<-function(tb,first_colname="first_column"){
  
  #delete first column that will become the column names
  trans_tb<-select(tb,-1)%>%
    t()
  
  #Clean rownames and set colnames correctly, then save as tibble and add 
  #sample name column, 
  row.names(trans_tb)<-NULL
  colnames(trans_tb)<-pull(tb,1)
  trans_tb<-as_tibble(trans_tb) %>% 
    mutate(!!first_colname := colnames(tb)[-1],.before=1)
}

#Extract rowwise isotopologue data from columnwise corrected isotopologue
#file, for faster fractional contribution calculation and generation of 
#summarized isotopologue text for later
#specify correct isotopologue suffix separator, character used to separate the 
#metabolite name from the isotopologue label in the input isotopologue
#column names. This character can be used in metabolite name withotu issue, but 
#not in the isotopologue label
extract_col_isotopologues<-function(iso_col_tb,iso_suffix_sep="_") {
  #Add column with metabolite name extracted from isotopologue name based on
  #given suffix, then rename Isotopologues from 0 to highest isotopologue per 
  #metabolite
  iso_col_tb %>% t_tibble(first_colname = "Isotopologue") %>%
    #required to apply all functions (esp max) to current row only
    rowwise() %>%    
    mutate(Metabolite=
             substr(Isotopologue,1,
                    max(gregexpr(iso_suffix_sep,
                                 Isotopologue,fixed = T)[[1]])-1),
           .before=1) %>%
    group_by(Metabolite) %>%
    #n() gives the current group size
    mutate(Isotopologue=seq(from=0,to=n()-1,by=1)) %>%
    ungroup()
}

#Extract abundance data in columns from Escher-Trace like corrected isotopologue
#file
extract_et_abund<-function(iso_et_tb,sample_colname="Sample"){
  abund_tb<-filter(iso_et_tb,!is.na(Metabolite)) %>%
    t_tibble(first_colname = sample_colname) %>%
    slice(-c(1)) %>%
    mutate(across(!(!!sample_colname),.fns= as.numeric))
}

#Extract isotopologue data from Escher-Trace like corrected isotopologue
#file, keeping them in rows for faster fractional contribution calculation
#and generation of summarized isotopologue text for later
extract_et_isotopologues<-function(iso_et_tb){
  #Prepare isotopologue column
  iso_et_tb<-add_column(iso_et_tb,Isotopologue=NA,.before = 2)
  
  
  #note metabolite and isotopologue name for every isotopologue
  for (i in 1:nrow(iso_et_tb)) {
    if (!is.na(iso_et_tb$Metabolite[i])) {
      
      #if encountering new metabolite, save name and (re)set isotopologue
      #to 0 to prepare for first isotopologue coming up
      metabolite<-iso_et_tb$Metabolite[i]
      isotopologue<-0
    } else {
      #add metabolite name for this isotopologue in metabolite column and
      #and set isotopologue to current isotopologue index
      #for ordering purposes
      iso_et_tb$Metabolite[i]<-metabolite
      iso_et_tb$Isotopologue[i]<-isotopologue
      
      #increase isotopologue count
      isotopologue<-isotopologue+1
    }
  }
  #delete Fragment column and abundance entries
  iso_tb<-select(iso_et_tb,!c(Fragment)) %>%
    filter(!is.na(Isotopologue))
}


#Calculates a columnwise FC table based on an extracted isotopologue table
#can specify sample column name
calculate_FC<-function(iso_tb,sample_colname="Sample"){
  #calculate FC table, then reformat to columnwise format
  iso_tb %>% group_by(Metabolite) %>%
    summarise(across(!Isotopologue,
                     .fns = ~ sum(.x*Isotopologue)/max(Isotopologue))) %>%
    t_tibble(first_colname = sample_colname)
}

#Function to check samples across meta, abundance and FC tibbles and check
# compounds present. Outputs a list noting whether an error message should
# be given, and a message containing the error message or in absence
# of the error any warning messages to display
check_samples_compounds<-function(meta_tb,abund_tb,frac_tb,sample_column,
                                        norm_column){
  #set error =T by default, will change once past all error checks
  outlist<-list(error=T,message=NULL)
  
  #Check if all samples in meta table are present abund table
  samples_miss_abund<-!all(pull(meta_tb,sample_column) %in%
                             pull(abund_tb,sample_column))
  if (samples_miss_abund) {
    outlist$message<-paste0("Sample from metadata file missing in abund file. ",
                            "If the correct sample column is chosen, verify samples ",
                            "and sample names in both files.")
    return(outlist)
  }
  
  #Check if all samples in meta table are present frac table
  samples_miss_frac<-!all(pull(meta_tb,sample_column) %in%
                            pull(frac_tb,sample_column))
  if (samples_miss_frac) {
    outlist$message<-paste0("Sample from metadata file missing in frac file. ",
                            "If the correct sample column is chosen, verify samples ",
                            "and sample names in both files.")
    return(outlist)
  }
  
  #check if normalisation column is numeric if there is a column specified
  if (norm_column != "None") {
    if (!is.numeric(pull(meta_tb,norm_column))) {
      outlist$message<-paste0("The chosen normalisation column does not contain ",
                              "numbers. Please pick the right column or check the ",
                              "input if this is it.")
      return(outlist)
    }
  }
  
  #no errors encountered, set outlist$error to false
  outlist$error<-F
  
  #Warn if more samples present in abund or frac file than in meta
  samples_ignored<-!all( abund_tb[,sample_column] %in%
                           meta_tb[,sample_column],
                         frac_tb[,sample_column] %in%
                           meta_tb[,sample_column])
  if (samples_ignored) {
    outlist$message<-
      c(outlist$message,
        paste0("Samples from abund and or frac file missing in ",
               "metadatafile. These samples will be removed from the ",
               "analysis."))
  }
  #Warn if compounds present in abund file not frac file 
  #will be 100% unlabeled
  comp_ab_only<-colnames(abund_tb)[which(!colnames(abund_tb)%in%
                                           colnames(frac_tb))]
  if (length(comp_ab_only)>0) {
    outlist$message<-
      c(outlist$message,
        paste0("Following compounds only in abundance file, will be ",
               "considered fully unlabeled: ",
               paste(comp_ab_only,collapse = ", ")))
  }
  
  #Warn if compounds present in frac file not abund file
  #will be removed
  comp_fc_only<-colnames(frac_tb)[which(!colnames(frac_tb)%in%
                                          colnames(abund_tb))]
  if (length(comp_fc_only)>0) {
    outlist$message<-
      c(outlist$message,
        paste0("Following compounds only in fractional contribution file, ",
               " will be removed: ",
               paste(comp_fc_only,collapse = ", ")))
  }
  
  #Warn if compounds have 0 abundance in every sample, they will be dropped
  compounds_notdetected<-colnames(
    select(abund_tb,-where(has_nonzero))
  )
  
  if (length(compounds_notdetected)>0) {
    outlist$message<-
      c(outlist$message,
        paste0("Following compounds are never detected (abundance always 0), ",
               " and will be removed: ",
               paste(compounds_notdetected,collapse = ", ")))
  }
  
  #return empty text if no warnings, else give them in single orange text (html)
  if (length(outlist$message)>0) {
    outlist$message<-paste("<b><p style='color:orange'>Warning: </b>",
                           outlist$message,
                           "</p>", sep = "<br/>")
  } else {
    outlist$message<-""
  }
  return(outlist)
}

#Summarize extracted isotopologue data as a single row per sample 
#containing for each metabolite a string with all contributions from lowest to 
#highest isotopologue in order separated by |
summarize_isotopologue<-function(iso_tb,sample_colname="Sample"){
  #remove isotopologue column, group per metabolite and generate single string
  #per metabolite for each sample
  #then transpose from rowwise to columnwise representation
  iso_tb %>% select(-Isotopologue) %>%
    group_by(Metabolite) %>%
    summarize(across(.fns = ~ paste0(.x,collapse = "|"))) %>%
    t_tibble(first_colname = sample_colname)
}

merge_input<-function(meta_tb,abund_tb,frac_tb,iso_tb=NULL,
                          sample_col="Sample",compounds) {
  print(meta_tb)
  #Per compound adapt FC's below 0 (artefacts due to natural abundance
  #correction) to be positive to avoid problems with the visualisations
  #later on.
  for (i in (2:ncol(frac_tb))) {
    if (any(frac_tb[,i]<0)) {
      FCs<-pull(frac_tb[,i])
      FCs[which(FCs<0)]<-FCs[which(FCs<0)]-min(FCs[which(FCs<0)])     
    }
  }
  
  #modify iso_tb if it exists before summarizing
  if (length(iso_tb)>0) {
    #Per compound adapt isotopologues's below 0 (artefacts due to natural abundance 
    #correction) to be positive to avoid problems with the visualisations
    #later on
    for (i in (2:nrow(iso_tb))) {
      if (any(iso_tb[i,]<0)) {
        #check if any value for this isotopologue below 0
        metabolite<-iso_tb$Metabolite[i]
        isos<-iso_tb[i,-c(1,2)]
        negisos<-which(isos<0)
        
        #if no values negative, skip this section to avoid empty reference  
        #warnings and useless computing. If negatives, no zero correction was done
        #before and should be done now
        if (length(negisos)>0) {
          #Make variable containing negative iso value and 0 for others
          #then overwrite negative iso values to 0
          toadd<-isos
          toadd[-negisos]<-0
          isos[negisos]<-0   
          iso_tb[i,-c(1,2)]<-isos
          
          #add negative iso values to parent to offset previous addition to 
          #parent to compensate negative values
          parent_index<-which(iso_tb$Metabolite==metabolite & 
                                iso_tb$Isotopologue==0)
          iso_tb[parent_index,-c(1,2)]<-iso_tb[parent_index,-c(1,2)]+toadd
          
          #if any parents became <0, set to 0 (likely parent was undetectable)
          iso_tb[parent_index,][which(iso_tb[parent_index,]<0)]<-0
        }
        
      }
    }
    
    #summarize isotopologue data with name sample column
    iso_tb<-summarize_isotopologue(iso_tb,sample_colname = sample_col)
  }
  
  #rename sample column in all inputs
  meta_tb<-rename(meta_tb,Sample=all_of(sample_col))
  abund_tb<-rename(abund_tb,Sample=all_of(sample_col))
  frac_tb<-rename(frac_tb,Sample=all_of(sample_col))
  
  #add metadata to abundance and fractional contribution data respectively
  #retaining only selected samples, and drop metabolites with 0 abundance
  #in every sample to avoid errors
  abund_tb<-left_join(meta_tb,abund_tb,by="Sample") %>%
    select(1:ncol(meta_tb),any_of(compounds)) %>%
    select_if(has_nonzero)
  
  frac_tb<-left_join(meta_tb,frac_tb,by="Sample") %>%
    select(1:ncol(meta_tb),any_of(colnames(abund_tb))) 
  
  if(!length(iso_tb)==0) {
    iso_tb<-left_join(meta_tb,iso_tb,by="Sample") %>%
      select(1:ncol(meta_tb),any_of(colnames(abund_tb)))
  }
  
  #add fractional contribution and isotopologues equal to 100% unlabeled to 
  #compounds in abundance but not fraction labeling table
  if (any(!colnames(abund_tb) %in% colnames(frac_tb))) {
    nolabnames<-colnames(abund_tb)[which(! colnames(abund_tb) %in%
                                           colnames(frac_tb))]
    for (i in nolabnames) {
      frac_tb$new<-0
      colnames(frac_tb)[ncol(frac_tb)]<-i
    }
    if(!length(iso_tb)==0) {
      for (i in nolabnames) {
        iso_tb$new<-"1"
        colnames(iso_tb)[ncol(iso_tb)]<-i
      }
    }
  }
  

  #prepare abundance data for joining: 
  #calculate normalized abundances if normalization column provided and add
  #to abund tb as different datatype. 
  #add as character as isotopologue summaries will be character too
  abund_tb <-abund_tb %>% add_column(datatype="Abund")
  
  if ("Normalisation" %in% colnames(meta_tb)) {
    abund_tb<-abund_tb %>% 
      mutate(across((ncol(meta_tb)+1):(ncol(abund_tb)-1),
                    function(x) x/Normalisation)) %>%
      mutate(datatype="NormAbund") %>%
      full_join(abund_tb,by=colnames(abund_tb)) %>%
      mutate(across(any_of(compounds),as.character)) 
  } else {
    abund_tb<-abund_tb %>%mutate(across(any_of(compounds),as.character)) 
  }
  
  #prepare labeling  data for joining: 
  #Add isotopologue data to fractional contribution data
  frac_tb <-frac_tb %>% mutate(across(any_of(compounds),as.character)) %>%
    add_column(datatype="FracCont")
 
  print(frac_tb)
  print(iso_tb)
  if(!length(iso_tb)==0) {
    iso_tb$datatype<-"Isotopologues"
    frac_tb<-full_join(frac_tb,iso_tb,by=colnames(frac_tb))
  }

  tb<-full_join(frac_tb,abund_tb,by=colnames(abund_tb))

  
  #join all tables then order and remove normalisation factor if present
  tb<-full_join(frac_tb,abund_tb,by=colnames(abund_tb)) %>%
    select(colnames(meta_tb),datatype,everything())
  
  if ("Normalisation" %in% colnames(meta_tb)) {
    tb<-select(tb,-Normalisation) 
  }

  return(tb)
}
  

#Extract data for one compound in merged input, only for desired factor
#levels and set factor order
#need to use !! for dynamic variable names in tidyverse selection
#see https://stackoverflow.com/questions/50537164/summarizing-by-dynamic-column-name-in-dplyr 
obtain_compounddata<-function(tb,compound,fact_name,
                              fact_order=unique(pull(tb,!!fact_name)),
                              normalize=F){
  #prepare factor name symbol to use as target column name for mutate
  fact_symbol<-rlang::sym(fact_name)
  
  #select only one compound, filter to include normalized or non normalized
  #abundances, then select only given factor levels, then drops unused levels
  compound_tb<-tb %>% select(!!fact_name,datatype,!!compound) %>%
    filter(datatype %in% c("FracCont","Isotopologues",
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

#this function parses an isotopologue pattern string entry in a standardized
#TraVis file tomultiple rows each containing the contribution of one isotopologue
parse_isos_torow<-function(tb,valuecolumn) {
  #Separate isotopologue string entries from other entries
  isostring_tb<-filter(tb,datatype=="Isotopologues")
  out_tb<-filter(tb,!datatype=="Isotopologues")
  
  #for each isotopologue string entry, parse to one per isotopologue
  #add that to the output tibble
  for (i in 1:nrow(isostring_tb)) {
    iso_string<-pull(isostring_tb[i,valuecolumn])
    isos<-unlist(strsplit(iso_string,split = "|",fixed=T))
    iso_labels<-paste0("M",seq(0,length(isos)-1,1))
    iso_tb<-tibble(datatype=iso_labels,!!valuecolumn:=isos)%>%
      add_column()
    
    iso_tb<-isostring_tb %>% 
      slice(rep(i,length(isos))) %>%
      select(-datatype,-!!valuecolumn) %>%
      #use by=character() to cross join tibbles with no common column
      add_column(datatype=iso_labels,!!valuecolumn:=isos)
    out_tb<-bind_rows(out_tb,iso_tb)
  }
  
  return(out_tb)
}

summarize_addP<-function(tb,cohortcolumn,valuecolumn,
                         data_type=c("checkColumn","Abundance","FracCont",
                                     "Isotopologue")){
  #prepare cohort factor symbol to use for ordering
  fact_symbol<-rlang::sym(cohortcolumn)
  
  #if datatype is provided in column, sort tb per datatype to make sure order is
  # ok for rest of function. Otherwise, check if datatype provided as variable,
  # and add column with only that type.
  #If so, set to that datatype, if not, error.
  if (data_type=="checkColumn"){
    if ("datatype" %in% colnames(tb)) {
      tb<-tb[order(tb$datatype),]
    } else {
      print(paste0("summarize_P function requested to check for datatype in ",
                   "tibble column called 'datatype' (default option), but no such column ",
                   "provided. Either provide column name or specify datatype in function",
                   "call"))
    }
  } else {
    tb$datatype<-data_type
  }
  
  
  #initialize tibble for output with one entry per factor level each for all
  #datatypes, with initialized column for p values, and an index noting
  #the last row in the P column that received data
  tb_out<-unique(tb[,-which(colnames(tb)==valuecolumn)])
  tb_out$P<-NA
  index<-0
  
  #if only 1 cohort is provided, set P to 1 for further checking
  if (!length(unique(pull(tb[,cohortcolumn])))>1) {
    tb_out$P<-99
  } else {
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
        #t.test for abundance data and kruskal wallis for fraccont or iso
        #Set P=1 if all values are the same(likely 0) resulting in NaN. Make 
        #string depending on datatype
        if (length(tgtvalues)==1|length(refvalues)==1) {
          tb_out$P[index+i]<-99
        } else {
          if (datatype_selected=="Abund"){
            p<-t.test(refvalues,tgtvalues,)$p.value
            if (is.nan(p)) p<-1                 
            tb_out$P[index+i]<-p
          } else {
            p<-kruskal.test(c(refvalues,tgtvalues),
                            c(rep("Reference",length(refvalues)),
                              rep("Target",length(tgtvalues))))$p.value             
            if (is.nan(p)) p<-1                 
            tb_out$P[index+i]<-p
          }
        }
      }
      #raise index by amount of cohorts in last set
      index<-index+i
    }
  }
  
  #set datatypes to P labels to be output
  tb_out <- tb_out %>% mutate(datatype=case_when(
    tb_out$datatype=="Abund"     ~"P.RA",
    tb_out$datatype=="FracCont"  ~"P.FC",
    TRUE                           ~paste0("p",tb_out$datatype))) 
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
                             FC_position=="slice" & FracCont==0),"",labFC))%>%
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
#per slice. The  average abundance is normalized to the largest average abundance 
#and if desired this ratio is log10 transformed. This is the pie radius. 
#The fractions of the above parameter multiplied with the 
#labeled and unlabeled fraction correspond to the desired slices of a pie with
#this radius. P values for relative abundance and fractional contribution
#are calculated and a label for these on the pie chart is generated
prepare_slicedata<-function(compound_tb,compound,fact_name,
                            label_decimals,percent_add,FC_position,min_lab_dist,
                            P_isotopologues){
  #factor and compound need to be symbolized to use in 
  #tidyverse grouping function
  fact_symbol<-rlang::sym(fact_name)
  comp_symbol<- rlang::sym(compound) 
  
  #Calculates p values of significance tests of both relative abundance, and
  #fractional contribution per cohort factor level.
  P_tb<-compound_tb %>% select(!!fact_name,datatype,
                               !!compound)%>%
    summarize_addP(cohortcolumn = fact_name,valuecolumn = compound,
                   data_type = "checkColumn")%>%
    pivot_wider(names_from=datatype,values_from=P) 
    
  
  #get mean abundance and fractional contribution per cohort factor level
  #then join with P data from above.
  #Then format as table with single entry per cohort with abundance and FC as
  #variables. Use to calculate abundance normalized to biggest abundance, and if
  #desired transform this ratio by log10. Then calculate the part of this value
  #that is labeled, the part that is unlabeled, and finally the fractional 
  #contribution of the unlabeled part. Format as table with two entries
  #factor level, one for the labeled part and one for the unlabeled part
  print(sum_tb<-group_by(compound_tb,!! fact_symbol,datatype)%>%
          summarise(!!compound := mean(!!comp_symbol),.groups = "drop")%>%
          left_join(P_tb))
  sum_tb<-group_by(compound_tb,!! fact_symbol,datatype)%>%
    summarise(!!compound := mean(!!comp_symbol),.groups = "drop")%>%
    left_join(P_tb)%>%
    pivot_wider(names_from=datatype,values_from=!!compound)%>% 
    #repeat log abund n() times so amount matches amount of data entries as 
    #required by if_else()
    mutate(Abund=Abund/max(Abund),
           Labeled=FracCont*Abund,
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
    mutate(
      P.FClab=case_when(
        is.na(P.FC)       ~ "",
        P.FC==99          ~ "N=1,P=NA",
        P.FC<0.05         ~ paste0("pFC=",round(P.FC,2),"*"),
        P.FC>=0.05        ~ paste0("pFC=",round(P.FC,2))),
      P.RAlab=case_when(
        is.na(P.RA)       ~ "",
        P.RA==99          ~ "N=1,P=NA",
        P.RA<0.05         ~ paste0("pRA=",round(P.RA,2),"*"),
        P.RA>=0.05        ~ paste0("pRA=",round(P.RA,2)))
    ) %>%
    ungroup()
  
  #If required, add * to cohort name if any isotopologue P < 0.05
  #don't add anything if all are NA (reference cohort)
  if (P_isotopologues) {
    #make variable to store factor levels that have significant isotopologue
    #difference
    levels_orig<-levels(pull(slice_tb[,fact_name]))
    iso_cols<-colnames(slice_tb)[which(substr(colnames(slice_tb),1,2)=="pM")]
    for (i in levels_orig) {
      #get P's of isotopologues from first entry, 
      #check if any significant, add 1 to vector to avoid warnings
      iso_Ps<-unlist(slice_tb[
        which(slice_tb[,fact_name]==i & slice_tb$Labeling=="Unlabeled"),
        iso_cols],use.names = F)
      
      #add * to level name if significant
      if (min(c(iso_Ps,1),na.rm = T)<0.05) {
        newname<-paste0(i,"*")
        slice_tb<-slice_tb %>%
          mutate(!!fact_name:=recode(!!fact_symbol,!!i := newname))
      }
    }
  }
  
  return(slice_tb)
}

#makes pie chart based on table with required data per pie slice
make_piechart<-function(slice_tb,compound,
                        fact_name=fact_name,log_abund=F,
                        circlelinecolor="gray",maxcol_facet=4,
                        circlelinetypes=c(1,1,1,1),
                        include_name=F,col_labeling=c("#bfbfbf","#ffd966"),
                        alpha=0.7,
                        otherfontsize=10,font="sans",legendtitlesize=10,
                        cohortsize=12,include_legend=T,show_P=T){
  
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
    slice_tb <- slice_tb %>% 
      rowwise() %>%
      mutate(modAbund=10^((log10(Abund)+log10(minvalue))/2),
             modAbund_width=-(log10(minvalue)-log10(Abund)),
             #reset label distance position to work on logscale, either to middle
             #or to intended distance depending on whether one was given.
             FClab_posDist=if_else(FClab_posDist>0,
                                   log10(modAbund),
                                   log10(minvalue)))
    
    plotrect<-slice_tb %>% ggplot(aes(x = modAbund, y = Fraction, fill = Labeling, 
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
  }
  
  
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
                       normalize=T,fact_name,fact_order,label_decimals,
                       percent_add,FC_position,min_lab_dist,P_isotopologues,
                       log_abund,circlelinecolor,circlelinetypes,
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
  
  #either remove isotopologues or parse them into one entry per isotopologue
  #then make sure the value column is numeric for further analysis
  if (!P_isotopologues) {
    compound_tb<-filter(compound_tb,!datatype=="Isotopologues") %>%
      mutate(across(!!compound,as.numeric))
  } else {
    compound_tb<-parse_isos_torow(compound_tb,valuecolumn = compound) %>%
      mutate(across(!!compound,as.numeric))
  }
  
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
                              FC_position = FC_position,
                              P_isotopologues=P_isotopologues)
  
  if (detail_charts) {
    #plot detailed chart based on information in slice table
    print(paste0("saving detailed chart"))
    
    pies<-make_piechart(slice_tb,fact_name = fact_name,
                        log_abund=log_abund,
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
    if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder),
                                                recursive = T)
    ggsave(plotfilepath,pies,width=24.6,height=16,units = "cm",
           device = format)
  }
  
  if (pathway_charts) {
    print(paste0("saving pathway chart"))
    
    #plot summary pie chart for pathway based on information in slice table
    pies<-make_piechart(slice_tb,fact_name = fact_name,
                        log_abund=log_abund,
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
    if (!dir.exists(plotfilefolder)) dir.create(paste0(plotfilefolder),
                                                recursive = T)
    ggsave(plotfilepath,pies,width=24.6,height=16,units = "cm",
           device = format)
  }
  
}

# Generate pie chart plot for each compound and save if requested
generate_multiple_pies<-
  function(tb,compounds,detail_charts,pathway_charts,savepath,
           normalize=T,fact_name,fact_order,label_decimals,percent_add,
           FC_position,min_lab_dist,P_isotopologues,log_abund,circlelinecolor,
           circlelinetypes,maxcol_facet,include_name,col_labeling,
           alpha,otherfontsize,
           font,legendtitlesize,cohortsize,include_legend,
           mapotherfontsize=16,mapcohortsize=18,format="png",
           show_P=T) {
  
  #loop over each compound in input tibble
  for (compound in compounds) {
    print(paste0("Processing compound ",which(compounds==compound),
                 " of ",length(compounds)))
    generate_pie(tb=tb,compound=compound,detail_charts=detail_charts,
                 pathway_charts=pathway_charts,savepath=savepath,
                 normalize=normalize,fact_name=fact_name,fact_order=fact_order,
                 label_decimals=label_decimals,percent_add=percent_add,
                 FC_position=FC_position,min_lab_dist=min_lab_dist,
                 P_isotopologues=P_isotopologues,log_abund=log_abund,
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
create_caption<-function(fact_order,log_abund,circlelinetypes,FC_position,show_P,
                         P_isotopologues) {
  #start and add factor level order
  caption<-
    paste0("Pie chart visualizations by Travis Pies applied to ",
    length(fact_order),
    " cohorts: ",
           replace_lastchar(
             paste0(fact_order,collapse = ", "),
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
      fact_order[1],
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


# To keep for later --------------------------------------------------------------------
####salvage for function for generating isotopologue labels fit for ordering
# for (i in 1:nrow(iso_et_tb)) {
#   if (!is.na(iso_et_tb$Metabolite[i])) {
#     
#     #if encountering new metabolite, save name and (re)set isotopologue
#     #to 0 to prepare for first isotopologue coming up
#     metabolite<-iso_et_tb$Metabolite[i]
#     isotopologue<-0
#     
#     #count how many isotopologues for this metabolite, note how many digits
#     #for unanimous counting
#     count<-0
#     while(is.na(iso_et_tb$Metabolite[i+count+1]) & 
#           (i+count+1)<=nrow(iso_et_tb)) {
#       count<-count+1
#     }
#     maxdigit<-floor(log10(count))+1  
#   } else {
#     #add metabolite name for this isotopologue in metabolite column and
#     #and set isotopologue suffix to M0 if first isotopologue, and otherwise to 
#     #MX with X the current isotope number having the same amount of digits as 
#     #the highest isotopologue adding preceding 0's as neccesary
#     #for ordering purposes
#     iso_et_tb$Metabolite[i]<-metabolite
#     if (isotopologue==0) {
#       isotext<-paste0("_M",isotopologue)
#       
#     } else {
#       isotext<-paste0("_M",
#                       paste0(rep(0,maxdigit-(floor(log10(isotopologue))+1)),
#                              collapse = ""),
#                       isotopologue)
#     }
#     #join metabolite name and isotopologue suffix
#     iso_et_tb$Isotopologue[i]<- paste0(metabolite,isotext)
#     
#     #increase isotopologue count
#     isotopologue<-isotopologue+1
#   }
# }
# Test --------------------------------------------------------------------
#test crash when certain characters in names
#load library and metadata file, get metadata samples
#or spikes, and save additional variables as symbols for plot
# rawfolderpath<-r"(C:\GBW_MyPrograms\R_Apps_MEC\Apps\TraVis_Pies_v1.3\Rscripts\Example_data\Standardized input)" #folder with files, need this command to properly read in backslashes
# inputfile<-"Input_Example_standardized w isotopologues.csv"
# rawfolderpath<-r"(C:\GBW_MyPrograms\R_Apps_MEC\Apps\TraVis_Pies_v1.4\Rscripts\Example_data\Marco input failing)" #folder with files, need this command to properly read in backslashes
# inputfile<-"SDC_merged data.csv"
# folderpath<-gsub("\\\\", "/", rawfolderpath)         #get correct filepath from raw reference in input
# tb<-read_csv(paste0(folderpath,"/",inputfile))
# 
# compound<-colnames(tb)[4]
# fact_name<-"Cohort"
# fact_order<-pull(unique(tb[,fact_name]))
# normalize<-F
# P_isotopologues<-F
# 
# # first compound in inputtb
# v_settings<-list(compound=compound,
#                  fact_name=fact_name,
#                  fact_order=fact_order,
#                  norm=normalize,
#                  percent_add=F,
#                  FC_position="center",
#                  label_decimals=1,
#                  min_lab_dist=0.42,
#                  P_isotopologues=P_isotopologues,
#                  log_abund=T,
#                  circlelinecolor="gray",
#                  circlelinetypes=c(1,1,1,1),
#                  maxcol_facet=4,
#                  include_name=F,
#                  show_P=T,
#                  col_labeling=c("#bfbfbf","#ffd966"),
#                  alpha=0.7,
#                  otherfontsize=10,
#                  font="sans",
#                  legendtitlesize=10,
#                  cohortsize=12,
#                  include_legend=T)
# 
# out_settings<-list(plottype = "Stand-alone",
#                    format = "png",
#                    compounds = colnames(example_tb)[-c(1:3)])
# 
# format<-out_settings$format
# label_decimals<-v_settings$label_decimals
# min_lab_dist<-v_settings$min_lab_dist
# percent_add<-v_settings$percent_add
# FC_position<-v_settings$FC_position
# 
# #prepare filename
# if (normalize) {
#   plotfilename<-paste0("pies normalized ",compound,".",format)
# } else {
#   plotfilename<-paste0("pies ",compound,".",format)
# }
# 
# 
# 
# 
# #get table with only measured compound data, then a table summarizing
# #derived means and p values per cohort for abundance and one for fractional
# #contribution, then put together table with inputformat for pie function
# print(paste0("extracting compounddata"))
# 
# compound_tb<-obtain_compounddata(tb,compound,fact_name,
#                                  fact_order = fact_order,
#                                  normalize = normalize)
# 
# #either remove isotopologues or parse them into one entry per isotopologue
# #then make sure the value column is numeric for further analysis
# if (!P_isotopologues) {
#   compound_tb<-filter(compound_tb,!datatype=="Isotopologues") %>%
#     mutate(across(!!compound,as.numeric))
# } else {
#   compound_tb<-parse_isos_torow(compound_tb,valuecolumn = compound) %>%
#     mutate(across(!!compound,as.numeric))
# }
# 
# #make table with summarized data in the right format for pie creation
# #each entry containing the needed info for one slice of one of the pie 
# #charts.The average abundance normalized to the largest average abundance 
# #is the pie radius. The fractions of the above parameter multiplied with the
# #labeled and unlabeled fraction correspond to the desired slices of a pie 
# #with this radius 
# print(paste0("preparing slice data"))
# 
# slice_tb<-prepare_slicedata(compound_tb,fact_name = fact_name,
#                             compound=compound,label_decimals = label_decimals,
#                             min_lab_dist = min_lab_dist,
#                             percent_add = percent_add,
#                             FC_position = FC_position,
#                             P_isotopologues=P_isotopologues)
# 
# 

# #test isotopologue adapated functions
# example_tb<-read_csv(
#   file = "~/GitHub/mec-shiny-apps/shiny-server/TraVisPies/Example_data/Standardized input/Input_Example_standardized w isotopologues.csv")
# inputtb<-example_tb[,c(2:4)] %>%
#   filter(!datatype=="Abund")
# head(inputtb)
# 
# # Coenzyme_A
# v_settings<-list(compound="Coenzyme_A",
#                  fact_name=colnames(inputtb)[1],
#                  fact_order=pull(unique(inputtb[,1])),
#                  norm=T,
#                  percent_add=F,
#                  FC_position="center",
#                  label_decimals=1,
#                  min_lab_dist=0.42,
#                  P_isotopologues=T,
#                  log_abund=T,
#                  circlelinecolor="gray",
#                  circlelinetypes=c(1,1,1,1),
#                  maxcol_facet=4,
#                  include_name=F,
#                  show_P=T,
#                  col_labeling=c("#bfbfbf","#ffd966"),
#                  alpha=0.7,
#                  otherfontsize=10,
#                  font="sans",
#                  legendtitlesize=10,
#                  cohortsize=12,
#                  include_legend=T)
# 
# out_settings<-list(plottype = "Stand-alone",
#                    format = "png",
#                    compounds = colnames(example_tb)[-c(1:3)])
# 
# compound<-v_settings$compound
# 
# comptb_sumiso<-obtain_compounddata(
#   example_tb,compound=v_settings$compound,fact_name = v_settings$fact_name,
#   normalize=v_settings$norm)
# 
# comptb<-parse_isos_torow(comptb_sumiso,valuecolumn = v_settings$compound) %>%
#   mutate(across(v_settings$compound,as.numeric))
# 
# # head(comptb)
# # comptb[comptb$Cohort=="10uM AMA"&!comptb$datatype %in% c("Abund","FracCont"), "Phosphoenolpyruvic_acid"]<-
# #   comptb[comptb$Cohort=="NT"&!comptb$datatype %in% c("Abund","FracCont"), "Phosphoenolpyruvic_acid"]
# 
# 
# slicetb<-prepare_slicedata(comptb,fact_name = v_settings$fact_name,
#                       compound=v_settings$compound,
#                       label_decimals = v_settings$label_decimals,
#                       min_lab_dist = v_settings$min_lab_dist,
#                       percent_add = v_settings$percent_add,
#                       FC_position = v_settings$FC_position,
#                       P_isotopologues=v_settings$P_isotopologues)
# 
# #makes pie chart based on table with required data per pie slice
# (test<-make_piechart(slicetb,fact_name = v_settings$fact_name,
#                      log_abund = v_settings$log_abund,
#                      compound=v_settings$compound,
#                      circlelinecolor=v_settings$circlelinecolor,
#                      maxcol_facet=v_settings$maxcol_facet,
#                      circlelinetypes=v_settings$circlelinetypes,
#                      include_name=v_settings$include_name,
#                      col_labeling=v_settings$col_labeling,
#                      alpha=v_settings$alpha,
#                      otherfontsize=v_settings$otherfontsize,
#                      font=v_settings$font,
#                      legendtitlesize=v_settings$legendtitlesize,
#                      cohortsize=v_settings$cohortsize,
#                      include_legend=v_settings$include_legend,
#                      show_P=v_settings$show_P))
# 
# 
# #prepare folder to save to, and set which charts to generate depending on
# #requested plottype
# target_savepath<-paste0(getwd(),"/temp")
# detail_charts<-F
# pathway_charts<-F
# filelist<-NULL
# 
# if (out_settings$plottype %in% c("Stand-alone","Both")) {
#   detail_charts<-T
#   filelist<-c(filelist,"Pie charts/")
# }
# if (out_settings$plottype %in% c("Pathway-compatible","Both")) {
#   pathway_charts<-T
#   filelist<-c(filelist,"Pie charts pathway/")
# }
# 
# generate_pie(example_tb,compound=compound,detail_charts=detail_charts,
#              pathway_charts=pathway_charts,savepath=target_savepath,
#              normalize=v_settings$norm,fact_name=v_settings$fact_name,
#              fact_order=v_settings$fact_order,
#              P_isotopologues=v_settings$P_isotopologues,
#              log_abund=v_settings$log_abund,
#              label_decimals=v_settings$label_decimals,
#              percent_add = v_settings$percent_add ,
#              FC_position=v_settings$FC_position,
#              min_lab_dist =v_settings$min_lab_dist,
#              circlelinecolor=v_settings$circlelinecolor,
#              circlelinetypes=v_settings$circlelinetypes,
#              maxcol_facet=v_settings$maxcol_facet,
#              include_name=v_settings$include_name,
#              show_P=v_settings$show_P,
#              col_labeling=v_settings$col_labeling,
#              alpha=v_settings$alpha,
#              otherfontsize=v_settings$otherfontsize,
#              font=v_settings$font,
#              legendtitlesize=v_settings$legendtitlesize,
#              cohortsize=v_settings$cohortsize,
#              include_legend=v_settings$include_legend,
#              format=out_settings$format)
# 
# generate_multiple_pies(example_tb,compounds=out_settings$compounds,
#                        detail_charts=detail_charts,
#                        pathway_charts=pathway_charts,savepath=target_savepath,
#                        normalize=v_settings$norm,fact_name=v_settings$fact_name,
#                        fact_order=v_settings$fact_order,
#                        P_isotopologues=v_settings$P_isotopologues,
#                        log_abund=v_settings$log_abund,
#                        label_decimals=v_settings$label_decimals,
#                        percent_add = v_settings$percent_add ,
#                        FC_position=v_settings$FC_position,
#                        min_lab_dist =v_settings$min_lab_dist,
#                        circlelinecolor=v_settings$circlelinecolor,
#                        circlelinetypes=v_settings$circlelinetypes,
#                        maxcol_facet=v_settings$maxcol_facet,
#                        include_name=v_settings$include_name,
#                        show_P=v_settings$show_P,
#                        col_labeling=v_settings$col_labeling,
#                        alpha=v_settings$alpha,
#                        otherfontsize=v_settings$otherfontsize,
#                        font=v_settings$font,
#                        legendtitlesize=v_settings$legendtitlesize,
#                        cohortsize=v_settings$cohortsize,
#                        include_legend=v_settings$include_legend,
#                        format=out_settings$format)

# old tests ---------------------------------------------------------------
# print(create_caption(fact_order = v_settings$fact_order,
#                log_abund = v_settings$log_abund,
#                FC_position = v_settings$FC_position,
#                show_P = v_settings$show_P,
#                P_isotopologues = v_settings$P_isotopologue))

#test pie chart function for logscale abundance ratios
# fact_name<-"Cohort"
# circlelinecolor="gray"
# maxcol_facet=4
# circlelinetypes=c(2,1,1,1)
# include_name=F
# col_labeling=c("#bfbfbf","#ffd966")
# alpha=0.7
# otherfontsize=10
# font="sans"
# legendtitlesize=10
# cohortsize=12
# include_legend=T
# show_P=T
# min_lab_dist=0.42
# 
# #set minimal value to include on log axis, changing not recommended
# #and calculate minimal position distance on new sclae
# minvalue<-0.0001
# lab_dist_mod<-10^((log10(min_lab_dist)+log10(minvalue))/2)
# 
# 
# 
# if_else(slice_tb$FClab_posDist>0,
#         max(log10(lab_dist_mod),log10(slice_tb$modAbund)),
#         minvalue)
# 
# #width can only be symmetric, so modify abundance to the value on normal 
# #scale to correspond to the value in the middle 
# #between log scale minimal limit and logscale actual abundance, 
# #and abundance width to be the double of this modified abundance
# slice_tb <- tb %>% 
#   rowwise() %>%
#   mutate(modAbund=10^((log10(Abund)+log10(minvalue))/2),
#          modAbund_width=-(log10(minvalue)-log10(Abund)),
#          #reset label distance position to work on logscale, either to middle
#          #or to intended distance depending on whether one was given.
#          FClab_posDist=if_else(FClab_posDist>0,
#                                max(log10(lab_dist_mod),log10(slice_tb$modAbund)),
#                                log10(minvalue)))
# 
# 
# plotrect<-slice_tb %>% ggplot(aes(x = modAbund, y = Fraction, fill = Labeling, 
#                                   width = modAbund_width)) + 
#   scale_x_log10(limits= c(minvalue, 1)) +
#   geom_vline(xintercept=c(0.001),colour=circlelinecolor,
#              linetype=circlelinetypes[1])+ 
#   geom_vline(xintercept=c(0.01),colour=circlelinecolor,
#              linetype=circlelinetypes[2])+ 
#   geom_vline(xintercept=c(0.1),colour=circlelinecolor,
#              linetype=circlelinetypes[3])+ 
#   geom_vline(xintercept=c(1),colour=circlelinecolor,
#              linetype=circlelinetypes[4])+ 
#   geom_bar(stat = "identity", position = "fill",alpha=alpha) 
# 
# 
# 
# #add name of compound if desired, and the assign colors and their legend order
# if (include_name) plotrect<-plotrect+ggtitle(compound)
# plotrect<-plotrect  +
#   scale_fill_manual(values=col_labeling,guide=guide_legend(reverse=T))
# 
# #positions of text at specified locations, if desired.
# #Fontsize needs to be adjusted for reasons:
# #https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
# plotrect<-plotrect  +
#   geom_text(aes(label=labFC),x = slice_tb$FClab_posDist,
#             y=slice_tb$FClab_posAngle,size=otherfontsize*5/14)
# if (show_P) {
#   plotrect<-plotrect  +
#     geom_text(aes(label=P.RAlab),x=1.6,y=7/8,size=otherfontsize*5/14,
#               hjust="inward",vjust="inward") +
#     geom_text(aes(label=P.FClab),x=1.6,y=5/8,size=otherfontsize*5/14,
#               hjust="inward",vjust="inward")
# }
# 
# #transform bar to pie chart and plot pies on grid.
# piebasic<-plotrect+
#   facet_wrap(vars(!!rlang::sym(fact_name)),ncol=maxcol_facet) +
#   coord_polar("y", start = 0, direction = 1)
# 
# 
# #apply final formatting to pie plots. Removes x and y labels entirely,
# #including the space reserved for them on the plot
# #sets relative abundance p values in upper right corner of pie plots
# pies<-piebasic +
#   labs(x=NULL, y=NULL)+
#   #Change plots to black on white, remove text axes (fraction) that interfere
#   #with circles, axis ticks, fraction grid lines. Set text font,
#   #set legend title size,
#   #remove rectangles and background around factor levels, set factor levels
#   #to right text size
#   theme_bw(base_size = otherfontsize) +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid = element_blank(),
#         text=element_text(family = font),
#         plot.title = element_text(size = cohortsize, face = "bold"),
#         legend.title = element_text(size = legendtitlesize),
#         strip.background = element_rect(fill = NA, colour = NA),
#         strip.text = element_text(size = cohortsize))
# 
# #removes legend if desired
# if (!include_legend) pies <-pies + theme(legend.position = "none")
# 
# pies






#test input merging -----------------------------------------------------------------------
meta_tb<-read_csv_clean(paste0(getwd(),
                               "/Example_data/Original input/Input_Example_metadata.csv"),
                        remove_empty = T,perc_to_num = F)
# meta_tb[,3]<-1
abund_tb<-read_csv_clean(paste0(getwd(),
                                "/Example_data/Original input/Input_Example_RA.csv"),
                         remove_empty = T,perc_to_num = F)
iso_et_tb<-read_csv_clean(paste0(getwd(),
                                 "/Example_data/Original input/Input_Example_RA.csv"),
                          remove_empty = F,perc_to_num = F)
iso_col_tb<-read_csv_clean(paste0(getwd(),
                                  "/Example_data/Original input/Input_Example_isotopologues.csv"),
                           remove_empty = T,perc_to_num = T)
iso_tb<-extract_col_isotopologues(iso_col_tb) %>%
    slice(-c(5,6,7,8))

# iso_tb<-extract_et_isotopologues(iso_et_tb) %>%
#   slice(-c(5,6,7,8))

# abund_tb<-extract_et_abund(iso_et_tb,sample_colname = "Sample")

frac_tb<-calculate_FC(iso_tb,sample_colname = "Sample")

sample_col<-"Sample"

compounds<-colnames(abund_tb)[-1]
head(meta_tb)
(meta_formatted_tb<-format_metadata(meta_tb,sample_column = "Sample",
                                   factor_column = "Cohort",
                                   norm_column = "None"))

test<-merge_input(meta_tb = meta_formatted_tb,
                  abund_tb = abund_tb,
                  frac_tb = frac_tb,
                  iso_tb=iso_tb,
                  sample_col = sample_col,
                  compounds = compounds)

abund_tb<-abund_tb %>%mutate(across(-c(1:3),as.character)) 




