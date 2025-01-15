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


# Functions for data curation and merging ---------------------------------------------

# function for checking if any column cell is different from 0
has_nonzero <- function(x) { any(x != 0)}         

#replaces all occurences of a character in a string except last.
#Useful to make sure compound names match between isotopologue data and other types
replace_except_last <- function(input_strings, to_replace = "_", replacement = " ") {
  # Function to handle a single string
  replace_single_string <- function(input_string) {
    # Find all positions of the character to replace
    positions <- gregexpr(to_replace, input_string, fixed = TRUE)[[1]]
    
    # If there is only one or no occurrences, return the string as is
    if (length(positions) <= 1) {
      return(input_string)
    }
    
    # Replace all occurrences except the last one with the replacement character
    for (i in seq_along(positions)[-length(positions)]) {
      substr(input_string, positions[i], positions[i]) <- replacement
    }
    
    return(input_string)
  }
  
  # Apply the function to each string in the vector
  sapply(input_strings, replace_single_string)
}


#function to check escher-trace like corrected isotopologue input data
#returns "OK" if all checks are passed, error message otherwise
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

#function that return a symbol or a symbol list from a value, length >1, or 
#NULL if NULL was provided
sym_or_null<-function(vect,returnlist=F,allownull=T){
  if(length(vect)>1) {
    if(!returnlist){
      warning(paste0("Returning symbollist instead of the requested single ",
                   "symbol as vect contains multiple elements: ",
                   paste0(vect, collapse= ", ")))
      }
    symbol<-rlang::syms(vect)
  } else if(length(vect)==1) {
    if (returnlist) symbol<-rlang::syms(vect) else symbol<-
        rlang::sym(vect)
  } else if (allownull) symbol<-NULL else {
    stop("NULL input not allowed, check input or set allownull=T")
  }
  
  return(symbol)
}

#extract metadata from an excel sheet
extract_excelsheet_tb<-function(excelpath,sheetnamestring,datatype_intended,
                                samplename="Sample"){
  sample_symbol<-sym_or_null(samplename,allownull = F)
  
  sheetname<-excel_sheets(excelpath)[
    which(grepl(sheetnamestring,tolower(excel_sheets(excelpath))))]
  if(length(sheetname)==0) {
    stop(paste0("No sheets with sheetnamestring in its name."))
  }
  if(length(sheetname)>1) {
    print(paste0("Multiple sheets with ",
                 sheetnamestring," in their name, took first."))
    sheetname<-sheetname[1]
  }
  if(length(datatype_intended)>0) {
    print(paste0("Name of sheet used for ",datatype_intended,
                 ": ",sheetname))
  } else print(paste0("Name of sheet used: ",sheetname))
  
  
  #read excelsheet, rename sample column
  read_excel(excelpath,sheetname)%>%
    rename(!!sample_symbol:=1)
}

#function to prepare metadata to uniform format
format_metadata<-function(meta_tb,sample_column,factor_columns=NULL,
                          norm_column=NULL,tracer_column=NULL,
                          sampletype_column=NULL) {
  #make symbols for dplyr pipelines
  sample_symbol<-sym_or_null(sample_column,allownull = F)
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  tracer_symbol<-sym_or_null(tracer_column)
  norm_symbol<-sym_or_null(norm_column)
  sampletype_symbol<-sym_or_null(sampletype_column)
  
  #set right column types and keep only desired columns, if existing
  meta_tb<-mutate(meta_tb,
                  across(
                    c(!!sample_symbol,!!tracer_symbol,!!sampletype_symbol),
                    as.character),
                  across(
                    !!norm_symbol,
                    as.numeric),
                  across(
                    c(!!!factor_symbols),
                    as.factor)
                  )%>%
    select(where(~ !all(is.na(.))),!!sample_symbol,!!!factor_symbols,
           !!tracer_symbol,!!sampletype_symbol,!!norm_symbol)
}

#look for .csv file containing string in folder
loadfile_stringmatch<-function(path,filestring){
  #get .csv files that contain filestring
  filelist<-list.files(path)
  filelist<-filelist[which(grepl(".csv",filelist,fixed=T))]
  matchedfiles<-filelist[which(grepl(filestring,filelist,fixed=T))]
  if(length(matchedfiles)>1) {
    warning(paste0("Multiple .csv files contain the string `",filestring,
                 "`. Only the first of these wille be used: ",
                 paste0(matchedfiles,collapse=", ")))
  } else if (length(matchedfiles)==0) {
    stop(paste0("None of the .csv files in this folder contain the string `",
                filestring,"`. .csv files in this folder: ",
                paste0(filelist,collapse=", ")))
  }
  
  vroom::vroom(file = paste0(path,"/",matchedfiles[1]),
               delim = ",",show_col_types = FALSE)
}

#detect if labeling is isotopologues or fractional contribution,
#sanitazing metabolite names as required for later functions depending on
#labeling type
is_isodata<-function(label_tb,isostring=NULL){
  if(length(isostring)==0){
    print(paste0("No isostring provided, assuming labeling data is fractional",
    " contribution."))
    return(F)
  }
  isocols<-colnames(label_tb)[
    grepl(tolower(isostring),tolower(colnames(label_tb)))]
  if(length(isocols)>0){
    print(paste0("The isotopologuestring `",isostring,"` was detected in ",
                 "these columns (max first 10 shown): ",
                 paste0(isocols[1:min(10,length(isocols))],collapse = ", "),
                 "."))
    print(paste0("If this is fractional contribution data, set isostring to ",
                 "a different string not present in metabolite names to avoid ",
                 "errors."))
    return(T)
  } else {
    print(paste0("The isostring `",isostring,"` is not found ",
                 "in any column names in the labeling tb, ",
                 "assuming this is fractional contribution ",
                 "data. If this is isotoplogue data, set ",
                 "isostring to a different string unique ",
                 "to isotopologue names to avoid errors"))
    return(F)
  }
}

#rename compounds in colnames (incl isotopologue data) based on library
#ideal to rename derivatised compounds to the original compound
rename_lib<-function(data_tb,lib_tb=NULL,currentnamecol="compound",
                     newnamecol="Orig_name"){
  if(length(lib_tb)==0) {
    print("No library provided, using compound names as they are in input")
    return(data_tb)
  }
  
  if(any(colnames(lib_tb)==currentnamecol)&
     any(colnames(lib_tb)==newnamecol)) {
    for(i in 1:nrow(lib_tb)){
      colnames(data_tb)<-sub(lib_tb[i,currentnamecol],
                              lib_tb[i,newnamecol],
                              colnames(data_tb),
                              fixed = T)
    }  
  } else warning(paste0("Renaming requested by naming columns not found in ",
                        "library, continue without renaming"))
  
  return(data_tb)
}

#check inputpath to see if leading to excel or folder, import meta, abundance,
#and labeling data accordingly, and return a list of these three datasets as
#tibble
list_inputdata_tbs<-function(inputpath,metastring="meta",abundstring="abund",
                             labelstring="iso",isostring="parent",sample_column,
                             factor_columns=NULL,
                             norm_column=NULL,tracer_column=NULL,
                             sampletype_column=NULL,
                             lib_tb=NULL){
  
  
  #correctly load in data depending on path being excel or folder
  if(grepl(".xls",inputpath)) {
    print(
      paste0("String `.xls` detected in ",inputpath,", expecting excelfile"))
    if (!file.exists(inputpath)) {
      stop(paste0(inputpath," is not an existing excel file, correct or remove ",
                  "`.xls` string from inputpath if it is a path to a folder"))
    }
    # debug(format_metadata)
    meta_tb<-extract_excelsheet_tb(inputpath,
                                             sheetnamestring = metastring,
                                             datatype_intended = "metadata",
                                             samplename = sample_column)
    
    abund_tb<-extract_excelsheet_tb(inputpath,
                                    sheetnamestring = abundstring,
                                    datatype_intended = "abundance data",
                                    samplename = sample_column)
    
    label_tb<-extract_excelsheet_tb(inputpath,
                                    sheetnamestring = labelstring,
                                    datatype_intended = "labeling data",
                                    samplename = sample_column)
    
    
    
    
  } else if(dir.exists(inputpath)) {
    meta_tb<-loadfile_stringmatch(inputpath,metastring)
    abund_tb<-loadfile_stringmatch(inputpath,abundstring)
    label_tb<-loadfile_stringmatch(inputpath,labelstring)
    
  } else stop(paste0(inputpath ,"not found. Please specify an existing folder or ",
                     "excel file. For the latter, make sure the path contains ",
                     "the string `.xls` somehwere, like ",
                     "examplefolder/examplefile.xlsx"))
  
  #start list with formatted metadata and abundance data with sanitized
  #metabolite names (e.g. removing "_" from metabolite names to avoid 
  #confusion with isotopologues). Remove empty columns and rows too.
  #rename compounds by library if provided
  output_list<-list(meta_tb=meta_tb%>%
                      format_metadata(sample_column = sample_column,
                                      factor_columns = factor_columns,
                                      norm_column = norm_column,
                                      tracer_column=tracer_column,
                                      sampletype_column=sampletype_column),
                    abund_tb=abund_tb%>%
                      rename_lib(lib_tb)%>%
                      rename_with(~gsub("_"," ",.x))%>%
                      select(where(~ !all(is.na(.))))%>%
                      na.omit()
  )
  
  #clean labeling data to remove empty rows and columns, and make sure
  #labeling is a numerical fraction in case input is a % character
  #rename compounds by library if provided
  label_tb<-label_tb%>%
    select(where(~ !all(is.na(.))))%>%
    na.omit()%>%
    mutate(across(where(~any(grepl("%",.x,fixed=T))),function(x) 
      as.numeric(sub(pattern="%", replacement = "",x,fixed = T))/100))%>%
    rename_lib(lib_tb)
  
  #rename labeling data to FC or iso and clean names accordingly depending on
  #whether an isotopologue string is found in the column names
  if(is_isodata(label_tb,isostring = isostring)) {
    output_list$iso_tb<-label_tb %>%
      rename_with(replace_except_last)
  } else {
    output_list$frac_tb<-label_tb %>%
      rename_with(~gsub("_"," ",.x))
  }
  return(output_list)
}

#reads in an abundance, fractional contribution or isotopologue sheet of a 
#dolly excel file by name,and prepares the desired table from it
curate_abundancedata<-function(abund_tb,meta_tb,sample_column,sampletype_column,
                               norm_column) {
  sample_symbol<-sym_or_null(sample_column)
  sampletype_symbol<-sym_or_null(sampletype_column)
  norm_symbol<-sym_or_null(norm_column)
  
  #detect internal standards in sheet
  headers<-abund_tb %>% 
    select(where(~ all(is.na(.)))) %>%
    colnames()
  
  if (length(headers[which(grepl("internal",tolower(headers)))])>0) {
    intstdfirstcol<-which(colnames(abund_tb) ==
                            headers[which(grepl("internal",
                                                tolower(headers)))][1])+1
    intstdlastcol<-which(colnames(abund_tb) ==
                           headers[which(grepl("internal",
                                               tolower(headers)))+1][1])-1
    intstds<-colnames(abund_tb)[intstdfirstcol:intstdlastcol]
  } else {
    intstds<-NULL
  }  
  
  #remove empty rows and columns and check if abundance sheet for specific actions to 
  #take with it. Convert data to long format if not done before.
  abund_cleantb <- abund_tb %>% 
    select(where(~ !all(is.na(.))))%>%
    na.omit()
  
  #if blanks are present in metadata
  blankspresent<-F
  if (length(sampletype_column)>0) {
    if(any(tolower(pull(meta_tb,!!sampletype_symbol))=="blank")) blankspresent<-T
  }
  if(blankspresent) {
    
    #calculate LOD and average blank from blanks samples, Set blank and LOD to 0
    #for internal standards
    blanks_tb<-abund_cleantb %>%
      filter(!!sample_symbol %in% c(meta_tb %>%
                                      filter(tolower(
                                        pull(meta_tb,!!sampletype_symbol))=="blank")%>%
                                      pull(!!sample_symbol))) %>%
      pivot_longer(cols = 2:ncol(.),
                   names_to = "compound", 
                   values_to= "Abund")%>%
      group_by(compound)%>%
      summarise(av_blank=mean(Abund),
                LOD=av_blank+3*sd(Abund))%>%
      mutate(
        av_blank=if_else(compound %in% intstds,0,av_blank),
        LOD=if_else(compound %in% intstds,0,LOD)
      )
    
    #join blank and LOD data to the abundance data, do blank correction
    #and check if compound above LOD. Then remove rows for blanks and LOD from 
    #abundance data.
    abund_cleantb <- abund_cleantb %>%
      pivot_longer(cols = 2:ncol(.),
                   names_to = "compound", 
                   values_to= "Abund")%>%
      left_join(blanks_tb, by=join_by(compound))%>%
      right_join(meta_tb%>%
                   select(!!sample_symbol,!!sampletype_symbol,!!norm_symbol),
                 join_by(!!sample_symbol))%>%
      mutate(BlankcorAbund=Abund-av_blank,
             LOD_blankcor=LOD-av_blank,
             detected=Abund>LOD)%>% 
    pivot_longer(cols = c("Abund","BlankcorAbund"),
                 names_to = "datatype", 
                 values_to= "value") %>%
      select(!!sample_symbol,!!sampletype_symbol,!!norm_symbol,compound,datatype,value,LOD,LOD_blankcor,
             detected)%>%
      filter(!!sampletype_symbol!="blank",
             !grepl("lod",tolower(!!sample_symbol)))%>%
      pivot_wider(names_from = datatype,values_from = value)
  } else {
    #if no mocks, don't do blank correction
    print(paste0("No samples indicated as blank in metadata column sampletype ",
                 ". Assumed no LOD calculation or blank correction needed"))
    abund_cleantb <- abund_cleantb %>%
      pivot_longer(cols = 2:ncol(.),
                   names_to = "compound", 
                   values_to= "Abund")%>%
      right_join(meta_tb%>%
                   select(!!sample_symbol,!!norm_symbol),
                 join_by(!!sample_symbol))
  }
  
  #if norm_column is not null, add normalized data based on blank corrected
  #if available, not blank corrected if not
  if(length(norm_column)>0) {
    abund_cleantb<-abund_cleantb%>%
      {
        if("BlankcorAbund" %in% colnames(.)){
          mutate(.,NormAbund=BlankcorAbund/!!norm_symbol)
        } else {
          mutate(.,NormAbund=Abund/!!norm_symbol)
        }
      }
  }
  return(abund_cleantb)
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
#file, to use for calculating fraccon and making long tibble
#specify correct isotopologue suffix separator, character used to separate the 
#metabolite name from the isotopologue label in the input isotopologue
#column names. This character can no be present in the metabolite name or  
#isotopologue label, only to separate the name and the isotopologue label
extract_col_isotopologues<-function(iso_col_tb,iso_suffix_sep="_") {
  #Add column with metabolite name extracted from isotopologue name based on
  #given suffix, then rename Isotopologues from 0 to highest isotopologue per 
  #metabolite
  iso_col_tb %>% 
    select(where(~ !all(is.na(.))))%>%
    t_tibble(first_colname = "Isotopologue") %>%
    rowwise() %>%    
    mutate(compound=
             substr(Isotopologue,1,
                    max(gregexpr(iso_suffix_sep,
                                 Isotopologue,fixed = T)[[1]])-1),
           datatype=paste0("Isotopologue_",
                           substr(Isotopologue,
                                  max(gregexpr(iso_suffix_sep,
                                               Isotopologue,fixed = T)[[1]])+1,
                                  nchar(Isotopologue))),
           .before=1) %>%
    group_by(compound) %>%
    #n() gives the current group size
    mutate(Isotopologue=seq(from=0,to=n()-1,by=1)) %>%
    ungroup()%>%
    pivot_longer(4:ncol(.),names_to = "Sample",values_to = "value") %>%
    select(Sample,everything())
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
  sample_symbol<-rlang::sym(sample_colname)
  
  #calculate FC table, dropping values lower than 0
  # then reformat to columnwise format
  iso_tb %>%
    filter(!value<0) %>%
    group_by(compound,!!sample_symbol) %>%
    summarise(value = sum(value*Isotopologue)/max(Isotopologue))%>%
    select(!!sample_symbol,everything())%>%
    pivot_wider(names_from = compound,values_from = value)
    
  # %>%
  #   t_tibble(first_colname = sample_colname)
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

# merge_input<-function(meta_tb,abund_tb,frac_tb,iso_tb=NULL,
#                                       sample_col="Sample",compounds) {
#   #Per compound adapt FC's below 0 (artefacts due to natural abundance
#   #correction) to be positive to avoid problems with the visualisations
#   #later on.
#   for (i in (2:ncol(frac_tb))) {
#     if (any(frac_tb[,i]<0)) {
#       FCs<-pull(frac_tb[,i])
#       FCs[which(FCs<0)]<-FCs[which(FCs<0)]-min(FCs[which(FCs<0)]) 
#       frac_tb[,i]<-FCs    
#     }
#   }
#   
#   #modify iso_tb if it exists before summarizing
#   if (length(iso_tb)>0) {
#     iso_tb<-iso_tb %>%
#       pivot_wider(names_from = any_of(sample_col),values_from = value)
#       
#     #Per compound adapt isotopologues's below 0 (artefacts due to natural abundance 
#     #correction) to be positive to avoid problems with the visualisations
#     #later on
#     for (i in (2:nrow(iso_tb))) {
#       if (any(iso_tb[i,]<0)) {
#         #check if any value for this isotopologue below 0
#         metabolite<-iso_tb$Metabolite[i]
#         isos<-iso_tb[i,-c(1,2)]
#         negisos<-which(isos<0)
#         
#         #if no values negative, skip this section to avoid empty reference  
#         #warnings and useless computing. If negatives, no zero correction was done
#         #before and should be done now
#         if (length(negisos)>0) {
#           #Make variable containing negative iso value and 0 for others
#           #then overwrite negative iso values to 0
#           toadd<-isos
#           toadd[-negisos]<-0
#           isos[negisos]<-0   
#           iso_tb[i,-c(1,2)]<-isos
#           
#           #add negative iso values to parent to offset previous addition to 
#           #parent to compensate negative values
#           parent_index<-which(iso_tb$Metabolite==metabolite & 
#                                 iso_tb$Isotopologue==0)
#           iso_tb[parent_index,-c(1,2)]<-iso_tb[parent_index,-c(1,2)]+toadd
#           
#           #if any parents became <0, set to 0 (likely parent was undetectable)
#           iso_tb[parent_index,][which(iso_tb[parent_index,]<0&
#                                         is.numeric(iso_tb[parent_index,]))]<-0
#         }
#         
#       }
#     }
#   }
#   
#   #rename sample column in all inputs
#   meta_tb<-rename(meta_tb,Sample=all_of(sample_col))
#   abund_tb<-rename(abund_tb,Sample=all_of(sample_col))
#   frac_tb<-rename(frac_tb,Sample=all_of(sample_col))
#   
#   #add metadata to abundance and fractional contribution data respectively
#   #retaining only selected samples, and drop metabolites with 0 abundance
#   #in every sample to avoid errors
#   abund_tb<-left_join(meta_tb,abund_tb,by="Sample") %>%
#     select(1:ncol(meta_tb),any_of(compounds)) %>%
#     select_if(has_nonzero)
#   
#   frac_tb<-left_join(meta_tb,frac_tb,by="Sample") %>%
#     select(1:ncol(meta_tb),any_of(colnames(abund_tb))) 
#   
#   if(!length(iso_tb)==0) {
#     iso_tb<-iso_tb%>%
#       pivot_longer(3:ncol(.),names_to = "Sample",values_to = "value")%>%
#       pivot_wider(names_from = Metabolite,values_from = value)%>%
#       left_join(meta_tb,by="Sample") %>%
#       mutate(datatype=paste0("Isotopologue_",as.character(Isotopologue)),
#              across(any_of(colnames(abund_tb)),as.character)) %>%
#       select(any_of(colnames(meta_tb)),datatype,any_of(colnames(abund_tb)))%>%
#       filter(Sample %in% meta_tb$Sample)
#   }
#   
#   #add fractional contribution and isotopologues equal to 100% unlabeled to 
#   #compounds in abundance but not fraction labeling table
#   if (any(!colnames(abund_tb) %in% colnames(frac_tb))) {
#     nolabnames<-colnames(abund_tb)[which(! colnames(abund_tb) %in%
#                                            colnames(frac_tb))]
#     for (i in nolabnames) {
#       frac_tb$new<-0
#       colnames(frac_tb)[ncol(frac_tb)]<-i
#     }
#     if(!length(iso_tb)==0) {
#       for (i in nolabnames) {
#         iso_tb$new<-"1"
#         colnames(iso_tb)[ncol(iso_tb)]<-i
#       }
#     }
#   }
#   
#   
#   #prepare abundance data for joining: 
#   #calculate normalized abundances if normalization column provided and add
#   #to abund tb as different datatype. 
#   #add as character as isotopologue summaries will be character too
#   abund_tb <-abund_tb %>% add_column(datatype="Abund")
#   
#   if ("Normalisation" %in% colnames(meta_tb)) {
#     abund_tb<-abund_tb %>% 
#       mutate(across((ncol(meta_tb)+1):(ncol(abund_tb)-1),
#                     function(x) x/Normalisation)) %>%
#       mutate(datatype="NormAbund") %>%
#       full_join(abund_tb,by=colnames(abund_tb)) %>%
#       mutate(across(any_of(compounds),as.character)) 
#   } else {
#     abund_tb<-abund_tb %>%mutate(across(any_of(compounds),as.character)) 
#   }
#   
#   #prepare labeling  data for joining: 
#   #Add isotopologue data to fractional contribution data
#   frac_tb <-frac_tb %>% mutate(across(any_of(compounds),as.character)) %>%
#     add_column(datatype="FracCont")
#   
#   if(!length(iso_tb)==0) {
#     # iso_tb$datatype<-"Isotopologues"
#     frac_tb<-full_join(frac_tb,iso_tb,by=colnames(frac_tb))
#   }
#   
#   #join all tables then order and put in long format
#   #remove normalisation factor if present
#   tb<-full_join(frac_tb,abund_tb,by=colnames(abund_tb)) %>%
#     select(colnames(meta_tb),datatype,everything())%>%
#     pivot_longer(-c(any_of(colnames(meta_tb)),datatype),names_to = "compound",
#                  values_to = "value")%>%
#     na.omit()%>%
#     mutate(value=as.numeric(value))
#     
#     if (any("Normalisation" %in% colnames(meta_tb))) {
#       tb<-select(tb,-Normalisation) 
#     }
#   
#   return(tb)
# }  

#make list of factor orders per factor in data to allow multiple factors
#apply preselection of factor levels if desired, and order as supplied
clean_order_factors<-function(tb,factor_columns,factor_order=NULL) {
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  
  #if factor order is supplied as a vector, assume it is only the order
  #of the first factor
  if(length(factor_order)>0 & !is.list(factor_order)){
    factor_order<-list(factor_order)
  }
    
  if(length(factor_symbols)>0) {
    for(i in 1:length(factor_symbols)) {
      #if factor order is not given, just assume levels and order already
      #in data
      if(length(factor_order)<i){
          factor_order[[i]]<-unique(pull(tb,!!factor_symbols[[i]]))
      }
      # if(length(factor_order[[i]])==0){
      # } 
      # select only given factor levels
      tb<-tb %>%
        filter(!!factor_symbols[[i]] %in% factor_order[[i]]) %>%

        #Just in case factor is still in a character variable,
        #Change factor variable from text into actual factor for visualisation and
        #significance testing, then arrange data order to match the factor levels
        mutate(!!factor_symbols[[i]]:=factor(!!factor_symbols[[i]],
                                             levels = factor_order[[i]]))%>%
        arrange(!!factor_symbols[[i]])
    }
  }
  
  #return result after dropping unused factor levels
  return(tb%>%
           droplevels())
}

#join all data, keeping meta for last, then  drop unused factor
# levels and reorder them like input if specified(like those of blanks!)
join_metabo_longdata<-function(abund_longtb,frac_longtb,iso_longtb=NULL,meta_tb,
                    sample_column){
  tb<-full_join(abund_longtb,frac_longtb)%>%
    {
      if(length(iso_longtb)>0){
        full_join(.,iso_longtb)
      } else .
    }%>%
    left_join(meta_tb)%>%
    select(!!sample_symbol,any_of(colnames(meta_tb)),everything())
}

#Function that takes input specifying polly/dolly output data location and type, 
#and which metadata variables have to be taken along 
dolly_to_longtibble<-function(path,inputpath,metastring="meta",
                              abundstring="abund",
                              labelstring="iso",
                              isostring="iso",
                              minfract_detected=0,
                              sample_column,factor_column=NULL,
                              comparative_factor_column=NULL,
                              factor_levels_ordered=NULL,
                              tracer_column=NULL,
                              norm_column=NULL,
                              sampletype_column=NULL,
                              lib_tb=NULL){
  
  #make symbols for dplyr pipelines
  factor_columns <- c(factor_column,comparative_factor_column)
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  factor_symbol<-factor_symbols[[1]]
  if(length(factor_symbols)==2) compar_factor_symbol<-factor_symbols[[2]] else {
    compar_factor_symbol<-NULL
  }
  sample_symbol<-sym_or_null(sample_column,allownull = F)
  tracer_symbol<-sym_or_null(tracer_column)
  norm_symbol<-sym_or_null(norm_column)
  sampletype_symbol<-sym_or_null(sampletype_column)
  
  input_list<-
    list_inputdata_tbs(inputpath,metastring = metastring,
                       abundstring = abundstring,labelstring = labelstring,
                       isostring = isostring,sample_column = sample_column,
                       factor_columns = factor_columns,
                       norm_column = norm_column,
                       tracer_column = tracer_column,
                       sampletype_column = sampletype_column,lib_tb = lib_tb)
  
  #if no fractional contribution data present, 
  #calculate from isotopologue data
  if(!"frac_tb"%in% names(input_list)) {
    #calculate fractional contribution
    frac_worktb<-extract_col_isotopologues(input_list$iso_tb,
                                           iso_suffix_sep = "_")%>%
      select(-datatype)%>%
      calculate_FC()
    
  } else frac_worktb<-input_list$frac_tb
  
  #do checks on input data
  #generate error or warning messages if any
  check_output<-check_samples_compounds(
    meta_tb = input_list$meta_tb,
    abund_tb = input_list$abund_tb,
    frac_tb = frac_worktb,
    sample_column = sample_column,
    norm_column = norm_column)
  
  if (check_output$error) {
    validate(check_output$message)
  } else {
    outputtext<-check_output$message
  }
  
  #curate abundance data, save LOD data if present and note which compounds
  #are detected in less samples than required
  #detected too little
  abund_worktb<-input_list$abund_tb%>%
    curate_abundancedata(meta_tb=input_list$meta_tb,sample_column = sample_column,
                         norm_column = norm_column,
                         sampletype_column = sampletype_column)%>%
    select(-any_of(c(sampletype_column,norm_column)))
  
  compounds_toomany_undetected<-character(0)
  if("detected" %in% colnames(abund_worktb)){
    abund_LODtb<-abund_worktb%>%
      select(!!sample_symbol,compound,detected,any_of(c("LOD","LOD_blankcor")))%>%
      group_by(compound) %>%
      summarise(detected_fraction=length(which(detected))/n())%>%
      left_join(
        abund_worktb%>%
          select(compound,any_of(c("LOD","LOD_blankcor")))%>%
          unique(),
        by="compound")
    
    compounds_toomany_undetected<-abund_LODtb%>%
      filter(detected_fraction<minfract_detected)%>%
      pull(compound)
    
    write_csv(abund_LODtb,paste0(path,"/LOD and compound detection table.csv"))
  }
  
  if(length(compounds_toomany_undetected)>0){
    print(paste0("Following compounds are below LOD in more than ",
                 (1-minfract_detected)*100,"% of the samples."))
  }
  #transform abundance data to long format, saving all found types of abundance
  #with the common name
  abund_longtb<-abund_worktb%>%
    filter(!compound %in% compounds_toomany_undetected)%>%
    select(!!sample_symbol,
           compound,any_of(c("Abund","BlankcorAbund","NormAbund")))%>%
    pivot_longer(any_of(c("Abund","BlankcorAbund","NormAbund")),
                 names_to = "datatype",values_to = "value")
  
  #transform isotopologue and fraccon data to long format, saving all found types
  #of abundance with the common name, and dropping any samples or compounds not
  #in abund_longtb
  if("iso_tb"%in% names(input_list)) {
    iso_longtb<-extract_col_isotopologues(input_list$iso_tb,
                                          iso_suffix_sep = "_")%>%
      select(-Isotopologue)%>%
      filter(compound %in% abund_longtb$compound,
             !!sample_symbol %in% pull(abund_longtb,sample_column))
  } else iso_longtb<-NULL
  
  frac_longtb<-frac_worktb %>% 
    pivot_longer(2:ncol(.),names_to = "compound",values_to = "value")%>%
    mutate(datatype="FracCont",.before = 3)%>%
    filter(compound %in% abund_longtb$compound,
           !!sample_symbol %in% pull(abund_longtb,sample_column))
  
  
  #transform metatb for joining to long tibble
  input_list$meta_tb%>%
    select(-any_of(c(sampletype_column,norm_column)))
  
  
  meta_tb<-input_list$meta_tb%>%
    {
      if(length(sampletype_column)>0){
        if(sampletype_column %in% colnames(.)) {
          filter(.,!!sampletype_symbol != "blank")
        } else .
      } else .
    }%>%
    select(-any_of(c(sampletype_column,norm_column)))
  
  #join all data, keeping meta for last, then  drop unused factor
  # levels and reorder them like input if specified(like those of blanks!)
  #todo turn joinin series into function
  tb<-join_metabo_longdata(abund_longtb,frac_longtb,iso_longtb,meta_tb,
                           sample_column = sample_column)%>%
    clean_order_factors(factor_columns,factor_levels_ordered)
}
# Functions for generating pie charts -------------------------------------


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
  kwtest_results<-tibble(!!factor_symbol:=factor_order[[1]][-1],
                         p.value=NA)
  for(i in 2:length(factor_order[[1]])){
    tgt_cohort<-factor_order[[1]][i]
    partdata<-data %>% filter(!!factor_symbol %in% c(ref_cohort,tgt_cohort))
    if(length(unique(pull(partdata,factor_column)))<2) next
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
    factor_order[[1]]<-unique(pull(tb,!!factor_symbols[[1]]))
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
  
  #make summary tibble of fractional contribution and isotopologue data,
  #then add abundance data
  sum_tb<- prepare_tb %>%
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
    rowwise()%>%
    mutate(mod=list(kruskal_piedata(
      data=data,
      test_formula= as.formula(paste0("value ~ ",factor_column)),
      factor_column = factor_column,
      factor_order = factor_order)))%>%
    reframe(mod)%>%
    mutate(p.value=if_else(is.nan(p.value),
                           NA,
                           p.value))%>%
    right_join(
      prepare_tb %>%
        filter(grepl("frac",tolower(datatype))|
                 grepl("iso",tolower(datatype))) %>%
        group_by(.,compound,!!!factor_symbols,!!tracer_symbol,datatype)%>%
        summarise(average = mean(value))
    )%>%
    full_join(ab_sum_tb)%>%
    rename(P=p.value)
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
  
  
  if(length(tracers)>1){
    FC_position="slice"
  }
  
  rowwise(sum_tb) %>%    #to apply following functions per row  
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
           labFC=if_else(Abund==0,"ND",labFC)
    )%>%
    group_by(!!!rlang::syms(factor_columns)) %>%
    
    #get labeling positions on FC and abundance axes. Depends if in
    #center or in slice. Center if not detected (label ND)
    #If slice, set posFC as sum of current and all 
    #preceding FC's-half the current FC. Set posAb in slice at min_lab_dist radius 
    #if abundance smaller than twice min_lab_dist. 
    mutate(FClab_posAngle=if_else(FC_position=="center"|Abund==0,0,
                                  -cumsum(!!fraction_symbol)+!!fraction_symbol/2),
  FClab_posDist=if_else(FC_position=="center"|Abund==0|!!fraction_symbol==1,
                                 as.double(0),
                                 if_else(Abund<min_lab_dist*2,
                                         as.double(min_lab_dist),
                                         as.double(Abund/2))))%>%
    ungroup()                   #undo grouping
}

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
make_FC_slices<-function(sum_tb,factor_columns,tracer_column){
  
  #get symbols factor and tracer; and each tracer nutrient used
  #prepare factor and tracer symbols. Need a placeholder at the start if no
  #tracer information was provided
  factor_symbols<-sym_or_null(factor_columns,returnlist = T)
  if(length(tracer_column)==0){
    tracer_column<-"Labeling"
    sum_tb<-sum_tb %>%
      mutate(Labeling="Labeled")
  }
  tracer_symbol<-sym_or_null(tracer_column)
  nutrient_symbols<-sum_tb%>%
    filter(grepl("frac",tolower(datatype),fixed = T)) %>%
    pull(tracer_column)%>%
    unique()%>%
    rlang::syms()

  
  #take fractional contribution data only, keep only slices above 0
  # join abundances or if desired
  #normalized abundances to it, calculate the fractions as abundance*fraccon
  #for each tracer, and assume a single unnamed tracer when no
  #tracer column is provided. Rejoin the fraccon p values and add labels for
  #for the plots
  sum_tb%>%
    filter(grepl("frac",tolower(datatype),fixed = T),
           average>0) %>%
    rename(FracCont=average,P_FC=P)%>%
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
    ) %>%
    group_by(compound)%>%
    mutate(Abund=Abund/max(Abund),
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
    mutate(FracCont=Fraction/Abund)%>%
    left_join(
      sum_tb%>%
        rename(P_FC=P) %>%
        select(P_FC,compound,datatype,!!tracer_symbol,!!!factor_symbols)
    )%>%
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
    ungroup()
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
  
  #turn tracer column into factor, save original levels in order of appearance
  #extract data of selected compound only, keep only fractions above 0
  #keep only tracer colors linked to existing levels
  slice_mod_tb <- slice_tb %>%
    mutate(!!tracer_symbol:=factor(!!tracer_symbol,
                                   levels=unique(!!tracer_symbol))
           )
  orig_levels<-levels(slice_mod_tb$Tracer)
    
  slice_mod_tb <- slice_mod_tb %>%
    filter(compound==selected_compound,
           Fraction>0)
  col_labeling<-col_labeling[which(orig_levels %in% 
                                     as.character(unique(slice_mod_tb$Tracer)))]
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
    
    if (length(unique(pull(slice_mod_tb[,tracer_column])))>2) {
      plotrect<-plotrect  +
        geom_text_repel(data=slice_mod_tb,
                        aes(label=paste0(labFC,"\n",P_FClab)),
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
  
  plotrect
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

make_slicetb<-function(tb,detail_charts,pathway_charts,savepath,
                       normalize=T,factor_columns,tracer_column,factor_order,label_decimals,
                       percent_add,FC_position,min_lab_dist,P_isotopologues,
                       log_abund,circlelinecolor,circlelinetypes,
                       maxcol_facet,include_name,col_labeling,
                       alpha,otherfontsize,
                       font,legendtitlesize,cohortsize,include_legend,
                       mapotherfontsize=16,mapcohortsize=18,format="png",
                       show_P=T) {
  
  #Add dummy tracer column called labeling,
  #if tracer column is missing from dataframe
  if(!tracer_column %in% colnames(tb)) {
    tracer_symbol<-rlang::sym(tracer_column)
    tb<-tb %>%mutate(!!tracer_symbol:="Labeled")
  }
  
  print(paste0("Processing ",compound))
  
  #prepare filename, remove problematic characters
  if (normalize) {
    plotfilename<-paste0("pies normalized ",compound,".",format)
  } else {
    plotfilename<-paste0("pies ",compound,".",format)
  }
  plotfilename<-gsub("/","-",plotfilename)
  
  #get table with only measured compound data, then a table summarizing
  #derived means and p values per cohort for abundance and one for fractional
  #contribution, then put together table with inputformat for pie function
  print(paste0("extracting compounddata"))
  
  compound_tb<-obtain_compounddata(tb,compound,factor_columns = factor_columns,
                                   tracer_column = tracer_column,
                                   factor_order = factor_order,
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
  
  slice_tb<-prepare_slicedata(compound_tb,factor_columns = factor_columns,
                              tracer_column=tracer_column,
                              compound=compound,label_decimals = label_decimals,
                              min_lab_dist = min_lab_dist,
                              percent_add = percent_add,
                              FC_position = FC_position,
                              P_isotopologues=P_isotopologues)
  
  if (detail_charts) {
    #plot detailed chart based on information in slice table
    print(paste0("saving detailed chart"))
    
    pies<-make_piechart(slice_tb,factor_columns = factor_columns,
                        tracer_column = tracer_column,
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
    pies<-make_piechart(slice_tb,factor_columns = factor_columns,
                        tracer_column=tracer_column,
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

generate_pie_old<-function(tb,compound,detail_charts,pathway_charts,savepath,
                       normalize=T,factor_columns,tracer_column,factor_order,label_decimals,
                       percent_add,FC_position,min_lab_dist,P_isotopologues,
                       log_abund,circlelinecolor,circlelinetypes,
                       maxcol_facet,include_name,col_labeling,
                       alpha,otherfontsize,
                       font,legendtitlesize,cohortsize,include_legend,
                       mapotherfontsize=16,mapcohortsize=18,format="png",
                       show_P=T) {
  
  #Add dummy tracer column called labeling,
  #if tracer column is missing from dataframe
  if(!tracer_column %in% colnames(tb)) {
    tracer_symbol<-rlang::sym(tracer_column)
    tb<-tb %>%mutate(!!tracer_symbol:="Labeled")
  }
  
  print(paste0("Processing ",compound))
  
    #prepare filename, remove problematic characters
  if (normalize) {
    plotfilename<-paste0("pies normalized ",compound,".",format)
  } else {
    plotfilename<-paste0("pies ",compound,".",format)
  }
  plotfilename<-gsub("/","-",plotfilename)
  
  #get table with only measured compound data, then a table summarizing
  #derived means and p values per cohort for abundance and one for fractional
  #contribution, then put together table with inputformat for pie function
  print(paste0("extracting compounddata"))
  
  compound_tb<-obtain_compounddata(tb,compound,factor_columns = factor_columns,
                                     tracer_column = tracer_column,
                                   factor_order = factor_order,
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
  
  slice_tb<-prepare_slicedata(compound_tb,factor_columns = factor_columns,
                                tracer_column=tracer_column,
                              compound=compound,label_decimals = label_decimals,
                              min_lab_dist = min_lab_dist,
                              percent_add = percent_add,
                              FC_position = FC_position,
                              P_isotopologues=P_isotopologues)
  
  if (detail_charts) {
    #plot detailed chart based on information in slice table
    print(paste0("saving detailed chart"))
    
    pies<-make_piechart(slice_tb,factor_columns = factor_columns,
                            tracer_column = tracer_column,
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
    pies<-make_piechart(slice_tb,factor_columns = factor_columns,
                            tracer_column=tracer_column,
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
           normalize=T,factor_columns,tracer_column,factor_order,label_decimals,percent_add,
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
                 normalize=normalize,factor_columns=factor_columns,
                 tracer_column=tracer_column,factor_order=factor_order,
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


# #test isotopologue adapated functions
# example_tb<-read_csv(
#   file = "~/GitHub/mec-shiny-apps/shiny-server/TraVisPies/Example_data/Standardized input/Input_Example_standardized w isotopologues.csv")
# inputtb<-example_tb[,c(2:4)] %>%
#   filter(!datatype=="Abund")
# head(inputtb)
# 
# # Coenzyme_A
# v_settings<-list(compound="Coenzyme_A",
#                  factor_columns=colnames(inputtb)[1],
#                  factor_order=pull(unique(inputtb[,1])),
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
#   example_tb,compound=v_settings$compound,factor_columns = v_settings$factor_columns,
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
# slicetb<-prepare_slicedata(comptb,factor_columns = v_settings$factor_columns,
#                       compound=v_settings$compound,
#                       label_decimals = v_settings$label_decimals,
#                       min_lab_dist = v_settings$min_lab_dist,
#                       percent_add = v_settings$percent_add,
#                       FC_position = v_settings$FC_position,
#                       P_isotopologues=v_settings$P_isotopologues)
# 
# #makes pie chart based on table with required data per pie slice
# (test<-make_piechart(slicetb,factor_columns = v_settings$factor_columns,
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
#              normalize=v_settings$norm,factor_columns=v_settings$factor_columns,
#              factor_order=v_settings$factor_order,
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
#                        normalize=v_settings$norm,factor_columns=v_settings$factor_columns,
#                        factor_order=v_settings$factor_order,
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
# print(create_caption(factor_order = v_settings$factor_order,
#                log_abund = v_settings$log_abund,
#                FC_position = v_settings$FC_position,
#                show_P = v_settings$show_P,
#                P_isotopologues = v_settings$P_isotopologue))

#test pie chart function for logscale abundance ratios
# factor_columns<-"Cohort"
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
#     geom_text(aes(label=P_RAlab),x=1.6,y=7/8,size=otherfontsize*5/14,
#               hjust="inward",vjust="inward") +
#     geom_text(aes(label=P_FClab),x=1.6,y=5/8,size=otherfontsize*5/14,
#               hjust="inward",vjust="inward")
# }
# 
# #transform bar to pie chart and plot pies on grid.
# piebasic<-plotrect+
#   facet_wrap(vars(!!rlang::sym(factor_columns)),ncol=maxcol_facet) +
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






# old discontinued functions ----------------------------------------------

# #Summarize extracted isotopologue data as a single row per sample 
# #containing for each metabolite a string with all contributions from lowest to 
# #highest isotopologue in order separated by |
# summarize_isotopologue<-function(iso_tb,sample_colname="Sample",){
#   #remove isotopologue column, group per metabolite and generate single string
#   #per metabolite for each sample
#   #then transpose from rowwise to columnwise representation
#   iso_tb %>% select(-Isotopologue) %>%
#     group_by(Metabolite) %>%
#     summarize(across(everything(),.fns = ~ paste0(.x,collapse = "|"))) %>%
#     t_tibble(first_colname = sample_colname)
# }
# 
# merge_input_allow_summarize<-function(meta_tb,abund_tb,frac_tb,iso_tb=NULL,
#                                       sample_col="Sample",compounds,summarize_isos=F) {
#   #Per compound adapt FC's below 0 (artefacts due to natural abundance
#   #correction) to be positive to avoid problems with the visualisations
#   #later on.
#   for (i in (2:ncol(frac_tb))) {
#     if (any(frac_tb[,i]<0)) {
#       FCs<-pull(frac_tb[,i])
#       FCs[which(FCs<0)]<-FCs[which(FCs<0)]-min(FCs[which(FCs<0)]) 
#       frac_tb[,i]<-FCs    
#     }
#   }
#   
#   #modify iso_tb if it exists before summarizing
#   if (length(iso_tb)>0) {
#     #Per compound adapt isotopologues's below 0 (artefacts due to natural abundance 
#     #correction) to be positive to avoid problems with the visualisations
#     #later on
#     for (i in (2:nrow(iso_tb))) {
#       if (any(iso_tb[i,]<0)) {
#         #check if any value for this isotopologue below 0
#         metabolite<-iso_tb$Metabolite[i]
#         isos<-iso_tb[i,-c(1,2)]
#         negisos<-which(isos<0)
#         
#         #if no values negative, skip this section to avoid empty reference  
#         #warnings and useless computing. If negatives, no zero correction was done
#         #before and should be done now
#         if (length(negisos)>0) {
#           #Make variable containing negative iso value and 0 for others
#           #then overwrite negative iso values to 0
#           toadd<-isos
#           toadd[-negisos]<-0
#           isos[negisos]<-0   
#           iso_tb[i,-c(1,2)]<-isos
#           
#           #add negative iso values to parent to offset previous addition to 
#           #parent to compensate negative values
#           parent_index<-which(iso_tb$Metabolite==metabolite & 
#                                 iso_tb$Isotopologue==0)
#           iso_tb[parent_index,-c(1,2)]<-iso_tb[parent_index,-c(1,2)]+toadd
#           
#           #if any parents became <0, set to 0 (likely parent was undetectable)
#           iso_tb[parent_index,][which(iso_tb[parent_index,]<0&
#                                         is.numeric(iso_tb[parent_index,]))]<-0
#         }
#         
#       }
#     }
#     
#     #if desired (needed in travis pies) summarize isotopologue data with name
#     #sample column, otherwise make row per sample isotopologue combo with
#     #isotopologue_isotopologueNR as datatype
#     if (summarize_isos) {
#       iso_tb<-summarize_isotopologue(iso_tb,sample_colname = sample_col)
#     } else {
#       iso_tb<-iso_tb %>% 
#         mutate(datatype=paste0("Isotopologue_",as.character(Isotopologue))) %>%
#         select(Metabolite,datatype,everything(),-Isotopologue) %>%
#         pivot_longer(3:ncol(.),names_to = "Sample",values_to = "value") %>%
#         pivot_wider(names_from = Metabolite,values_from = value)%>%
#         select(Sample,everything())%>%
#         mutate(across(any_of(compounds),as.character))
#     }
#   }
#   
#   #rename sample column in all inputs
#   meta_tb<-rename(meta_tb,Sample=all_of(sample_col))
#   abund_tb<-rename(abund_tb,Sample=all_of(sample_col))
#   frac_tb<-rename(frac_tb,Sample=all_of(sample_col))
#   
#   #add metadata to abundance and fractional contribution data respectively
#   #retaining only selected samples, and drop metabolites with 0 abundance
#   #in every sample to avoid errors
#   abund_tb<-left_join(meta_tb,abund_tb,by="Sample") %>%
#     select(1:ncol(meta_tb),any_of(compounds)) %>%
#     select_if(has_nonzero)
#   
#   frac_tb<-left_join(meta_tb,frac_tb,by="Sample") %>%
#     select(1:ncol(meta_tb),any_of(colnames(abund_tb))) 
#   
#   if(!length(iso_tb)==0) {
#     iso_tb<-left_join(iso_tb,meta_tb,by="Sample") %>%
#       select(any_of(colnames(meta_tb)),datatype,any_of(colnames(abund_tb)))%>%
#       filter(Sample %in% meta_tb$Sample)
#   }
#   
#   #add fractional contribution and isotopologues equal to 100% unlabeled to 
#   #compounds in abundance but not fraction labeling table
#   if (any(!colnames(abund_tb) %in% colnames(frac_tb))) {
#     nolabnames<-colnames(abund_tb)[which(! colnames(abund_tb) %in%
#                                            colnames(frac_tb))]
#     for (i in nolabnames) {
#       frac_tb$new<-0
#       colnames(frac_tb)[ncol(frac_tb)]<-i
#     }
#     if(!length(iso_tb)==0) {
#       for (i in nolabnames) {
#         iso_tb$new<-"1"
#         colnames(iso_tb)[ncol(iso_tb)]<-i
#       }
#     }
#   }
#   
#   
#   #prepare abundance data for joining: 
#   #calculate normalized abundances if normalization column provided and add
#   #to abund tb as different datatype. 
#   #add as character as isotopologue summaries will be character too
#   abund_tb <-abund_tb %>% add_column(datatype="Abund")
#   
#   if ("Normalisation" %in% colnames(meta_tb)) {
#     abund_tb<-abund_tb %>% 
#       mutate(across((ncol(meta_tb)+1):(ncol(abund_tb)-1),
#                     function(x) x/Normalisation)) %>%
#       mutate(datatype="NormAbund") %>%
#       full_join(abund_tb,by=colnames(abund_tb)) %>%
#       mutate(across(any_of(compounds),as.character)) 
#   } else {
#     abund_tb<-abund_tb %>%mutate(across(any_of(compounds),as.character)) 
#   }
#   
#   #prepare labeling  data for joining: 
#   #Add isotopologue data to fractional contribution data
#   frac_tb <-frac_tb %>% mutate(across(any_of(compounds),as.character)) %>%
#     add_column(datatype="FracCont")
#   
#   if(!length(iso_tb)==0) {
#     # iso_tb$datatype<-"Isotopologues"
#     frac_tb<-full_join(frac_tb,iso_tb,by=colnames(frac_tb))
#   }
#   
#   #join all tables then order and put in long format
#   #remove normalisation factor if present
#   tb<-full_join(frac_tb,abund_tb,by=colnames(abund_tb)) %>%
#     select(colnames(meta_tb),datatype,everything())%>%
#     pivot_longer(-c(any_of(colnames(meta_tb)),datatype),names_to = "compound",
#                  values_to = "value")%>%
#     na_omit()%>%
#     
#     if ("Normalisation" %in% colnames(meta_tb)) {
#       tb<-select(tb,-Normalisation) 
#     }
#   
#   return(tb)
# }
#test input merging -----------------------------------------------------------------------
# meta_tb<-read_csv_clean(paste0(getwd(),
#                                "/Example_data/Original input/Input_Example_metadata.csv"),
#                         remove_empty = T,perc_to_num = F)
# # meta_tb[,3]<-1
# abund_tb<-read_csv_clean(paste0(getwd(),
#                                 "/Example_data/Original input/Input_Example_RA.csv"),
#                          remove_empty = T,perc_to_num = F)
# iso_et_tb<-read_csv_clean(paste0(getwd(),
#                                  "/Example_data/Original input/Input_Example_RA.csv"),
#                           remove_empty = F,perc_to_num = F)
# iso_col_tb<-read_csv_clean(paste0(getwd(),
#                                   "/Example_data/Original input/Input_Example_isotopologues.csv"),
#                            remove_empty = T,perc_to_num = T)
# iso_tb<-extract_col_isotopologues(iso_col_tb) %>%
#     slice(-c(5,6,7,8))
# 
# # iso_tb<-extract_et_isotopologues(iso_et_tb) %>%
# #   slice(-c(5,6,7,8))
# 
# # abund_tb<-extract_et_abund(iso_et_tb,sample_colname = "Sample")
# 
# frac_tb<-calculate_FC(iso_tb,sample_colname = "Sample")
# 
# sample_col<-"Sample"
# 
# compounds<-colnames(abund_tb)[-1]
# head(meta_tb)
# (meta_formatted_tb<-format_metadata(meta_tb,sample_column = "Sample",
#                                    factor_columns = "Cohort",
#                                    norm_column = "None"))
# 
# test<-merge_input(meta_tb = meta_formatted_tb,
#                   abund_tb = abund_tb,
#                   frac_tb = frac_tb,
#                   iso_tb=iso_tb,
#                   sample_col = sample_col,
#                   compounds = compounds)
# 
# abund_tb<-abund_tb %>%mutate(across(-c(1:3),as.character)) 
# 
# 
# 
# 
