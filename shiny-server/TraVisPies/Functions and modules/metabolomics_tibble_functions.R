# Description ---------------------------------------------------

###Author: Sam De Craemer
#Vlaams Instituut voor Biotechnologie (VIB) and KULeuven
#Metabolomics Expertise Center (MEC)

###Summary: Functions used to transform metabolomics data in different input 
#formats to the same type of long tibble fit for ggplot2.

#Libraries ---------------------------------------------------------------
#libraries for UI
library(dplyr)        #for faster.easier manipulation of data
library(tibble)       #for manipulating tibbles
library(vroom)        #for easier file loading
library(forcats)      #for factor manipulation
library(readr)        #for writing .csv file of merged output
library(tidyr)        #for restructuring data tibbles
library(stringr)      #for padding leading zeros to isotopologue strings
library(shiny)        #for error outputs to UI

# Functions ---------------------------------------------

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
  #confusion with isotopologues), but not renaming metadata columns like 
  #samplecolumn if present elsewhere. Remove empty columns and rows too.
  #rename compounds by library if provided
  output_list<-list(meta_tb=meta_tb%>%
                      format_metadata(sample_column = sample_column,
                                      factor_columns = factor_columns,
                                      norm_column = norm_column,
                                      tracer_column=tracer_column,
                                      sampletype_column=sampletype_column),
                    abund_tb=abund_tb%>%
                      rename_lib(lib_tb)%>%
                      rename_with(~gsub("_"," ",.x),
                                  .cols=-any_of(colnames(meta_tb)))%>%
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
      rename_with(replace_except_last,
                  .cols=-any_of(colnames(meta_tb)))
      
  } else {
    output_list$frac_tb<-label_tb %>%
      rename_with(~gsub("_"," ",.x),
                  .cols=-any_of(colnames(meta_tb)))
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
extract_col_isotopologues<-function(iso_col_tb,iso_suffix_sep="_",
                                    sample_column="Sample") {
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
    pivot_longer(4:ncol(.),names_to = sample_column,values_to = "value") %>%
    select(all_of(sample_column),everything())
}

#Extract abundance data in columns from Escher-Trace like corrected isotopologue
#file
extract_et_abund<-function(iso_et_tb,sample_column="Sample"){
  abund_tb<-filter(iso_et_tb,!is.na(Metabolite)) %>%
    t_tibble(first_colname = sample_column) %>%
    slice(-c(1)) %>%
    mutate(across(!(!!sample_column),.fns= as.numeric))
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
calculate_FC<-function(iso_tb,sample_column="Sample"){
  sample_symbol<-rlang::sym(sample_column)
  
  #calculate FC table, replacing levels lower than 0 with 0 as all isotopologues
  #must be kept, then reformat to columnwise format
  iso_tb %>%
    mutate(value=if_else(value<0,
                         0,
                         value)) %>%
    group_by(compound,!!sample_symbol) %>%
    summarise(value = sum(value*Isotopologue)/max(Isotopologue))%>%
    select(!!sample_symbol,everything())%>%
    pivot_wider(names_from = compound,values_from = value)
    
  # %>%
  #   t_tibble(first_colname = sample_column)
}

#Function to check samples across meta, abundance and FC tibbles and check
# compounds present. Outputs a list noting whether an error message should
# be given, and a message containing the error message or in absence
# of the error any warning messages to display
check_samples_compounds<-function(meta_tb,abund_tb,frac_tb,sample_column,
                                        norm_column){
  #set error =T by default, will change once past all error checks
  outlist<-list(error=T,message=NULL)
  
  #check if the sample column is present in all tables
  if(!sample_column %in% colnames(meta_tb) | 
     !sample_column %in% colnames(abund_tb)|
     !sample_column %in% colnames(frac_tb)) {
    outlist$message<-paste0("The sample column name `",sample_column,"` was ",
                            "not present in at least one of the inputdata ",
                            "files: metadata, abundance data or labeling data.")
    return(outlist)
  }
  
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
  if (length(norm_column)>0) {
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
  
  #return empty text if no warnings, else print (for use in non-UI apps) and 
  #give them in single orange text (html) for shiny apps
  if (length(outlist$message)>0) {
    print(outlist$message)
    outlist$message<-paste("<b><p style='color:orange'>Warning: </b>",
                           outlist$message,
                           "</p>", sep = "<br/>")
  } else {
    outlist$message<-""
  }
  return(outlist)
}

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
  sample_symbol<-sym_or_null(sample_column)
  
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
                              lib_tb=NULL,
                              savedata=T){
  
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
                                           iso_suffix_sep = "_",
                                           sample_column=sample_column)%>%
      select(-datatype)%>%
      calculate_FC(sample_column=sample_column)
    
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
  #in abund_longtb. 
  if("iso_tb"%in% names(input_list)) {
    iso_longtb<-extract_col_isotopologues(input_list$iso_tb,
                                          iso_suffix_sep = "_",
                                          sample_column = sample_column)%>%
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
  
  if(savedata) {
    write_csv(tb,paste0(path,"/merged data long tibble.csv"))
    saveRDS(tb,paste0(path,"/merged data long tibble.Rdata"))
  }
  
  
  return(tb)
}


