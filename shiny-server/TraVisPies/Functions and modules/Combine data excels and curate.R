# Libraries and functions -------------------------------------------------
library(here)    #to make r source from file location instead of magic stuff
library(tibble)       #for manipulating tibbles
library(tidyr)        #for restructuring data tibbles to tidy format (eg pivot)
library(readxl)
library(dplyr)
library(readr)        #for writing .csv file of merged output
library(ggplot2)      #for generating the pie chart plots


# function for checking if any column cell contains non-NA data
has_data <- function(x) { sum(!is.na(x)) > 0 } 

# function for loading and cleaning abundance and FC sheets, add file as variable
read_excelsheet_clean<- function(file,path,sheet=NULL,datatype,remove_empty=FALSE,perc_to_num=T,
                          remove_rowempty=FALSE){
  #read in file in path
  filepath<-paste0(path,"/",file)
  totalcol<-ncol(read_excel(path = filepath,sheet=sheet))
  coltypes<-c("text",rep("numeric",times=totalcol-1))
  input_tb<-read_excel(path = filepath,sheet=sheet,col_types = coltypes)
  
  #drop empty columns and rows if desired
  if (remove_empty) {
    input_tb<-select_if(input_tb,has_data)          
  }
  
  if (remove_rowempty) {
    input_tb<-input_tb %>% na.omit()          
  }
  
  #set percentage strings to fractions if desired
  if (perc_to_num){
    percolumns<-grep("%",input_tb)
    
    
    input_tb<-mutate(input_tb,across(all_of(percolumns),function(x) 
      as.numeric(sub(pattern="%", replacement = "",x,fixed = T))/100))
  }
  
  #turn into long format, add filename as variable
  input_tb<-input_tb %>%
    pivot_longer(cols = 2:last_col(),names_to = "compound",
                                  values_to = "value") %>%
    mutate(file=file,datatype=datatype)%>%
    select(file,everything())
  
  return(input_tb)
}

# function for loading and cleaning metadata sheets, add file as variable
read_excelsheet_meta<-function(file,path,sheet=NULL,datatype,remove_empty=FALSE,perc_to_num=T,
                              remove_rowempty=FALSE){
  #read in file in path, add file name column
  filepath<-paste0(path,"/",file)
  input_tb<-read_excel(path = filepath,sheet=sheet)%>%
    mutate(file=file)
    
  
  #drop empty columns and rows if desired
  if (remove_empty) {
    input_tb<-select_if(input_tb,has_data)          
  }
  
  if (remove_rowempty) {
    input_tb<-input_tb %>% na.omit()          
  }
  
  return(input_tb)
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

#Standardize tibble
standardize_sampletb <- function(tb) { 
  if ("abundance" %in% tolower(colnames(tb))) {
    colnames(tb)[which(tolower(colnames(tb))=="abundance")]<-"Abundance"
  }
  
  #replace NA by 0, then long format with one value
  tb <- tb %>%
    replace(is.na(.), 0)
  return(tb)
} 

#Make standardized input table
join_standardize_folderexcels<-function(path,metasheetname,abundsheetname,
                                        FCsheetname) {
  #get existing xlsx and xls files in path
  excelfiles<-list.files(path=path,pattern=".xls")
  
  #read in desired sheets in excel files one by one
  #todo replace join_by to exclude references to abundance or FC
  for (i in 1:length(excelfiles)) {
    #initialise tibble in first file, add to it afterwards
    if(i ==1) {
      print(excelfiles[i])
      #get metadata
      meta_tb<-read_excelsheet_meta(excelfiles[i],path,metasheetname,
                                    remove_empty = TRUE,
                                    remove_rowempty = TRUE)
      
      #get and join abundance and fc data, drop entries not in metadata
      value_tb<-read_excelsheet_clean(excelfiles[i],path,abundsheetname,
                                      datatype = "Abundance",
                                      remove_empty = TRUE,
                                      remove_rowempty = TRUE) %>%
        full_join(
          read_excelsheet_clean(excelfiles[i],path,FCsheetname,
                                datatype = "FC",
                                remove_empty = TRUE,
                                remove_rowempty = TRUE),
          by = join_by(file, Sample_Names, compound,value,datatype)
        )  %>%
        filter(Sample_Names %in% meta_tb$Sample_Names)
      
    } else {
      print(excelfiles[i])

      #get metadata
      meta_tb<-meta_tb %>%
        full_join(
          read_excelsheet_meta(excelfiles[i],path,metasheetname,
                               remove_empty = TRUE,
                               remove_rowempty = TRUE),
          by = join_by(Sample_Names, sample_type, file)
        )
      #get and join abundance and fc data, drop entries not in metadata
      value_tb<-value_tb %>%
        full_join(
          read_excelsheet_clean(excelfiles[i],path,abundsheetname,
                                datatype = "Abundance",
                                remove_empty = TRUE,
                                remove_rowempty = TRUE),
          by = join_by(file, Sample_Names, compound,value,datatype)
        ) %>%
        full_join(
          read_excelsheet_clean(excelfiles[i],path,FCsheetname,
                                datatype = "FC",
                                remove_empty = TRUE,
                                remove_rowempty = TRUE),
          join_by(file, Sample_Names, compound,value,datatype)
        ) %>%
        filter(Sample_Names %in% meta_tb$Sample_Names)
      
    }
    
    if(i==length(excelfiles)) {
      print("Join values and metadata, standardizing file")
      all_tb<-value_tb %>%
        left_join(meta_tb,by = join_by(file, Sample_Names))%>%
        standardize_sampletb()
    }
  }
  
  all_tb<-rename(all_tb,sample_name=Sample_Names)
  return(all_tb)
}

#replace derivatized compound vector names by originals if compound found
#in library
rename_origcompound<-function(compoundcol_tb,lib_tb,
                              compoundcol="compound",
                              origcompoundcol="Orig_name"){
  compound_symbol<-rlang::sym(colnames(compoundcol_tb)[1])
  
  origcompoundcol_tb <- compoundcol_tb %>%
    mutate(
      !!compound_symbol:=if_else(
        !!compound_symbol %in% lib_tb[,compoundcol],
        lib_tb[which(lib_tb[,compoundcol]==compound),origcompoundcol],
        !!compound_symbol
      )
    )
  # if (compound %in% lib_tb[,compoundcol]) {
  #   origcompound<-lib_tb[which(lib_tb[,compoundcol]==compound),
  #                        origcompoundcol]
  # } else {
  #   origcompound<-compound
  # }
  
  return(origcompound)
}

#Curate joined data-------------------------------------------------
#read joined data
path<-here::here("input excel")
file<-"combined_excel_data.csv"
input_tb<-read_csv(paste0(path,"/",file))
libfile<-"lib_compound_deri_and_stdisation.csv"

#input for normalisation

normfiles<-c("MCF001613 sicrit.xlsx","MCF001712 sicrit.xlsx",NULL)
normcomps<-c("Methyl_myristate","Methyl_myristate",NULL)

#rename compounds if references to pH from hilic library are present
# then if library given that contains compound and orig compound
#column. Then remove _pH9, _pH11 or similar characters from compound name
#add columns needed later
curated_tb<-curatedsavetb<- input_tb %>%
  mutate(compound = sub("_pH\\d+","",compound),
         IS=NA,
         IS_abund=NA)

if (libfile %in% list.files(path=path,pattern=".csv")) {
  lib_tb<-read_csv(paste0(path,"/",libfile))
  if("compound" %in% colnames(lib_tb)& "Orig_name" %in% colnames(lib_tb)) {
    curated_tb <- curated_tb %>% 
      left_join(lib_tb[,c("compound","Orig_name")],
                by = join_by(compound)) %>%
      mutate(orig_compound=if_else(is.na(Orig_name),compound,Orig_name)) %>%
      select(-Orig_name)
  } else {
    warning("No compounds were renamed as specified library ",file,
            " does not contain the column compound and/or Orig_name")
  }
} else {
  warning("No compounds were renamed as specified library ",file,
        " was not found in ",path)
}

#add internal standard data if any, then perform correction
for (i in 1:length(normfiles)) {

  normabund_tb<-curated_tb %>% filter(file==normfiles[i],
                                      datatype=="Abundance")%>%
    left_join(curated_tb %>% filter(file==normfiles[i],
                                    compound==normcomps[i],
                                    datatype=="Abundance") %>% 
                mutate(IS_abund=value,IS=normcomps[i])%>%
                select(file,sample_name,IS,IS_abund),
              by=join_by(file,sample_name))%>%
    mutate(IS=coalesce(IS.x,IS.y),
           IS_abund=coalesce(IS_abund.x,IS_abund.y),
           .keep="unused")%>%
    group_by(compound)%>%
    mutate(datatype="NormAbundance",value=value/IS_abund*mean(IS_abund))%>%
    ungroup() 
  
  curated_tb <-curated_tb %>%
    left_join(normabund_tb %>% select(file,sample_name,IS,IS_abund) %>% 
                distinct(),
              by = join_by(file, sample_name)) %>%
    mutate(IS=coalesce(IS.x,IS.y),
           IS_abund=coalesce(IS_abund.x,IS_abund.y),
           .keep="unused")
}

curated_tb<-curated_tb %>% 
  mutate(IS_abund=if_else(is.na(IS_abund),1,IS_abund))%>%
  full_join(curated_tb %>% mutate(IS_abund=if_else(is.na(IS_abund),1,IS_abund))%>%
              filter(datatype=="Abundance") %>%
              mutate(datatype="NormAbundance",
                     value=value/IS_abund*mean(IS_abund))) 

write_csv(curated_tb,paste0(path,"/curated_excel_data.csv"))

#Make tibble joined excelfiles SICRIT comparison, and save as csv -------------------------------------------------
#one entry per valueper compound per sample per file

#read in all excel data in a folder. These sheets will be read
#metadata => all metadata (cohort, labelled, normalisation value)
path<-here::here("input excel")

write_csv(join_standardize_folderexcels(path,metasheetname = "metadata",
                                        abundsheetname = "correctedAbundances",
                                        FCsheetname = "fracContribution_C13"),
          paste0(path,"/combined_excel_data.csv"))

