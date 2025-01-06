# Description --------------------------------------------------
#Written by Sam De Craemer

#Aim: convert a dolly excel with added metadata sheet to a long .csv file
#compatible with ggplot and dplyr for R data analysis and plotting.  

# Packages and functions --------------------------------------------------
library(dplyr)
library(readr)
library(readxl)
library(tidyr) #pivot function

#Chatgpt: replace all but the last occurrence of a character in elements of 
#a string vector with a different character
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



#reads in an abundance, fractional contribution or isotopologue sheet of a 
#dolly excel file by name,and prepares the desired table from it
read_dollysheet_to_long<-function(excelfile,sheetname,datatypename,lib_tb=NULL,meta_tb) {
  #read excel sheet data, if derivatized change names to underivatized
  excel_tb <- read_excel(excelfile,sheetname)
  # excel_tb <- read_excel(excelfile,"correctedIsotopologues_C13")
  
  if(any(colnames(lib_tb)=="Orig_name")) {
    for(i in 1:nrow(lib_tb)){
      colnames(excel_tb)<-sub(lib_tb$compound[i],
                              lib_tb$Orig_name[i],
                              colnames(excel_tb),
                              fixed = T)
    }  
  }
  
  #replace all but last _ in isotopologe compound names for easy 
  #substringing later. To match names, make sure to replace all _ in 
  #non-isotopologue compound names
  if (grepl("isotopologue",tolower(sheetname))) {
    colnames(excel_tb)<-replace_except_last(colnames(excel_tb))
  } else {
    colnames(excel_tb)<-gsub("_"," ",colnames(excel_tb))
  }

  
  #detect internal standards in sheet
  headers<-excel_tb %>% 
    select(where(~ all(is.na(.)))) %>%
    colnames()
  
  if (length(headers[which(grepl("internal",tolower(headers)))])>0) {
    intstdfirstcol<-which(colnames(excel_tb) ==
                            headers[which(grepl("internal",
                                                tolower(headers)))][1])+1
    intstdlastcol<-which(colnames(excel_tb) ==
                           headers[which(grepl("internal",
                                               tolower(headers)))+1][1])-1
    intstds<-colnames(excel_tb)[intstdfirstcol:intstdlastcol]
  } else {
    intstds<-NULL
  }  
  
  #remove empty rows and columns and check if abundance sheet for specific actions to 
  #take with it. Convert data to long format if not done before.
  excelclean <- excel_tb %>% 
    select(where(~ !all(is.na(.))))%>%
    na.omit()%>%
    rename(sample=1)
  
  long_format<-F
  if (grepl("abundance",tolower(sheetname))){
    
    #if LOD in sheet calculate LOD and average blank from supposed mock samples,
    #and remove these samples from the abundance sheet. Set blank and LOD to 0
    #for internal standards.
    #then remove rows that contain blanks, lod from sheet, subtract 
    #blank and add LOD and above LOD term to each compound-sample combination, 
    #and keep both blank uncorrected and corrected abundances and LOD.
    if (any(tolower(meta_tb$sample_type)=="blank")) {
      blanks_tb<-excelclean %>%
        filter(sample %in% c(meta_tb %>%
                               filter(tolower(meta_tb$sample_type)=="blank")%>%
                               pull(sample))) %>%
        pivot_longer(cols = 2:ncol(.),
                     names_to = "compound", 
                     values_to= "abundance")%>%
        group_by(compound)%>%
        summarise(av_blank=mean(abundance),
                  LOD=av_blank+3*sd(abundance))%>%
        mutate(
          av_blank=if_else(compound %in% intstds,0,av_blank),
          LOD=if_else(compound %in% intstds,0,LOD)
        )
      excelclean <- excelclean %>%
        pivot_longer(cols = 2:ncol(.),
                     names_to = "compound", 
                     values_to= "abundance")%>%
        left_join(blanks_tb, by=join_by(compound))%>%
        mutate(ab_blankcor=abundance-av_blank,
               LOD_blankcor=LOD-av_blank,
               detected=abundance>LOD) %>%
        pivot_longer(cols = c("abundance","ab_blankcor"),
                     names_to = "datatype", 
                     values_to= "value") %>%
        select(sample,compound,datatype,value,LOD,LOD_blankcor,detected)%>%
        filter(!sample %in% c(meta_tb %>%
                                filter(tolower(meta_tb$sample_type)=="blank")%>%
                                pull(sample)),
               !grepl("lod",tolower(sample)))
      
      long_format<-T
    } else {
      #if no mocks, don't do blank correction
      print(paste0("No samples indicated as blank in metadata column sample , ",
                   "type. Assumed no LOD calculation or blank correction needed"))
    }
  } else {
    #remove internal standard columns and mock samples for data other than
    #abundance
    excelclean<-excelclean %>%
      select(-any_of(intstds)) %>%
      filter(!sample %in% c(meta_tb %>%
                             filter(tolower(meta_tb$sample_type)=="blank")%>%
                             pull(sample)))
  }
  
  #get data in long format if still needed, for corrected isotopologues
  #make one entry per isotopologue
  if (!long_format) {
    excelclean<-excelclean%>%
      pivot_longer(2:ncol(.),names_to = "compound")%>%
      mutate(datatype=datatypename) %>%
      mutate(
        datatype=if_else(
          grepl("isotopologue",tolower(datatype)),
          paste0(datatype,substr(compound,regexpr("_",compound,fixed = T),
                                 nchar(compound))),
          datatype),
        compound=if_else(
          grepl("isotopologue",tolower(datatype)),
          substr(compound,1,regexpr("_",compound,fixed = T)-1),
          compound)
      )
    
    long_format<-T
  }
  
  return(excelclean)
}


#Update library from el-maven file based upon it, with up to date retention times-----------------------------------------
#specify input settings
rawfolderpath<-r"(F:\Projects\Own projects\_Method dev\SICRIT\SICRIT opt after Deans\20241115 MCF002044-47 homogenates\QEx)" #folder with files, need this command to properly read in backslashes
excelfilename<-"MCF002044-47 homogenates QEx labeling results"
libfilename<-NULL  #needed if renaming excel compound names to something else

rawfolderpath<-r"(F:\Projects\Own projects\_Method dev\SICRIT\SICRIT opt after Deans\20241115 MCF002044-47 homogenates\SICRIT)" #folder with files, need this command to properly read in backslashes
excelfilename<-"MCF002044-47_homogenates_SICRIT labeling"
libfilename<-"libTMS_Routine_HP-5MS_deri_SDC_20241129"  #needed if renaming excel compound names to something else

#load data and library
folderpath<-gsub("\\\\", "/", rawfolderpath)         #get correct filepath from raw reference in input
excelfile<-paste0(folderpath,"/",excelfilename,".xlsx")
all_sheets<-excel_sheets(excelfile)
print(paste0("The excel files contains the following sheets: ",
             paste0(all_sheets,collapse = ", ")))

#read library
if (length(libfilename)>0) {
  libfile<-paste0(folderpath,"/",libfilename,".csv")
  lib_tb<-read_csv(libfile)
} else lib_tb<-NULL

#get metadata from excel
metasheetname<-all_sheets[which(grepl("meta",tolower(all_sheets)))]
if(length(metasheetname)==0) {
  print("No sheets with metadata, add a metadata sheet with meta in its name.")
}
if(length(metasheetname)>1) {
  print("Multiple potential metadata sheets, took first sheet with meta in name.")
  metasheetname<-metasheetname[1]
}
print(paste0("Name of sheet used for metadata: ",metasheetname))
long_tb<-meta_tb <- read_excel(excelfile,metasheetname)%>%
  rename(sample=1)


#get abundance data from excel and join to metadata
abundsheetname<-all_sheets[which(grepl("abundance",tolower(all_sheets)))]
if(length(abundsheetname)>0) {
  if(length(abundsheetname)>1) {
    print("Multiple potential abundance sheets, took first sheet with abundance in name.")
    abundsheetname<-abundsheetname[1]
  }
  print(paste0("Name of sheet used for abundance: ",abundsheetname))

  abundance_tb<-read_dollysheet_to_long(excelfile,abundsheetname,
                                        datatypename = "Abund",
                                        lib_tb,
                                        meta_tb = meta_tb)
  long_tb<-abundance_tb %>%
    inner_join(meta_tb,by = join_by(sample))
    # right_join(meta_tb %>%
    #              filter(!sample %in% c(meta_tb %>%
    #                                      filter(tolower(meta_tb$sample_type)=="blank")%>%
    #                                      pull(sample))),
    #            by = join_by(sample))
  
} else {
  print("No sheets with abundance, add an abundance sheet with abundance in its name.")
}

#get fraccon data from excel and join to existing data
fracconsheetname<-all_sheets[which(grepl("contribution",tolower(all_sheets)))]
if(length(fracconsheetname)>0) {
  if(length(fracconsheetname)>1) {
    print("Multiple potential fracContribution sheets, took first sheet with contribution in name.")
    fracconsheetname<-fracconsheetname[1]
  }
  print(paste0("Name of sheet used for fracContribution: ",fracconsheetname))

  fraccon_tb<-read_dollysheet_to_long(excelfile,fracconsheetname,
                                      datatypename = "FracCont",
                                      lib_tb,meta_tb = meta_tb)

  long_tb<-fraccon_tb %>%
    { if(nrow(long_tb)>nrow(meta_tb)) {
      inner_join(.,long_tb %>% 
                  select(-value,-datatype)%>%
                  unique(), by = join_by(sample,compound))%>%
        bind_rows(long_tb)
      } else {
        inner_join(.,meta_tb,by = join_by(sample))
      }
    }
} else {
  print("No sheets with fracContribution, add an fracContribution sheet with contribution in its name.")
}

#get isotopologue data from excel and join to existing data
isotopologuesheetname<-all_sheets[which(grepl("isotopologue",tolower(all_sheets)))]
if(length(isotopologuesheetname)>0) {
  if(length(isotopologuesheetname)>1) {
    print("Multiple potential isotopologue sheets, took first sheet with isotopologue in name.")
    isotopologuesheetname<-isotopologuesheetname[1]
  }
  print(paste0("Name of sheet used for isotopologue: ",isotopologuesheetname))
  
  isotopologue_tb<-read_dollysheet_to_long(excelfile,isotopologuesheetname,lib_tb,
                                           datatypename = "Isotopologues",
                                           meta_tb = meta_tb)

  long_tb<-isotopologue_tb %>%
    { if(nrow(long_tb)>nrow(meta_tb)) {
      inner_join(.,long_tb %>% 
                  select(-value,-datatype)%>%
                  unique(), by = join_by(sample,compound))%>%
        bind_rows(long_tb)
    } else {
      inner_join(.,meta_tb,by = join_by(sample))
    }
    }
} else {
  print("No sheets with isotopologue, add an isotopologue sheet with isotopologue in its name.")
}

write_csv(long_tb,
          paste0(folderpath,"/",excelfilename,"_long.csv"))
# Old function with mock sample recognition for abundances based o --------
#reads in an abundance, fractional contribution or isotopologue sheet of a 
#dolly excel file by name,and prepares the desired table from it

# read_dollysheet_to_long_oldlayoutmocks<-function(excelfile,sheetname,lib_tb,meta_tb) {
#   #read excel sheet data, if derivatized change names to underivatized
#   excel_tb <- read_excel(excelfile,sheetname)
#   
#   if(any(colnames(lib_tb)=="Orig_name")) {
#     for(i in 1:nrow(lib_tb)){
#       colnames(excel_tb)<-sub(lib_tb$compound[i],
#                               gsub("_"," ",lib_tb$Orig_name[i]),
#                               colnames(excel_tb),
#                               fixed = T)
#     }  
#   }
#   
#   #detect internal standards in sheet
#   headers<-excel_tb %>% 
#     select(where(~ all(is.na(.)))) %>%
#     colnames()
#   
#   if (length(headers[which(grepl("internal",tolower(headers)))])>0) {
#     intstdfirstcol<-which(colnames(excel_tb) ==
#                             headers[which(grepl("internal",
#                                                 tolower(headers)))][1])+1
#     intstdlastcol<-which(colnames(excel_tb) ==
#                            headers[which(grepl("internal",
#                                                tolower(headers)))+1][1])-1
#     intstds<-colnames(excel_tb)[intstdfirstcol:intstdlastcol]
#   } else {
#     intstds<-NULL
#   }  
#   
#   #clean excel data and check if abundance sheet for specific actions to 
#   #take with it. Convert data to long format if not done before.
#   excelclean <- excel_tb %>% 
#     select(where(~ !all(is.na(.))))%>%
#     rename(sample=1)
#   
#   long_format<-F
#   if (grepl("abundance",tolower(sheetname))){
#     
#     #if LOD in sheet calculate LOD and average blank from supposed mock samples,
#     #and remove these samples from the abundance sheet. Set blank and LOD to 0
#     #for internal standards
#     #then subtract blank and add LOD and above LOD term to each compound-sample
#     #combination, and keep both blank uncorrected and corrected abundances and
#     #corresponding LODs
#     lod_row<-which(grepl("lod",tolower(excelclean$sample)))
#     if (length(lod_row) >0) {
#       blanks_tb<-excelclean %>%
#         slice(1:(lod_row-1))%>%
#         pivot_longer(cols = 2:ncol(.),
#                      names_to = "compound", 
#                      values_to= "abundance")%>%
#         group_by(compound)%>%
#         summarise(av_blank=mean(abundance),
#                   LOD=av_blank+3*sd(abundance))%>%
#         mutate(
#           av_blank=if_else(compound %in% intstds,0,av_blank),
#           LOD=if_else(compound %in% intstds,0,LOD)
#         )
#       
#       print(paste0("LOD found in abundance sheet, sample preceding it used as",
#                    " blanks: ",
#                    paste0(slice(excelclean,1:(lod_row-1))%>%
#                             pull(sample),collapse = ", ")))
#       
#       firstdatarow<-which(!is.na(excelclean$sample))[which(
#         which(!is.na(excelclean$sample))>lod_row)][1]
#       
#       excelclean <- excelclean %>%
#         slice(firstdatarow:n()) %>%
#         pivot_longer(cols = 2:ncol(.),
#                      names_to = "compound", 
#                      values_to= "abundance")%>%
#         left_join(blanks_tb)%>%
#         mutate(ab_blankcor=abundance-av_blank,
#                LOD_blankcor=LOD-av_blank,
#                detected=abundance>LOD) %>%
#         pivot_longer(cols = c("abundance","ab_blankcor"),
#                      names_to = "datatype", 
#                      values_to= "value") %>%
#         select(sample,compound,datatype,value,LOD,detected)
#       
#       long_format<-T
#     } else {
#       #if no LOD row, don't do blank correction and put in long format
#       print(paste0("No sample with LOD in name found in abundance sheet, ",
#                    "assumed no LOD calculation or blank correction needed"))
#     }
#   } else {
#     #remove internal standard columns for data other than abundance
#     excelclean<-excelclean %>%
#       select(-any_of(intstds))
#   }
#   
#   #get data in long format if still needed, for corrected isotopologues
#   #make one entry per isotopologue
#   if (!long_format) {
#     excelclean<-excelclean%>%
#       pivot_longer(2:ncol(.),names_to = "compound")%>%
#       mutate(datatype=sheetname) %>%
#       mutate(
#         datatype=if_else(
#           grepl("correctedIsotopologues",datatype),
#           paste0(datatype,substr(compound,regexpr("_",compound,fixed = T),
#                                  nchar(compound))),
#           datatype),
#         compound=if_else(
#           grepl("correctedIsotopologues",datatype),
#           substr(compound,1,regexpr("_",compound,fixed = T)-1),
#           compound)
#       )
#     
#     long_format<-T
#   }
#   
#   return(excelclean)
# }

