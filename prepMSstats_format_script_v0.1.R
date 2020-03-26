## Script for prepping two-groups comparisons from measurements with complex experimetal designs ----
## For MSstats (label-free) or MSstatsTMT 
## Miguel Cosenza. version 0.1 ----
##
## README section ----
# Use this script if your experimental design consists of more than 2 condtions and you want to 
# select two for comparison using ROAST enrichment 
#
#
## PLEASE MODIFY DE INPUT IN THE SCRIPT BY ANSWERING THE QUESTIONS ----

## *-1. DID YOU USE "MSstats" OR "MSstatsTMT"?

MSstatsType <- "MSstats"

## *-2 WHICH TWO EXPERIMENTAL CONDITIONS DO YOU WANT TO COMPARE WITH ROAST ENRICHMENT? ----
# Please write them exactly as they appear in your annotation file 

Filt_Cond <- c("less18m", "more24m")

## Load required packages ----

library(dplyr)

## Load or create the MSstats-processed data ----

if(file.exists("Data/ProcessedMSstatsData.Rda") | file.exists("Data/ProcessedMSstatsData.rds") == TRUE){
            
            if(file.exists("Data/ProcessedMSstatsData.Rda") == TRUE){
                        loadRData <- function(fileName){
                                    load(fileName)
                                    get(ls()[ls() != "fileName"])
                        }
                        
                        suppressWarnings(
                                    processed_msstats <- loadRData("Data/ProcessedMSstatsData.Rda"))           
                        
            } else if (file.exists("Data/ProcessedMSstatsData.rds") == TRUE){
                        
                        suppressWarnings(
                                    processed_msstats <-  readRDS(here::here("Data/ProcessedMSstatsData.rds")))
                        
                        
            }
            
            # Load annotation data for MSstats and creation of the design output. 
            
            annotation <- read.delim(file = here::here("Data/annotation.tsv"),
                                     stringsAsFactors = FALSE,
                                     sep = "\t")
            
} else {
            # Load raw data for processing ----
            
            annotation <- read.delim(file = here::here("Data/annotation.tsv"),
                                     stringsAsFactors = FALSE,
                                     sep = "\t")
            
            proteinGroups <- read.delim(file = here::here("Data/proteinGroups.txt"),
                                        stringsAsFactors = FALSE,
                                        sep = "\t")
            
            evidence <- read.delim(file = here::here("Data/evidence.txt"),
                                   stringsAsFactors = FALSE,
                                   sep = "\t")
            
            if(MSstatsType == "MSstatsTMT"){
            require(MSstatsTMT)            
                        ## Run MSstatsTMT formating protein summarization ----
                        
                        processed_msstats <- MaxQtoMSstatsTMTFormat(evidence = evidence,
                                                                    proteinGroups = proteinGroups,
                                                                    annotation = annotation,
                                                                    which.proteinid = "Leading.razor.protein",
                                                                    rmProt_Only.identified.by.site = TRUE)
                        
                        processed_msstats <- proteinSummarization(processed_msstats)
                        
                        saveRDS(object = processed_msstats,
                                file = here::here("Data/ProcessedMSstatsData.rds"))
                        
                        
            } else if (MSstatsType == "MSstats"){
            require(MSstats)
                        ## Run MSstats formating and normalization ----
                        
                        processed_msstats <- MaxQtoMSstatsFormat(evidence = evidence,
                                                                 annotation = annotation,
                                                                 proteinGroups = proteinGroups,
                                                                 proteinID = "Leading.razor.protein",
                                                                 useUniquePeptide = TRUE)
                        
                        processed_msstats <- dataProcess(processed_msstats)
                        
                        saveRDS(object = processed_msstats,
                                file = here::here("Data/ProcessedMSstatsData.rds"))
                        
            }

}

## Filter your annotation file and MSstats-processed data ----
# This is for having only two conditions to evaluate by ROAST 

if(MSstatsType == "MSstats"){
            
            annotation1 <- annotation %>%
                        dplyr::filter(Condition %in% Filt_Cond) %>% 
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Filt_Cond)) 
            write.table(x = annotation1,
                        file = here::here("Data/annotation_to_ROAST.tsv"),
                        sep = "\t")
            
            lfproc <- processed_msstats$RunlevelData %>% 
                                    dplyr::filter(GROUP_ORIGINAL %in% Filt_Cond) %>% 
                                    dplyr::mutate(GROUP_ORIGINAL = factor(GROUP_ORIGINAL,
                                                                          levels = Filt_Cond))
            
            
            processed_msstats$RunlevelData <- lfproc
            
            saveRDS(object = processed_msstats,
                    file = here::here("Data/ProcessedMSstatsData_to_ROAST.rds"))
            
} else if (MSstatsType == "MSstatsTMT"){
            
            annotation1 <- annotation %>%
                        dplyr::filter(Condition %in% Filt_Cond) %>% 
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Conditions)) 
            write.table(x = annotation1,row.names = FALSE,
                        file = here::here("Data/annotation_to_ROAST.tsv"),
                        sep = "\t")
            
            processed_msstats <- processed_msstats %>% 
                        dplyr::filter(Condition %in% Filt_Cond) %>% 
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Conditions))
            
            saveRDS(object = processed_msstats,
                    file = here::here("Data/ProcessedMSstatsData_to_ROAST.rds"))
            
            
}

