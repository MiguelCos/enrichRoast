# Function to get the normalized data from MSstats/MSstatsTMT; format it ----
# to ROAST and generate the design file 

MSstatsToRoast <- function(data,
                           MSstatsAnnotation,
                           MSstatsType, # "MSstats" or "MSsstatsTMT"
                           Paired = FALSE, # if paired design, set to TRUE
                           Conditions # a character vector of two elements corresponding to the same experimental conditions of your setting in MSstats 
                           ) {
      
      if (MSstatsType == "MSstats"){
            if (Paired == FALSE){
                  annot1 <- MSstatsAnnotation %>% 
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Conditions)) %>%
                        dplyr::mutate(Sample_ID = paste0(Condition,"_",BioReplicate)) %>% 
                        dplyr::select(Sample_ID, BioReplicate, Condition) %>% 
                        dplyr::arrange(Condition) %>%
                        dplyr::distinct()
                  
                  lfproc1 <- data$RunlevelData %>% 
                        dplyr::mutate(GROUP_ORIGINAL = factor(GROUP_ORIGINAL,
                                                              levels = Conditions),
                                      Protein = as.character(Protein),
                                      Sample_ID = paste0(GROUP_ORIGINAL,"_",SUBJECT_ORIGINAL)) %>%
                        dplyr::select(Sample_ID, Protein, LogIntensities)
                  
                  widedata <- tidyr::pivot_wider(lfproc1,
                                                 names_from = Sample_ID,
                                                 values_from = LogIntensities,
                                                 values_fn = list(LogIntensities = median)) %>% 
                        dplyr::select(ID = Protein,
                                      annot1$Sample_ID)
                  
                  groups <- factor(annot1$Condition,
                                   levels = Conditions)
                  
                  design <- model.matrix(~groups)
                  
                  roast_input <- list(tabular_data = widedata,
                                      design = design)
                  
                  return(roast_input)
                  
            } else if (Paired == TRUE){
               annot1 <- MSstatsAnnotation %>% 
                  dplyr::mutate(Condition = factor(Condition,
                                                   levels = Conditions)) %>%
                  dplyr::mutate(Sample_ID = paste0(Condition,"_",BioReplicate)) %>% 
                  dplyr::select(Sample_ID, BioReplicate, Condition) %>% 
                  dplyr::arrange(Condition) %>%
                  dplyr::distinct()
               
               lfproc1 <- data$RunlevelData %>% 
                  dplyr::mutate(GROUP_ORIGINAL = factor(GROUP_ORIGINAL,
                                                        levels = Conditions),
                                Protein = as.character(Protein),
                                Sample_ID = paste0(GROUP_ORIGINAL,"_",SUBJECT_ORIGINAL)) %>%
                  dplyr::select(Sample_ID, Protein, LogIntensities)
               
               widedata <- tidyr::pivot_wider(lfproc1,
                                              names_from = Sample_ID,
                                              values_from = LogIntensities,
                                              values_fn = list(LogIntensities = median)) %>% 
                  dplyr::select(ID = Protein,
                                annot1$Sample_ID)
               
               cond <- factor(annot1$Condition,
                              levels = Conditions)
               
               biorep <- factor(annot1$BioReplicate)
               
               design <- model.matrix(~biorep+cond)
               
               roast_input <- list(tabular_data = widedata,
                                   design = design)
               
               return(roast_input)
               }
            
      } else if(MSstatsType == "MSstatsTMT"){
            if (Paired == FALSE) {
                  
                  annot1 <- dplyr::filter(MSstatsAnnotation,
                                          Condition != "Empty",
                                          Condition != "Norm") %>%
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Conditions)) %>% 
                        dplyr::mutate(Sample_ID = paste0(Condition, "_", BioReplicate)) %>% 
                        dplyr::select(Sample_ID, BioReplicate, Condition) %>% 
                        dplyr::arrange(Condition) %>% 
                        dplyr::distinct()
                  
                  tmtproc <- data %>% 
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Conditions),
                                      Protein = as.character(Protein),
                                      Sample_ID = paste0(Condition,"_",BioReplicate)) %>% 
                        dplyr::select(Sample_ID, Protein, Abundance)
                  
                  widedata <- tidyr::pivot_wider(tmtproc,
                                                 names_from = Sample_ID,
                                                 values_from = Abundance,
                                                 values_fn = list(Abundance = median)) %>%
                        dplyr::select(ID = Protein, annot1$Sample_ID)
                  
                  groups <- factor(annot1$Condition,
                                   levels = Conditions)
                  
                  design <- model.matrix(~groups)
                  
                  roast_input <- list(tabular_data = widedata,
                                      design = design)
                  
                  return(roast_input)
                  
            } else if(Paired == TRUE){
                  
                  annot1 <- dplyr::filter(MSstatsAnnotation,
                                          Condition != "Empty",
                                          Condition != "Norm") %>%
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Conditions)) %>%
                        dplyr::mutate(Sample_ID = paste0(Condition,"_",BioReplicate)) %>%
                        dplyr::select(Sample_ID, BioReplicate, Condition) %>%
                        dplyr::arrange(BioReplicate) %>%
                        dplyr::distinct()
                  
                  tmtproc <- data %>% 
                        dplyr::mutate(Condition = factor(Condition,
                                                         levels = Conditions),
                                      Protein, as.character(Protein),
                                      Sample_ID = paste0(Condition,"_",BioReplicate)) %>%
                        dplyr::select(Sample_ID, Protein, Abundance)
                  
                  widedata <- tidyr::pivot_wider(tmtproc,
                                                 names_from = Sample_ID,
                                                 values_from = Abundance,
                                                 values_fn = list(Abundance = median)) %>% 
                        dplyr::select(ID = Protein,
                                      annot1$Sample_ID)
                  
                  cond <- factor(annot1$Condition,
                                 levels = Conditions)
                  
                  biorep <- factor(annot1$BioReplicate)
                  
                  design <- model.matrix(~biorep+cond)
                  
                  roast_input <- list(tabular_data = widedata,
                                      design = design)
                  
                  return(roast_input)
                  
            } else {stop("Please check your input for the 'Paired' argument")}
            
      } else {
            stop("Error: check your MSstatsType input!")
      }
} 
