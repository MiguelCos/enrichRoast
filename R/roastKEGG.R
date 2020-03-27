data = roast_input$tabular_data
geneIDtype = geneIDtype
orgDB = orgDB
organism = organism # here the sintax should correspond with the  sintax
design = roast_input$design
n_rotations = n_rotations
minSetSize = minSetSize
maxSetSize = maxSetSize
pvalueCutoff = pvalueCutoff
exclusionList = exclusionList


roastKEGG <- function(data,
                      geneIDtype = "SYMBOL",
                      orgDB = "org.Hs.eg.db",
                      organism = "hsa", # here the sintax should correspond with the  sintax
                      design,
                      n_rotations = 999,
                      minSetSize = 1,
                      maxSetSize = 1000,
                      pvalueCutoff = 0.05,
                      exclusionList = TRUE) {
      
   ## Load required packages ----
      
      require(KEGGREST) || stop("Package REST is required")
      require(orgDB, character.only = TRUE) || stop(paste("package", orgDb, "is required", sep=" "))
      require(limma) || stop("Package limma is required")
      suppressMessages(require(clusterProfiler)) || stop("Package clusterProfiler is required")
      require(dplyr) || stop("Package dplyr is required")
      require(stringr) || stop("Package stringr is required")
      
      ## Generate matrix to roast ----
      
      
      if (geneIDtype != "ENTREZID"){
         
         genenames <- data$ID
         
         suppressMessages(suppressWarnings(
         mapentrez <- clusterProfiler::bitr(geneID = genenames,
                                            fromType = geneIDtype,
                                            toType = "ENTREZID",
                                            OrgDb = eval(as.name(orgDB))) %>% na.omit()
         ))
         
         names(data) <- c(geneIDtype, names(data)[2:eval(dim(data)[2])])
         
         suppressMessages(suppressWarnings(
         premat1 <- dplyr::left_join(data,
                                     mapentrez,
                                     by = geneIDtype)
         ))
         
         suppressMessages(suppressWarnings(
         premat2 <- dplyr::select(premat1,
                                  -all_of(geneIDtype)) %>%  
            dplyr::rename(Name = ENTREZID) %>% 
            dplyr::filter(is.na(Name) == FALSE) %>% 
            dplyr::distinct() %>% na.omit()
         ))
         
         genesindata <- premat2$Name
         
         premat3 <- dplyr::select(premat2,
                                  -Name) 
         
         matrix1 <- as.matrix(premat3)
         
         row.names(matrix1) <- genesindata
         
      } else {
         
         premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
            dplyr::distinct() %>% na.omit()
         
         genesindata <- premat2$ID
         
         premat3 <- dplyr::select(premat2,
                                  -ID) 
         
         matrix1 <- as.matrix(premat3)
         
         row.names(matrix1) <- genesindata
         
      }    
      
      ## Prep index for roast ----
      
      keggid_to_term <- keggList("pathway", organism)
      
      keggid <- str_remove_all(names(keggid_to_term), "^path:") %>%
         str_trim()
      
      keggterm <- str_remove_all(keggid_to_term, "\\ \\-.*") %>%
         str_trim()
      
      keggidtoterm_df <- data.frame(KEGGID = keggid,
                                    KEGGTERM = keggterm)
      
      if (exclusionList == TRUE){
         getClass <- function(pathway){
            qery <- KEGGREST::keggGet(pathway)
            
            class <- qery[[1]]$CLASS
            
            return(class)
         }
         
         pathclss <- sapply(keggid, getClass)
         
         pathclss <- qdapTools::list2df(pathclss,
                                        col1 = "Class",
                                        col2 = "PathwayID") %>% 
            tidyr::separate(col = Class, into = c("Class", "Subclass"), sep = "; ")
         
         
         exclusion_subclass <- c("Drug resistance: antineoplastic",
                                 "Endocrine system", "Aging",
                                 "Circulatory system", "Xenobiotics biodegradation and metabolism",
                                 "Drug resistance: antineoplastic", "Nervous system", "Sensory system",
                                 "Excretory system", "Digestive system")
         
         exclusion_class <- c("Human Diseases")
         
         excluded <- dplyr::filter(pathclss,
                                   Subclass %in% exclusion_subclass,
                                   Class %in% exclusion_class)
         
         aftrexclud <- dplyr::filter(pathclss,
                                     !Subclass %in% exclusion_subclass,
                                     !Class %in% exclusion_class)
         
         keggid <- subset(keggid, keggid %in% aftrexclud$PathwayID)
         
         keggidtoterm_df <- dplyr::filter(keggidtoterm_df,
                                          KEGGID %in% keggid)
         
         path2entrez <- KEGGREST::keggLink("pathway", organism)
         
         entrezcor <- str_remove_all(names(path2entrez), 
                                     paste0("^",organism,":")) %>% 
            str_trim()
         
         path2entrez <- str_remove_all(path2entrez, "^path:") %>% 
            str_trim()
         
         names(path2entrez) <- entrezcor
         
         path2entrez <- subset(path2entrez, path2entrez %in% keggidtoterm_df$KEGGID)
         
         list_path2entrez <- split(x = names(path2entrez),
                                   f = path2entrez)
         
      } else {
         
         path2entrez <- KEGGREST::keggLink("pathway", organism)
         
         entrezcor <- str_remove_all(names(path2entrez), 
                                     paste0("^",organism,":")) %>% 
            str_trim()
         
         path2entrez <- str_remove_all(path2entrez, "^path:") %>% 
            str_trim()
         
         names(path2entrez) <- entrezcor
         
         path2entrez <- subset(path2entrez, path2entrez %in% keggidtoterm_df$KEGGID)
         
         list_path2entrez <- split(x = names(path2entrez),
                                   f = path2entrez)
         
      }
      
         leindex <- sapply(list_path2entrez, length)
      
         sublogi1 <- between(leindex, minSetSize, maxSetSize) 
         list_path2entrez <- list_path2entrez[sublogi1] 
         
         
         index <- limma::ids2indices(gene.sets = list_path2entrez,
                                     identifiers = genesindata,
                                     remove.empty = TRUE)
         
         index2id <- tibble(index = seq_along(genesindata),
                            ID = genesindata)
         
         genesinterm <- qdapTools::list2df(index,
                                           col1 = "index",
                                           col2 = "KEGGID") %>% 
            left_join(., index2id, by = "index") %>% 
            dplyr::select(-index) %>% 
            rename(ENTREZID = ID)
         
      
            suppressWarnings(
               suppressMessages(
                  symb1 <- clusterProfiler::bitr(genesinterm$ENTREZID,
                                                 fromType = "ENTREZID",
                                                 toType = "SYMBOL",
                                                 OrgDb = eval(as.name(orgDB)),
                                                 drop = FALSE)
               ))
            
            suppressWarnings(
               suppressMessages(
                  genesinterm <- left_join(genesinterm, symb1,
                                           by = "ENTREZID") %>% 
                     left_join(., keggidtoterm_df,
                               by = "KEGGID")
               ))
       
      
      ## Run roast ----
      
      roast_out <- mroast(y = matrix1,
                         contrast= ncol(design),
                         design = design, 
                         nrot = n_rotations, 
                         index = index)
      
      
      ## Process output ----
         suppressWarnings(suppressMessages(
      roast_out2 <- dplyr::mutate(roast_out,
                                  KEGGID = row.names(roast_out)) %>% 
         dplyr::left_join(.,keggidtoterm_df,
                          by = "KEGGID") %>% 
         dplyr::rename(CategoryID = KEGGID, CategoryTerm = KEGGTERM)
         ))
      
      roast_out2 <- dplyr::filter(roast_out2,
                                  FDR <= pvalueCutoff)
      
      fdrnterm <- dplyr::select(roast_out2,
                                CategoryTerm, FDR, NGenes)
      
      genesinterm <- dplyr::filter(genesinterm,
                                       KEGGID %in% roast_out2$CategoryID)
      
      
      # Run limma to get log2FC values per protein and category ----
      
         suppressWarnings(suppressMessages(
      limma_out <- lmFit(object = matrix1,
                               design = design)
         ))
      
      limma_out2 <- eBayes(limma_out)
      
      suppressWarnings(
         suppressMessages(
            limma_tab <- topTable(limma_out2, number = dim(matrix1)[1])
         ))
      
      # Get log2FC information from Limma and reformat output ----
      
      if(Paired == TRUE){
      suppressWarnings(
         suppressMessages(
            log2FCs <- dplyr::mutate(limma_tab,
                                     ENTREZID = row.names(limma_tab)) %>%  
               dplyr::select(ENTREZID, log2FC = eval(dim(.)[2]-5)) %>% 
               dplyr::left_join(., genesinterm, by = "ENTREZID")  %>%
               dplyr::left_join(., keggidtoterm_df, by = c("KEGGID", "KEGGTERM")) %>%
               dplyr::rename(CategoryID = KEGGID, CategoryTerm = KEGGTERM) %>%
               dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
               dplyr::filter(is.na(FDR) == FALSE) %>% 
               dplyr::select(ENTREZID, SYMBOL, log2FC, CategoryID, CategoryTerm, FDR, NGenes)
         ))
      } else if(Paired == FALSE){
         suppressWarnings(
            suppressMessages(
               log2FCs <- dplyr::mutate(limma_tab,
                                        ENTREZID = row.names(limma_tab)) %>%  
                  dplyr::select(ENTREZID, log2FC = logFC) %>% 
                  dplyr::left_join(., genesinterm, by = "ENTREZID")  %>%
                  dplyr::left_join(., keggidtoterm_df, by = c("KEGGID", "KEGGTERM")) %>%
                  dplyr::rename(CategoryID = KEGGID, CategoryTerm = KEGGTERM) %>%
                  dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
                  dplyr::filter(is.na(FDR) == FALSE) %>% 
                  dplyr::select(ENTREZID, SYMBOL, log2FC, CategoryID, CategoryTerm, FDR, NGenes)
            ))
      }
         
         
      suppressWarnings(suppressMessages(
         
         meanFCTerm <- log2FCs %>% 
            dplyr::group_by(CategoryID, CategoryTerm) %>% 
            dplyr::summarise(meanlog2FC = mean(log2FC)) %>% dplyr::ungroup() %>% 
            dplyr::filter(CategoryID %in% roast_out2$CategoryID) %>% 
            dplyr::select(-CategoryTerm, CategoryID)
      ))
      
      suppressWarnings(
         suppressMessages(
            roast_out3 <- dplyr::left_join(roast_out2,
                                           meanFCTerm,
                                           by = "CategoryID")
         ))
      
      
      
      roastResult <- list(roastOutput = roast_out3,
                          GenesPerTerm = genesinterm,
                          log2FCs = log2FCs,
                          exclusionList = excluded)
      
      if(exclusionList == TRUE){
         message("Pathways associated with the next category classes were excluded from the analysis:",
                 exclusion_class,"; ",paste(exclusion_subclass, sep = " ", collapse = "; "))
         
         exclumess <- paste("Pathways associated with the next category classes were excluded from the analysis:",
                            exclusion_class,"; ",paste(exclusion_subclass, sep = " ", collapse = "; "))
         
         roastResult$exclusionList <- excluded
         roastResult$exclusionMessage <- exclumess
      }
      
      if (dim(roast_out2)[1] == 0){warning("No terms were enriched under your p.adjusted-value threshold")}
      return(roastResult)
      
}
