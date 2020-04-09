### ROAST - Reactome function ----

roastReactome <- function(data, 
                          geneIDtype = "SYMBOL", 
                          orgDB = "org.Hs.eg.db",
                          design,
                          n_rotations = 9999,
                          minSetSize = 1,
                          maxSetSize = 1506,
                          pvalueCutoff = 0.05,
                          exclusionList = TRUE,
                          species = "Homo sapiens",
                          Paired = FALSE){
      
      ## Load required packages
      require(orgDB, character.only = TRUE) || stop(paste("package", orgDB, "is required", sep=" "))
      require(limma) || stop("Package limma is required")
      require(reactome.db) || stop("Package reactome.db is required")
      require(clusterProfiler) || stop("Package clusterProfiler is required")
      require(dplyr)
      require(stringr)
   
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
            dplyr::filter(is.na(Name) == FALSE,
                          !Name %in% .$Name[which(duplicated(.$Name))]) %>% 
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
   
   suppressMessages(suppressWarnings(
      toexclude <- AnnotationDbi::select(reactome.db,
                                         keys = keys(reactome.db, keytype = "PATHID"),
                                         keytype =  "PATHID",
                                         columns = c("PATHID", "PATHNAME"))
   ))
   
   if(exclusionList == TRUE){
      
      
      exluded1 <-  toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "disease"))
      
      exluded2 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "Disease"))
      
      exluded3 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "disorder"))
      
      exluded4 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "infection"))
      
      exluded5 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "Infection"))
      
      exluded6 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "Viral"))
      
      exluded9 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "viral"))
      
      exluded7 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "Influenza"))
      
      exluded8 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, paste0(species,": "))) %>%
         dplyr::filter(str_detect(PATHNAME, "influenza"))
      
      reactome_excluded <- bind_rows(exluded1,
                                     exluded2,
                                     exluded3,
                                     exluded4,
                                     exluded5,
                                     exluded6,
                                     exluded7,
                                     exluded8,
                                     exluded9) 
      
      pathexcl <- reactome_excluded$PATHID
      
      reactlistentrez1 <- as.list(reactomePATHID2EXTID)
      
      reactlistentrez <- reactlistentrez1[which(!names(reactlistentrez1) %in% pathexcl)]
      
      lenindex <- sapply(reactlistentrez, length)
      
      sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
      reactlistentrez <- reactlistentrez[sublogi1]  
      
      index <- limma::ids2indices(gene.sets = reactlistentrez,
                                  identifiers = genesindata,
                                  remove.empty = TRUE)
      
      index2id <- tibble(index = seq_along(genesindata),
                         ID = genesindata)
      
      pathterm_n_iddf <- toexclude %>% dplyr::filter(str_detect(PATHNAME, species)) %>% 
                           dplyr::mutate(PATHNAME = str_remove(PATHNAME, ".*: "))
      
      
   } else if(exclusionList == FALSE){
   
      reactlistentrez <- as.list(reactomePATHID2EXTID)
      
      lenindex <- sapply(reactlistentrez, length)
      
      sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
      reactlistentrez <- reactlistentrez[sublogi1]  
      
      index <- limma::ids2indices(gene.sets = reactlistentrez,
                                identifiers = genesindata,
                                remove.empty = TRUE)
   
      index2id <- tibble(index = seq_along(genesindata),
                         ID = genesindata)
      
      pathterm_n_iddf <- toexclude %>% dplyr::filter(str_detect(PATHNAME, species)) %>% 
                           dplyr::mutate(PATHNAME = str_remove(PATHNAME, ".*: "))
   
   } else {stop("Please check your input for 'exclusionList'")}
   
   genesinterm <- qdapTools::list2df(index,
                                     col1 = "index",
                                     col2 = "PATHID") %>% 
                  dplyr::left_join(., index2id, by = "index") %>% 
                  dplyr::select(-index) %>% 
                  dplyr::rename(ENTREZID = ID)
   
   
   
   suppressWarnings(suppressMessages(
      
      symb1 <- clusterProfiler::bitr(genesinterm$ENTREZID,
                                        fromType = "ENTREZID",
                                        toType = "SYMBOL",
                                        OrgDb = eval(as.name(orgDB)),
                                        drop = FALSE) 
               
      ))
   
   suppressWarnings(suppressMessages(
         genesintermread <- left_join(genesinterm, symb1,
                                      by = "ENTREZID") %>% 
            left_join(., pathterm_n_iddf,
                      by = "PATHID")
      ))
   # Run roast ----
   
   roast_out <- mroast(y = matrix1,
                      contrast= ncol(design),
                      design = design, 
                      nrot = n_rotations, 
                      index = index)
   
   # Process ROAST output ----
   suppressWarnings(suppressMessages(
      
         roast_out2 <- dplyr::mutate(roast_out,
                                     PATHID = row.names(roast_out)) %>% 
            dplyr::left_join(.,pathterm_n_iddf,
                             by = "PATHID") %>% 
            dplyr::rename(CategoryID = PATHID, CategoryTerm = PATHNAME)
         
      ))
   
   roast_out2 <- dplyr::filter(roast_out2,
                               FDR <= pvalueCutoff)
   
   nterms <- dim(roast_out2)[1]
   
   fdrnterm <- dplyr::select(roast_out2,
                             CategoryTerm, FDR, NGenes)
   
   genesintermread <- dplyr::filter(genesintermread,
                                    PATHID %in% roast_out2$CategoryID)
   
   # Run limma to get log2FC values per protein and category ----
   
   suppressWarnings(
      suppressMessages(
        limma_out <- lmFit(object = matrix1,
                            design = design)
      ))
   
   limma_out <- eBayes(limma_out)
   
   suppressWarnings(
      suppressMessages(
         limma_tab <- topTable(limma_out, number = dim(matrix1)[1], coef = 2)
      ))
   
   # Get log2FC information from Limma and reformat output ----
   
   if(Paired == TRUE){
   suppressWarnings(
      suppressMessages(
         log2FCs <- dplyr::mutate(limma_tab,
                                  ENTREZID = row.names(limma_tab)) %>% 
            dplyr::filter(ENTREZID %in% symb1$ENTREZID) %>% 
            dplyr::select(log2FC = eval(dim(.)[2]-5),ENTREZID) %>% 
            dplyr::left_join(., genesintermread, by = "ENTREZID") %>% ## Look for a way to extract the z-score values as they are used by the ROAST algorithm
            dplyr::select(log2FC, ENTREZID, SYMBOL, CategoryID = PATHID, CategoryTerm = PATHNAME) %>% 
            dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
            dplyr::filter(is.na(FDR) == FALSE)
      ))
   } else if(Paired == FALSE){
      #suppressWarnings(
         #suppressMessages(
            log2FCs <- dplyr::mutate(limma_tab,
                                     ENTREZID = row.names(limma_tab)) %>% 
               dplyr::filter(ENTREZID %in% symb1$ENTREZID) %>% 
               dplyr::select(log2FC = logFC, ENTREZID) %>% 
               dplyr::left_join(., genesintermread, by = "ENTREZID") %>% ## Look for a way to extract the z-score values as they are used by the ROAST algorithm
               dplyr::select(log2FC, ENTREZID, SYMBOL, CategoryID = PATHID, CategoryTerm = PATHNAME) %>% 
               dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
               dplyr::filter(is.na(FDR) == FALSE)
         #))
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
                       GenesPerTerm = genesintermread,
                       log2FCs = log2FCs)
   
   if(exclusionList == TRUE){
      message("Pathways associated with human diseases, disorders and infections where excluded from de analysis")
      
      exclumess <- "Pathways associated with human diseases, disorders and infections where excluded from de analysis"
      
      roastResult$exclusionList <- reactome_excluded
      roastResult$exclusionMessage <- exclumess
      
   } else {
      exclumess <- "No pathways were excluded from the analysis"
      
      roastResult$exclusionList <- NULL
      roastResult$exclusionMessage <- exclumess
      
   }
   
   
   if (dim(roast_out2)[1] == 0){warning("No terms were enriched under your p.adjusted-value threshold")}
   return(roastResult)
}
