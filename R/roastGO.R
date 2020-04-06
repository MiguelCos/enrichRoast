## roastGO function ----

#data = tabular_data
#geneIDtype = geneIDtype
#orgDB = orgDB
#design = design
#n_rotations = 99
#minSetSize = minSetSize
#maxSetSize = maxSetSize
#pvalueCutoff = pvalueCutoff
#exclusionList = exclusionList
#Paired = FALSE
#ontology = "MF"

roastGO <- function(data,
                    geneIDtype = "SYMBOL",
                    ontology = "MF",
                    orgDB = "org.Hs.eg.db", 
                    design,
                    n_rotations = 9999,
                    minSetSize = 1,
                    maxSetSize = 1506,
                    pvalueCutoff = 0.05,
                    Paired = FALSE){
      
      ## Load required packages ----
      require(orgDB, character.only = TRUE) || stop(paste("package", organism, "is required", sep = " "))
      require(limma) || stop("Package limma is required")
      require(GO.db) || stop("Package reactome.db is required")
      require(AnnotationDbi)
      require(dplyr)
      
      ## Generate matrix for roast ----
   
   premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
      dplyr::distinct() %>% na.omit()
   
   genesindata <- premat2$ID
   
   premat3 <- dplyr::select(premat2,
                            -ID) 
   
   matrix1 <- as.matrix(premat3)
   
   row.names(matrix1) <- genesindata
   
   ## Prep index for roast ----
   
      suppressWarnings(
      suppressMessages(
   goterm_n_ontol <- AnnotationDbi::Ontology(GO.db::GOTERM)
      ))
   
   gos_to_test <- goterm_n_ontol[goterm_n_ontol == ontology]  
   
      suppressWarnings(
      suppressMessages(
   goterm_n_id <- AnnotationDbi::mapIds(GO.db,
                                        keys = names(gos_to_test),
                                        column = "TERM",
                                        keytype = "GOID")
      ))
      
   goterm_n_iddf <- data.frame(GOID = names(goterm_n_id),
                               GOTERM = goterm_n_id)
   
   
   golist <- suppressMessages(
      AnnotationDbi::mapIds(eval(as.name(orgDB)), keys=names(gos_to_test), 
                            column=geneIDtype,
                            keytype="GOALL", multiVals='list'))
   
   lenindex <- sapply(golist, length)
   
   sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
   golist <- golist[sublogi1]  
   
   index <- limma::ids2indices(gene.sets = golist,
                               identifiers = genesindata,
                               remove.empty = TRUE)
   
   index2id <- tibble(index = seq_along(genesindata),
                      ID = genesindata)
   
   genesinterm <- qdapTools::list2df(index,
                                     col1 = "index",
                                     col2 = "GOID") %>% 
                  left_join(., index2id, by = "index") %>% 
                  dplyr::select(-index)
   
   
   if(geneIDtype != "SYMBOL"){
   
   suppressWarnings(
   suppressMessages(
   symb1 <- clusterProfiler::bitr(genesinterm$ID,
                                  fromType = geneIDtype,
                                  toType = "SYMBOL",
                                  OrgDb = eval(as.name(orgDB)),
                                  drop = FALSE) %>% 
               dplyr::rename(ID = geneIDtype)
   ))
   
   suppressWarnings(
      suppressMessages(
   genesinterm <- left_join(genesinterm, symb1,
                                by = "ID") %>% 
      left_join(., goterm_n_iddf,
                by = "GOID")
      ))
   }
   # Run ROAST ----
   
   roast_out <- mroast(y = matrix1,
                    contrast= ncol(design),
                    design = design, 
                    nrot = n_rotations, 
                    index = index)
   
   # Process ROAST output ----
   suppressWarnings(
      suppressMessages(
   roast_out2 <- dplyr::mutate(roast_out,
                               GOID = row.names(roast_out)) %>% 
      dplyr::left_join(.,goterm_n_iddf,
                       by = "GOID") %>% 
      dplyr::rename(CategoryID = GOID, CategoryTerm = GOTERM)
   
      ))
   
   roast_out2 <- dplyr::filter(roast_out2,
                               FDR <= pvalueCutoff)
   
   fdrnterm <- dplyr::select(roast_out2,
                             CategoryTerm, FDR, NGenes)
   
   genesinterm <- dplyr::filter(genesinterm,
                                GOID %in% roast_out2$CategoryID)
   
   # Run limma to get log2FC values per protein and category ----
   
   suppressWarnings(
      suppressMessages(
   limma_out <- lmFit(object = matrix1,
                      design = design)
      ))
   
   limma_out2 <- eBayes(limma_out)
   
   suppressWarnings(
      suppressMessages(
   limma_tab <- topTable(limma_out2, number = dim(matrix1)[1])
      ))
   
   # Get log2FC information from Limma and reformat output ----
   if (Paired == TRUE){
   suppressWarnings(
      suppressMessages(
   log2FCs <- dplyr::mutate(limma_tab,
                            ID = row.names(limma_tab)) %>%  
      dplyr::select(ID, log2FC = eval(dim(.)[2]-5)) %>% 
      dplyr::left_join(., genesinterm, by = "ID")  %>%
      dplyr::left_join(., goterm_n_iddf, by = c("GOID", "GOTERM")) %>%
      dplyr::rename(CategoryID = GOID, CategoryTerm = GOTERM) %>%
      dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
      dplyr::filter(is.na(FDR) == FALSE)
      ))
   } else if (Paired == FALSE){
      suppressWarnings(
         suppressMessages(
            log2FCs <- dplyr::mutate(limma_tab,
                                     ID =row.names(limma_tab)) %>%  
               dplyr::select(ID, log2FC = logFC) %>% 
               dplyr::left_join(., genesinterm, by = "ID")  %>%
               dplyr::left_join(., goterm_n_iddf, by = c("GOID", "GOTERM")) %>%
               dplyr::rename(CategoryID = GOID, CategoryTerm = GOTERM) %>%
               dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
               dplyr::filter(is.na(FDR) == FALSE)
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
                       ontology = ontology,
                       organism = orgDB)
   
    
   if (dim(roast_out2)[1] == 0){warning("No terms were enriched under your p.adjusted-value threshold",
                                        .call = FALSE)}
   
   return(roastResult)
   
   
}

