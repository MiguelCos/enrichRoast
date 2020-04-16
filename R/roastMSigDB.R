## roastMSigDB function ----
roastMSigDB <- function(data,
                        geneIDtype = "SYMBOL",
                        orgDB = "org.Hs.eg.db",
                        species = "Homo sapiens", # this can be any resulting from calling msigdbr::msigdbr_show_species()
                        category = "H", # Any of the main categories presented here: https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
                        subcategory = NULL,
                        specific_category = NULL,
                        design,
                        n_rotations = 999,
                        minSetSize = 1,
                        maxSetSize = 200,
                        pvalueCutoff = 0.05,
                        cutoff_by = "FDR", # one of "FDR" or "PValue"
                        Paired = FALSE){
      
      ## Load required packages ----  
      require(orgDB, character.only = TRUE) || stop(paste("package", orgDB, "is required", sep=" "))
      require(limma) || stop("Package limma is required")
      require(dplyr) || stop("Package dplyr is required")
      require(msigdbr) || stop("Package msigdbr is required")
      require(stringr) || stop("Package stringr is required")

      
      ## Generate matrix for roast ----
      
      if(geneIDtype != "SYMBOL" & geneIDtype != "ENTREZID"){
      genenames <- data$ID
      
      suppressMessages(suppressWarnings(
      mapentrez <- clusterProfiler::bitr(geneID = genenames,
                                         fromType = geneIDtype,
                                         toType = "SYMBOL",
                                         OrgDb = eval(as.name(orgDB))) %>% na.omit()
      ))
      
      
      names(data) <- c(geneIDtype, names(data)[2:eval(dim(data)[2])])
      
      data <- dplyr::left_join(data,
                               mapentrez,
                               by = geneIDtype)
      
      premat2 <- dplyr::filter(data, is.na(eval(as.name(geneIDtype))) == FALSE) %>% 
         dplyr::distinct() %>% na.omit() %>% 
         dplyr::rename(ID = SYMBOL) %>% 
         dplyr::select(ID,2:eval(dim(data)[2]-1), -1) %>%
         dplyr::filter(!ID %in% .$ID[which(duplicated(.$ID))])
      
      genesindata <- premat2$ID
      
      premat3 <- dplyr::select(premat2,
                               -ID) 
      
      matrix1 <- as.matrix(premat3)
      
      row.names(matrix1) <- genesindata
      
      } else {
         
         premat2 <- dplyr::filter(data, is.na(eval(as.name(geneIDtype))) == FALSE) %>% 
            dplyr::distinct() %>% na.omit() %>% 
            dplyr::rename(ID = SYMBOL) %>% 
            dplyr::select(ID,2:eval(dim(data)[2]-1), -1) %>%
            dplyr::filter(!ID %in% .$ID[which(duplicated(.$ID))])
         
         
         genesindata <- premat2$ID
         
         premat3 <- dplyr::select(premat2,
                                  -ID) 
         
         matrix1 <- as.matrix(premat3)
         
         row.names(matrix1) <- genesindata
      }
      
      ## Prep index for roast ----
      
      msigtab <- msigdbr::msigdbr(species = species,
                                  category = category,
                                  subcategory = subcategory)
      
      if (isEmpty(specific_category) == FALSE){
         
         msigtab <- dplyr::filter(msigtab,
                                  str_detect(gs_name, paste0("^",specific_category,"_")))
         
      } 
      
      pathterm_n_iddf <- dplyr::select(msigtab,
                                       gs_name, gs_id) %>% distinct()
                        
      
      list_msig <- msigtab %>% split(x = .$gene_symbol, 
                                     f = .$gs_id)
      
      lenindex <- sapply(list_msig, length)
      
      sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
      list_msig <- list_msig[sublogi1]  
      
      index <- limma::ids2indices(gene.sets = list_msig,
                                  identifiers = genesindata,
                                  remove.empty = TRUE)
      
      
      index2id <- tibble(index = seq_along(genesindata),
                         ID = genesindata)
      
      genesinterm <- qdapTools::list2df(index,
                                        col1 = "index",
                                        col2 = "gs_id") %>% 
         dplyr::left_join(., index2id, by = "index") %>% 
         dplyr::select(-index) 
      
      ## Run roast ----
      
      roast_out <- mroast(y = matrix1,
                         contrast= ncol(design),
                         design = design, 
                         nrot = n_rotations, 
                         index = index)
      
      # Process ROAST output ----
      suppressWarnings(suppressMessages(
         
         roast_out2 <- dplyr::mutate(roast_out,
                                     gs_id = row.names(roast_out)) %>% 
            dplyr::left_join(.,pathterm_n_iddf,
                             by = "gs_id") %>% 
            dplyr::rename(CategoryID = gs_id, CategoryTerm = gs_name)
         
      ))
      
      roast_out2 <- dplyr::filter(roast_out2,
                                  eval(as.name(cutoff_by)) <= pvalueCutoff)
      
      fdrnterm <- dplyr::select(roast_out2,
                                CategoryTerm, FDR, PValue, NGenes)
      
      
      genesinterm <- dplyr::filter(genesinterm,
                                    gs_id %in% roast_out2$CategoryID) %>% 
                     dplyr::rename(ID = 2)
      
      suppressWarnings(suppressMessages(
         genesintermread <-  genesinterm %>% 
            left_join(., pathterm_n_iddf,
                      by = "gs_id") %>% dplyr::select(ID, gs_id, gs_name)
      ))
      
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
      
      if(Paired == TRUE){
      suppressWarnings(
         suppressMessages(
            log2FCs <- dplyr::mutate(limma_tab,
                                     ID = row.names(limma_tab))  %>% 
               dplyr::filter(ID %in% genesinterm$ID) %>% 
               dplyr::select(log2FC = eval(dim(.)[2]-5),ID) %>% 
               dplyr::left_join(., genesintermread, by = "ID")  %>%
               dplyr::select(ID, log2FC, CategoryID = gs_id, CategoryTerm = gs_name) %>% 
               dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
               dplyr::filter(is.na(FDR) == FALSE)
         ))
      } else if(Paired == FALSE){
         suppressWarnings(
            suppressMessages(
               log2FCs <- dplyr::mutate(limma_tab,
                                        ID = row.names(limma_tab))  %>% 
                  dplyr::filter(ID %in% genesinterm$ID) %>% 
                  dplyr::select(log2FC = logFC, ID) %>% 
                  dplyr::left_join(., genesintermread, by = "ID")  %>%
                  dplyr::select(ID, log2FC, CategoryID = gs_id, CategoryTerm = gs_name) %>% 
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
                          GenesPerTerm = genesintermread,
                          log2FCs = log2FCs)
      
      if (dim(roast_out2)[1] == 0){warning("No terms were enriched under your p.adjusted-value threshold")}
      return(roastResult)
      
}
                        