
enrichRoastGO = roast_result
cutoff = 0.7
by = "FDR"


simplifyGO <- function(enrichRoastGO, cutoff = 0.7, by = "FDR"){
      
      # Load required packages ----
      require(GOSemSim) || stop("Package GOSemSim is required")
      require(dplyr) || stop("Package dplyr is required")
      require(tidyr) || stop("Package tidyr is required")
      
      
      terms <- enrichRoastGO$roastOutput
      
      selterms <- dplyr::select(terms,
                                go1 = CategoryID, all_of(by))
      
      semdat <- godata(OrgDb = enrichRoastGO$organism, ont = enrichRoastGO$ontology)
      
      
      mat <- mgoSim(terms$CategoryID, terms$CategoryID,
             semData = semdat,
             measure = 'Wang',
             combine = NULL)
      
      matdf <- as.data.frame(mat)
      matdf <- dplyr::mutate(matdf,
                             go1 = row.names(matdf))

      matdf <- tidyr::pivot_longer(matdf,
                                   cols = 1:eval(dim(matdf)[2]-1),
                                   names_to = "go2",
                                   values_to = "similarity") %>% na.omit()
      
      matsim <- dplyr::left_join(matdf, selterms, by = "go1")
      
      
      
      ID <- terms$CategoryID
      
      term2remove <- character()
      
      for (i in seq_along(ID)) {
            ii <- which(matsim$go2 == ID[i] & matsim$similarity > cutoff)
            if (length(ii) < 2)
                  next
            
            sim_subset <- matsim[ii,]
            
            jj <- which(sim_subset[, by] == min(sim_subset[, by]))
            
            term2remove <- c(term2remove, sim_subset$go1[-jj]) %>% unique
      }
      
      enrichRoastGO$roastOutput <- dplyr::filter(terms,
                                                 !CategoryID %in% term2remove)
      
      enrichRoastGO$GenesPerTerm <- dplyr::filter(enrichRoastGO$GenesPerTerm,
                                                  !GOID %in% term2remove)
      
      enrichRoastGO$log2FCs <- dplyr::filter(enrichRoastGO$log2FCs,
                                             !CategoryID %in% term2remove)
      
      return(enrichRoastGO)
 
}
