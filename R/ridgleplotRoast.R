## ridgeplotRoast function ----

ridgeplotRoast <- function(roastOutput,
                           show_n_terms = 25,
                           colorby = "FDR"){
      
      # Load required packages 
      
      require(ggplot2) || stop("Package ggplot2 is required")
      require(dplyr) || stop("Package dplyr is required")
      require(tidyr) || stop("Package tidyr is required")
      require(ggridges) || stop("Package ggridges is required")
      
      # Prep data ----
      
      tofil <- roastOutput$roastOutput
      toridge <- roastOutput$log2FCs
      
      catid2PValue <- dplyr::select(tofil,
                                    CategoryID, PValue)
      
      toproplot <- dplyr::select(tofil,
                                 NGenes, PropDown, PropUp, Direction, CategoryTerm,
                                 FDR, PValue) %>%
         dplyr::top_n(n = show_n_terms,
                      wt = NGenes) %>%
         dplyr::mutate(DiffProp = abs(PropUp - PropDown),
                       PropDown = -PropDown) %>%
         tidyr::pivot_longer(cols = c(PropDown, PropUp),
                             names_to = "PropDirection",
                             values_to = "Proportion") %>%
         dplyr::group_by(CategoryTerm, NGenes) 
      
      datatab <- dplyr::filter(toridge,
                               CategoryTerm %in% unique(toproplot$CategoryTerm)) %>% 
         dplyr::left_join(.,catid2PValue, by = "CategoryID") %>%
         dplyr::arrange(-NGenes) 
   
      # Plot ----
      
      ridges <- ggplot(data = datatab, aes(x = log2FC, y = CategoryTerm, fill = eval(as.name(colorby))))+
         scale_fill_gradient(low = "#477af8", high = "#ff3333", name = colorby)+
         geom_density_ridges()+
         xlab("Log2(Fold-change)")+ 
         ylab("Biological Category")+
         theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
               panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size=1.5),
               axis.title=element_text(size=10,face="bold"))
      
      return(ridges)
}
