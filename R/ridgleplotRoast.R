## ridgeplotRoast function ----

show_n_terms <- 30
roastresult <- roast_result
colorby <- cutoff_by
top_n_by <- "NGenes" # "NGenes"

ridgeplotRoast <- function(roastresult,
                           show_n_terms = 25,
                           colorby = "FDR",
                           top_n_by = "NGenes"){ # one of "Difference" or "NGenes"
      
      # Load required packages 
      
      require(ggplot2) || stop("Package ggplot2 is required")
      require(dplyr) || stop("Package dplyr is required")
      require(tidyr) || stop("Package tidyr is required")
      require(ggridges) || stop("Package ggridges is required")
      
      # Prep data ----
      
      tofil <- roastresult$roastOutput
      toridge <- roastresult$log2FCs
      
      #catid2PValue <- dplyr::select(tofil,
      #                              CategoryID, PValue)
      
      toproplot <- dplyr::select(tofil,
                                 NGenes, PropDown, PropUp, Direction, CategoryTerm,
                                 FDR, PValue) %>%
         #dplyr::top_n(n = show_n_terms,
         #             wt = NGenes) %>%
         dplyr::mutate(DiffProp = abs(PropUp - PropDown),
                       PropDown = -PropDown) %>%
         dplyr::top_n(n = show_n_terms,
                      wt = if(top_n_by == "Difference"){DiffProp}
                      else if(top_n_by == "NGenes"){NGenes}
         ) %>%
         tidyr::pivot_longer(cols = c(PropDown, PropUp),
                             names_to = "PropDirection",
                             values_to = "Proportion") %>%
         dplyr::group_by(CategoryTerm, NGenes) 
      
      datatab <- dplyr::filter(toridge,
                               CategoryTerm %in% unique(toproplot$CategoryTerm)) #%>% 
         #dplyr::left_join(.,catid2PValue, by = "CategoryID") #%>%
      
      summtab <-  dplyr::group_by(datatab, CategoryTerm) %>%
         dplyr::summarise(meadianlo2FC = median(log2FC))
      
      datatab <- left_join(datatab, summtab, by = "CategoryTerm") %>% ungroup() %>%
         dplyr::arrange(-meadianlo2FC) %>% filter(!NGenes <= 2) %>%
         dplyr::mutate(FDR = round(FDR, 2),
                       PValue = round(PValue, 2))
      
      zero_range <- function(x) {
         if (length(x) == 1) return(TRUE)
         x <- range(x) / mean(x)
         isTRUE(all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5))
      }
      
      pvals <- dplyr::pull(datatab, eval(as.name(colorby)))
      
      if(isEmpty(pvals)){stop("Error: no terms to plot")}
      
      if (zero_range(pvals) == TRUE){
         maxpval <- max(pvals)
         limits <- c(0,maxpval)
         breaks <- round(seq(0, maxpval, length = 7), 2)
      } else if(zero_range(pvals) == FALSE){
         maxpval <- max(pvals)
         minpval <- min(pvals)
         limits <- c(minpval,maxpval)
         breaks <- round(seq(minpval, maxpval, length = 7), 2)
      }
      
      # Plot ----
      
      ridges <- ggplot(data = datatab, aes(x = log2FC, y = fct_reorder(CategoryTerm, meadianlo2FC), 
                                           fill = eval(as.name(colorby))))+
         scale_fill_gradient(high = "#0fabbc",
                             low = '#fa163f',
                             guide = guide_colourbar(reverse = TRUE),
                             name = colorby,
                             limits = limits,
                             breaks = breaks)+
         geom_density_ridges()+
         xlab("Log2(Fold-change)")+ 
         ylab("Biological Category")+
         labs(caption = if(top_n_by == "Difference"){paste0("Showing top ",show_n_terms," terms by |ProportionUp - ProportionDown|")}
              else if(top_n_by == "NGenes"){paste0("Showing top ",show_n_terms," terms by N Genes per set")}) +
         theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"),
            legend.justification = c(0, 1))
      
      return(ridges)
}
