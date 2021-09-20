### propChangePlot function ----
#show_n_terms = show_n_termsprop
#colorby = colorbyprop
#top_n_by = "Difference"
#at_least_n_genes = at_least_n_genesprop
#roastresult <- roast_result

propChangePlot <- function(roastresult,
                           show_n_terms = 25,
                           colorby = "FDR",
                           top_n_by = "NGenes", # one of "Difference", "NGenes", "PValue" or "FDR"
                           at_least_n_genes = 0){ 
   
   # Load required packages ----
   
   require(ggplot2) || stop("Package ggplot2 is required")
   require(dplyr) || stop("Package dplyr is required")
   require(tidyr) || stop("Package tidyr is required")
   require(forcats) || stop("Package forcats is required")
   require(stringr) || stop("Package stringr is required")
   
   # Wrangle/prep data -----
   
   roastOutput <- roastresult$roastOutput %>% filter(NGenes >= at_least_n_genes)
   
   toproplot <- dplyr::select(roastOutput,
                              NGenes, Direction, PropUp, PropDown, CategoryTerm,
                              FDR, PValue) %>%
      #dplyr::top_n(n = show_n_terms,
      #            wt = NGenes) %>%
      dplyr::mutate(DiffProp = abs(PropUp - PropDown),
                    PropDown = -PropDown) %>%
      dplyr::top_n(n = show_n_terms,
                   wt = if(top_n_by == "Difference"){DiffProp}
                   else if(top_n_by == "NGenes"){NGenes} 
                   else if(top_n_by == "PValue"){-PValue}
                   else if(top_n_by == "FDR"){-FDR}
      ) %>%
      #dplyr::mutate(FDR = round(FDR, 4),
      #PValue = round(PValue, 4)) %>%
      tidyr::pivot_longer(cols = c(PropDown, PropUp),
                          names_to = "PropDirection",
                          values_to = "Proportion") %>%
      #dplyr::top_n(n = show_n_terms,
      #             wt = NGenes)
      dplyr::group_by(CategoryTerm, NGenes) %>% ungroup() %>%
      mutate(CategoryTerm = stringr::str_wrap(CategoryTerm, width = 30))
   
   zero_range <- function(x) {
      if (length(x) == 1) return(TRUE)
      x <- range(x) / mean(x)
      isTRUE(all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5))
   }
   
   pvals <- dplyr::pull(toproplot, eval(as.name(colorby)))
   
   #if(isEmpty(pvals)){stop("Error: no terms to plot")}
   
   if (zero_range(pvals) == TRUE){
      maxpval <- max(pvals)
      limits <- c(0,maxpval)
      breaks <- round(seq(0, maxpval, length = 7), 4)
   } else if(zero_range(pvals) == FALSE){
      maxpval <- max(pvals)
      minpval <- min(pvals)
      limits <- c(minpval,maxpval)
      breaks <- round(seq(minpval, maxpval, length = 7), 4)
   }
   
   # Plot ----
   
   proplot <- ggplot(data = toproplot,
                     aes(x=fct_reorder(CategoryTerm, Proportion), y=Proportion))+
      coord_flip() +
      geom_line(aes(group = CategoryTerm))+
      geom_hline(yintercept = 0, color = "red")+
      geom_point(aes(color=eval(as.name(colorby)), size = NGenes))+
      scale_color_gradient(high = "#0fabbc",
                           low = '#fa163f',
                           guide = guide_colourbar(reverse = TRUE),
                           limits = limits,
                           breaks = breaks,
                           name = colorby)+
      geom_text(data = dplyr::filter(toproplot, PropDirection == "PropUp"), 
                aes(label=round(Proportion,2)),
                hjust = -0.85)+
      geom_text(data = dplyr::filter(toproplot, PropDirection == "PropDown"), 
                aes(label=round(Proportion,2)),
                hjust = +1.85)+
      scale_y_continuous(expand=c(0.2,0), limits=c(-1, 1))+
      labs(title = "Proportion of Up- or Down-regulated Proteins by Category", 
           subtitle= "Positive values = Proportion of up-regulated proteins in category",
           x="Biological category", y = "Proportion of Proteins",
           color = colorby,
           size = "N Genes per set",
           caption = if(top_n_by == "Difference"){paste0("Showing top ",show_n_terms," terms by |ProportionUp - ProportionDown|")}
           else if(top_n_by == "NGenes"){paste0("Showing top ",show_n_terms," terms by N Genes per set")})+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"),
            legend.justification = c(0, 1))
   
   return(proplot)
   
}
