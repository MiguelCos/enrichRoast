## propChangePlot function ----
propChangePlot <- function(roastresult,
                           show_n_terms = 25,
                           colorby = "FDR"){
      
      # Load required packages ----
      
      require(ggplot2) || stop("Package ggplot2 is required")
      require(dplyr) || stop("Package dplyr is required")
      require(tidyr) || stop("Package tidyr is required")
      require(forcats) || stop("Package forcats is required")
      require(stringr) || stop("Package stringr is required")

      # Wrangle/prep data -----

      roastOutput <- roastresult$roastOutput
      
      toproplot <- dplyr::select(roastOutput,
                                 NGenes, Direction, PropUp, PropDown, CategoryTerm,
                                 FDR, PValue) %>%
            dplyr::top_n(n = show_n_terms,
                         wt = NGenes) %>%
            dplyr::mutate(DiffProp = abs(PropUp - PropDown),
                          PropDown = -PropDown) %>%
            tidyr::pivot_longer(cols = c(PropDown, PropUp),
                                names_to = "PropDirection",
                                values_to = "Proportion") %>%
            dplyr::group_by(CategoryTerm, NGenes) 
      
      # Plot ----
      
      proplot <- ggplot(data = toproplot,
                        aes(x=fct_reorder(CategoryTerm, Proportion), y=Proportion))+
            coord_flip() +
            geom_line(aes(group = CategoryTerm))+
            geom_hline(yintercept = 0, color = "red")+
            #and the oints
            geom_point(aes(color=eval(as.name(colorby)), size = NGenes))+
            geom_text(data = dplyr::filter(toproplot, PropDirection == "PropUp"), 
                      aes(label=round(Proportion,2)),
                      hjust = -0.85)+
            geom_text(data = dplyr::filter(toproplot, PropDirection == "PropDown"), 
                      aes(label=round(Proportion,2)),
                      hjust = +1.85)+
            scale_y_continuous(expand=c(0.2,0), limits=c(-1, 1))+
            #scale_x_discrete(labels = function(x)  str_wrap(x, width = 40))+
            labs(title = "Proportion of Up- or Down-regulated Proteins by Category", 
                 subtitle= "Positive values = Proportion of up-regulated proteins in category",
                 x="Biological category", y = "Proportion of Proteins",
                 color = colorby,
                 size = "N Genes per set")+
            theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                  axis.title=element_text(size=10,face="bold"),
                  legend.justification = c(0, 1))
      
      return(proplot)
      
}
