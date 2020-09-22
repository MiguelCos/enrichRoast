### Generate pathview visualizations of interesting KEGG pathways ----

## NOTE: Only run this script if you have executed the enrichRoast script previously and you want to explore particularly interesting KEGG pathways 
##       enriched in your data 

## Check if the required package is installed ----

packages <- c("dplyr", "magrittr")

biopackgs <- c(orgDB, "pathview")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
            install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
            if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
            
            BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
            
}

## HOW MANY KEGG GENE SETS DO YOU WANT TO VISUALIZE? ----

show_n_terms <- 10

## HOW DO YOU WANT TO SELECT THE TOP N TERMS? (BY N GENES PER TERM "NGenes" OR by bigger differences in the proportions "Difference") ----

# "NGenes" or "Difference" / "NGenes" is default

top_n_by <- "NGenes" 

## Load required packages ----

library(pathview)
library(magrittr)
library(pathview)
library(dplyr)
library(stringr)

### SCRIPT EXECUTION ----

roastOutput <- roast_result$roastOutput

topNterms <- dplyr::select(roastOutput,
                           NGenes, Direction, PropUp, PropDown, CategoryTerm,
                           FDR, PValue) %>%
            dplyr::mutate(DiffProp = abs(PropUp - PropDown),
                          PropDown = -PropDown,
                          FDR = round(FDR, 4),
                          PValue = round(PValue, 4)) %>%
            dplyr::top_n(n = show_n_terms,
                         wt = if(top_n_by == "Difference"){DiffProp}
                         else if(top_n_by == "NGenes"){NGenes})

log2fcs <- roast_result$log2FCs %>% filter(CategoryTerm %in% topNterms$CategoryTerm)

expr_data <- log2fcs$log2FC
names(expr_data) <- log2fcs$ENTREZID

keggpaths <- log2fcs$CategoryID %>% str_sub(string = ., start = 4, end = 8) %>% unique()

for (i in 1:show_n_terms){
            pathview(gene.data = expr_data, pathway.id = keggpaths[i],
                     species = organism, kegg.dir = ".", kegg.native = TRUE,same.layer=TRUE)
}

if (dir.exists("KEGG_pathways") == FALSE){dir.create("KEGG_pathways")}

kegg_files <- list.files()
kegg_files <- kegg_files[c(str_which(kegg_files, pattern = ".png$"), 
                str_which(kegg_files, pattern = ".xml$"))]

file.copy(kegg_files, "KEGG_pathways")
unlink(kegg_files)
