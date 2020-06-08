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

## Load required packages ----

library(pathview)
library(magrittr)
library(pathview)
library(dplyrW)

# * 1. RUN THE NEXT COMMAND TO visualize the enriched KEGG gene sets after enrichRoast and to extract their respective pathway ID ----

View(roast_result$log2FCs %>% dplyr::select(CategoryID, CategoryTerm) %>% dplyr::distinct())

#### * 2. WHICH KEGG PATHWAY YOU WANT TO MAP TO YOUR PROTEIN IDS? ----

# Please provide the code of the KEGG pathway that you want to explore 
# The code should be supplied without the organism prefix.
# Example: 'Lysosome' gene set for mouse = mmu04142; then, you must provide "04142"

keggpath <- "04144"

### SCRIPT EXECUTION -----

expr_data <- roast_result$log2FCs$log2FC
names(expr_data) <- roast_result$log2FCs$ENTREZID


pathview(gene.data = expr_data, pathway.id = keggpath,
                   species = organism)
