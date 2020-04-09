## Stand-alone script for running limma::roast starting from a typical Limma input ----
## enrichRoast script. v 0.1 . Miguel Cosenza 24.03.2020 

## PLEASE RUN THE NEXT LINES OF CODE TO CORROBORATE IF YOU HAVE INSTALLED THE REQUIRED PACKAGES ----
# Note: If some installation is needed, it could take a few minutes to finish.

### Install required packages if necessary

packages <- c("dplyr", "here", "stringr", "tidyr", "ggplot2")

biopackgs <- c(orgDB, "limma", "reactome.db", "clusterProfiler",
               "msigdbr", "KEGGREST", "AnnotationDbi", "GO.db")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
        install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        
        BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
        
}

## PLEASE MODIFY THE INPUT IN THE SCRIPT BY ANSWERING THE QUESIONS

## PLEASE GIVE A CODE TO IDENTIFY YOUR DATA ----

datasetcode <- "no-neoadjuvant-treatment"

## *-1. WHICH DATABASE WOULD YOU LIKE TO EXPLORE? (one of "GO", "KEGG", "REACTOME" or "MSIGDB") ----

enrichFunc <- "REACTOME"

## *-2. ORGANISM DATABASE (Please input the name of the Bioconductor org.db you need: i.e. "org.Hs.eg.db" for human) ----

orgDB <- "org.Hs.eg.db"
species <- "Homo sapiens" # this can be any resulting from calling msigdbr::msigdbr_show_species()

## *-3. WHAT ID TYPE ARE YOU USING? (i.e. "SYMBOL", "UNIPROT", "ENTREZID") ----

geneIDtype <- "UNIPROT"

## *-4. MINIMUM AND MAXIMUM SIZE OF THE GENE SETS YOU WANT TO TEST ----
# This depends if you want to look at general or specific terms.
# For more general enrichment, set your set sizes >  100; 
# If you want evaluate very specific pathways, it might be a better idea to set them < 100
# Example: for specific pathway enrichment:
# You will look if your proteins only appear within pathways with between 5 and 80 components.
# minSetSize = 5
# maxSetSize = 80 

minSetSize = 15
maxSetSize = 600

## *-5. P-VALUE CUTOFF AFTER FDR CONTROL TO CONSIDER A GENE SET AS ENRICHED AND NUMBER OF ROTATIONS (set to 999 for exploration and 9999 for final p-value) ----

pvalueCutoff <- 1
n_rotations = 999

## *-6 EXPERIMENTAL DESIGN ----

### Define experimental design ####

condition1 <- 4 # number of samples associated to the first condition (treatment, stage, patient, etc...)
condition2 <- 7 # number of samples associated to the second condition 

Conditions <- c("less18", "more24") # condition1, condistion2

## *-7. VISUALIZATION PARAMTERS: ----

### * 7.1. PROPORTIONS PLOT ---

#### * 7.1.1 HOW MANY ENRICHED TERMNS DO YOU WANT TO PLOT?
show_n_termsprop <- 30 # how many enriched terms do you want to plot?

#### * 7.1.2 VISUALIZE COLOR-CODING FOR "FDR" OR "PVALUE"  

# Note: visualize with "PValue" is recomended when you set up FDR cutoff to 1 because of heterogeneos data
colorbyprop <- "FDR" 

### * 7.2. RIDGELINE DENSITY PLOTS ---

#### * 7.2.1 HOW MANY ENRICHED TERMNS DO YOU WANT TO PLOT?
show_n_termsdens <- 30 # how many enriched terms do you want to plot?

#### * 7.2.2 VISUALIZE COLOR-CODING FOR "FDR" OR "PVALUE"  

# Note: visualize with "PValue" is recomended when you set up FDR cutoff to 1 because of heterogeneos data
colorbydens <- "FDR" 

## *-8. OPTIONAL PARAMETERS: FILL THESE UP DEPENDING ON WHAT YOU CHOOSE IN SECTION *-1. ----

# Note: leave at NULL if not required

### * 8.1. IF "GO" ENRICHMENT WILL BE PERFORMED ----

#### * 8.1.1. WHICH GO ONTOLOGY YOU DO WANT TO EXPLORE? (one of: "MF", "CC" or "BP") ----

ontology = NULL

#### * 8.1.2. DO YOU WANT TO REMOVE REDUNDANT GO TERMS?

simplify <- NULL
cutoff <- NULL # how similar should be two GO terms to be considered redundant? (0.7 is recommended)
by = NULL # if two terms are equally similar, which condition you want to use to select between them ("FDR")

### * 8.2. IF "REACTOME" OR "KEGG" ENRICHMENT WILL BE PERFORMED ---- 

#### * 8.2.1 DO YOU WANT TO USE A EXCLUSION LIST TO REMOVE NON-INTERSTING TERMS FROM THE FINAL OUTPUT?

exclusionList <- TRUE

### * 8.3 IF "KEGG" ENRICHMENT WILL BE PERFORMED ----

#### * 8.3.2 ORGANISM NAME IN KEGG TERMS (i.e. human = "hsa"; mouse = "mmu"; ...)

organism <- NULL

### * 8.4 IF 'MSIGDB' ENRICHMENT WILL BE PERFORMED ----

#### * 8.4.1 SET SPECIFIC DATABASE PARAMETERS
# For example: if you want to check for Matrisome/ECM components:
# category = "C2"
# subcategory = "CP"
# specific_category = "NABA"

category = NULL # Any of the main categories presented here: https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
subcategory = NULL # Any subcategory within the main categories presented in the link above (i.e. "REACTOME", "BIOCARTA", "PID"...)
specific_category = NULL  # i.e. "NABA"... A string that can be used to subset your categories.


# SCRIPT EXECUTION ----

# Load packages  ----

library(dplyr)
library(ggplot2)

## Define experimental design ####

experiment <- c(rep(0,condition1),rep(1,condition2))
design <- model.matrix(~experiment)
Paired <- FALSE
## Load Limma input ####

tabular_data <- read.delim(file = here::here("Data/input_limma.txt"), header = TRUE, stringsAsFactors = FALSE,
           sep = "\t")

## Run roast function ----

if (enrichFunc == "GO"){
            source(file = "R/roastGO.R")
            
            
            roast_result <- roastGO(data = tabular_data,
                                    geneIDtype = geneIDtype,
                                    ontology = ontology,
                                    orgDB = orgDB, 
                                    design = design,
                                    n_rotations = n_rotations,
                                    minSetSize = minSetSize,
                                    maxSetSize = maxSetSize,
                                    pvalueCutoff = pvalueCutoff,
                                    Paired = Paired)
            
            if (simplify == TRUE){
                        source("R/simplifyGO.R")
                        roast_result <- simplifyGO(roast_result, cutoff = cutoff, by = by)
            }
            
            
} else if (enrichFunc == "REACTOME"){
            source(file = "R/roastReactome.R")
            
            roast_result <- roastReactome(data = tabular_data,
                                          geneIDtype = geneIDtype,
                                          orgDB = orgDB, 
                                          design = design,
                                          n_rotations = n_rotations,
                                          minSetSize = minSetSize,
                                          maxSetSize = maxSetSize,
                                          pvalueCutoff = pvalueCutoff,
                                          exclusionList = exclusionList,
                                          species = species,
                                          Paired = Paired)
} else if (enrichFunc == "KEGG"){
            source(file = "R/roastKEGG.R")
            
            roast_result <- roastKEGG(data = tabular_data,
                                      geneIDtype = geneIDtype,orgDB = orgDB,
                                      design = design,
                                      n_rotations = n_rotations,
                                      minSetSize = minSetSize,
                                      maxSetSize = maxSetSize,
                                      pvalueCutoff = pvalueCutoff,
                                      exclusionList = exclusionList,
                                      organism = organism,
                                      
                                      Paired = Paired)
} else if (enrichFunc == "MSIGDB"){
            source(file = "R/roastMSigDB.R")
            
            roast_result <- roastMSigDB(data = tabular_data,
                                        geneIDtype = geneIDtype,
                                        orgDB = orgDB,
                                        design = design,
                                        n_rotations = n_rotations,
                                        minSetSize = minSetSize,
                                        maxSetSize = maxSetSize,
                                        pvalueCutoff = pvalueCutoff,
                                        species = species,
                                        category = category, 
                                        subcategory = subcategory,
                                        specific_category = specific_category)
}

## Visualization ----

# Proportions-change plot ----
source("R/propChangePlot.R")

prochangeplot <- propChangePlot(roast_result,
                                show_n_terms = show_n_termsprop,
                                colorby = colorbyprop)

# Ridge-line density plots -----

source("R/ridgleplotRoast.R")

ridgelineroast <- ridgeplotRoast(roast_result,
                                 show_n_terms = show_n_termsdens,
                                 colorby = colorbydens)


## Generate outputs ----

### Outputs for tabular data 

if (dir.exists(here::here("Outputs")) == FALSE){
            dir.create(here::here("Outputs"))
            
            if (dir.exists(here::here("Outputs/Tabular_data")) == FALSE){
                        dir.create(here::here("Outputs/Tabular_data"))
            }
            
            if (dir.exists(here::here("Outputs/Figures")) == FALSE){
                        dir.create(here::here("Outputs/Figures"))
            }
} 

### Outputs for tabular data ----

if (enrichFunc == "GO"){
            
            write.table(x = roast_result$roastOutput,
                        file = here::here(paste0("Outputs/Tabular_data/Roast_Output","_",enrichFunc,"_",ontology,"min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$GenesPerTerm,
                        file = here::here(paste0("Outputs/Tabular_data/GenesPerEnrichTerm","_",enrichFunc,"_",ontology,"min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$log2FCs,
                        file = here::here(paste0("Outputs/Tabular_data/CombinendRoastNLimma_wFCs","_",enrichFunc,"_",ontology,"min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
} else if (enrichFunc == "REACTOME"){
            write.table(x = roast_result$roastOutput,
                        file = here::here(paste0("Outputs/Tabular_data/Roast_Output","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$GenesPerTerm,
                        file = here::here(paste0("Outputs/Tabular_data/GenesPerEnrichTerm","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$log2FCs,
                        file = here::here(paste0("Outputs/Tabular_data/CombinendRoastNLimma_wFCs","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$ExclusionList,
                        file = here::here(paste0("Outputs/Tabular_data/Blacklist_of_terms","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
} else if (enrichFunc == "KEGG"){
            
            write.table(x = roast_result$roastOutput,
                        file = here::here(paste0("Outputs/Tabular_data/Roast_Output","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$GenesPerTerm,
                        file = here::here(paste0("Outputs/Tabular_data/GenesPerEnrichTerm","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$log2FCs,
                        file = here::here(paste0("Outputs/Tabular_data/CombinendRoastNLimma_wFCs","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$exclusionList,
                        file = here::here(paste0("Outputs/Tabular_data/Blacklist_of_terms","_",enrichFunc,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
} else if (enrichFunc == "MSIGDB"){
            
            write.table(x = roast_result$roastOutput,
                        file = here::here(paste0("Outputs/Tabular_data/Roast_Output","_",enrichFunc,"_",category,"_",subcategory,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$GenesPerTerm,
                        file = here::here(paste0("Outputs/Tabular_data/GenesPerEnrichTerm","_",enrichFunc,"_",category,"_",subcategory,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
            write.table(x = roast_result$log2FCs,
                        file = here::here(paste0("Outputs/Tabular_data/CombinendRoastNLimma_wFCs","_",enrichFunc,"_",category,"_",subcategory,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tsv")),
                        sep = "\t",
                        row.names = FALSE)
            
}

### outputs for figures ----


prochttofig  <- prochangeplot + 
            labs(caption = datasetcode,
                 subtitle = paste("Positive values = Proportion of up-regulated proteins in",
                                  Conditions[2]))+
            theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 11),
                  axis.text.y = element_text(size = 11),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                  axis.title=element_text(size=13, face="bold"),
                  legend.justification = c(0, 1),
                  plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 12, face = "plain"))

prochttofig

ridgelinetofig <- ridgelineroast +
            labs(caption = datasetcode,
                 subtitle = paste("> 0 indicates positive regulation in",
                                  Conditions[2]))+
            theme(axis.text.x = element_text(hjust = 0.5, size = 11),
                  axis.text.y = element_text(size = 11),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                  axis.title=element_text(size=13, face="bold"),
                  plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 12, face = "plain"))

ridgelinetofig

# Change plot height according to the number of terms shown ----
# For the density plot: 

if (show_n_termsdens == 25){
            heightdens <- 180
            heightdensin <- heightdens/25.4
} else if (show_n_termsdens <= 15){
            heightdens <- 160
            heightdensin <- heightdens/25.4
} else if (show_n_termsdens >= 30 & show_n_termsdens <= 40){
            heightdens <- 190
            heightdensin <- heightdens/25.4
} else if (show_n_termsdens >= 41 & show_n_termsdens <= 50){
            heightdens <- 200
            heightdensin <- heightdens/25.4
} else if (show_n_termsdens >= 51 & show_n_termsdens <= 60){
            heightdens <- 210
            heightdensin <- heightdens/25.4
} else if (show_n_termsdens >= 61 & show_n_termsdens <= 70){
            heightdens <- 220
            heightdensin <- heightdens/25.4
} else if (show_n_termsdens >= 71 & show_n_termsdens <= 100){
            heightdens <- 235
            heightdensin <- heightdens/25.4
} else if (show_n_termsdens >= 101){
            heightdens <- 250
            heightdensin <- heightdens/25.4
}

# For the proportion plot: 

if (show_n_termsprop == 25){
        heightprop <- 180
        heightpropin <- heightprop/25.4
} else if (show_n_termsprop <= 15){
        heightprop <- 160
        heightpropin <- heightprop/25.4
} else if (show_n_termsprop >= 30 & show_n_termsprop <= 40){
        heightprop <- 190
        heightpropin <- heightprop/25.4
} else if (show_n_termsprop >= 41 & show_n_termsprop <= 50){
        heightprop <- 200
        heightpropin <- heightprop/25.4
} else if (show_n_termsprop >= 51 & show_n_termsprop <= 60){
        heightprop <- 210
        heightpropin <- heightprop/25.4
} else if (show_n_termsprop >= 61 & show_n_termsprop <= 70){
        heightprop <- 220
        heightpropin <- heightprop/25.4
} else if (show_n_termsprop >= 71 & show_n_termsprop <= 100){
        heightprop <- 235
        heightpropin <- heightprop/25.4
} else if (show_n_termsprop >= 101){
        heightprop <- 250
        heightpropin <- heightprop/25.4
}

# Generate Ridgeline plot ----
ggsave(filename = here::here(paste0("Outputs/Figures/Ridgeline_plot","_",enrichFunc,ontology,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tiff")),
       plot = ridgelinetofig,
       device = 'tiff',
       width = 297,
       height = heightdens,
       units = 'mm',
       dpi = 300,
       compression = "lzw")


# Generate figure for prop-change plot ----
ggsave(filename = here::here(paste0("Outputs/Figures/Prop_Change_Plot","_",enrichFunc,ontology,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tiff")),
       plot = prochttofig,
       device = 'tiff',
       width = 297,
       height = heightprop,
       units = 'mm',
       dpi = 300,
       compression = "lzw")


# Generate report ----
excl <- roast_result$exclusionMessage

rmarkdown::render(input = here::here("R/renderReport.R"),
                  output_file = here::here(paste0("Outputs/Analysis_report","_",enrichFunc,ontology,category,subcategory,specific_category,
                                                  "_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".html")))


