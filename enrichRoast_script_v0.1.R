## Stand-alone script for running limma::roast from a Max Quant output ----
## enrichRoast script. v 0.1 . Miguel Cosenza 24.03.2020 

## PLEASE MODIFY THE INPUT IN THE SCRIPT BY ANSWERING THE QUESIONS

## PLEASE GIVE A CODE TO IDENTIFY YOUR DATA ----

datasetcode <- "MFA254-266"

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
# If you want evaluate very specific pathways, it might be a better idea te set them < 100
# Example: for specific pathway enrichment:
# You will look if your proteins only appear within pathways with between 5 and 80 components.
# minSetSize = 5
# maxSetSize = 80 

minSetSize = 100
maxSetSize = 500

## *-5. P-VALUE CUTOFF AFTER FDR CONTROL TO CONSIDER A GENE SET AS ENRICHED AND NUMBER OF ROTATIONS (set to 999 for exploration and 9999 for final p-value) ----

pvalueCutoff <- 0.9
n_rotations = 9999

## *-6 EXPERIMENTAL DESIGN ----

Conditions <- c("less18m", "more24m") # character vector with the name of the two conditions to compare as they appear in the annotation file for MSstats 
Paired <- FALSE # if paired set to TRUE

### * 6.1. DID YOU USE "MSstats" OR "MSstatsTMT"?

MSstatsType <- "MSstats"

## *-7. VISUALIZATION PARAMTERS: ----

### * 7.1. PROPORTIONS PLOT ---

#### * 7.1.1 HOW MANY ENRICHED TERMNS DO YOU WANT TO PLOT?
show_n_terms <- 30 # how many enriched terms do you want to plot?

#### * 7.1.2 VISUALIZE COLOR-CODING FOR "FDR" OR "PVALUE"  

# Note: visualize with "PValue" is recomended when you set up FDR cutoff to 1 because of heterogeneos data
colorby <- "PValue" 

### * 7.2. RIDGELINE DENSITY PLOTS ---

#### * 7.2.1 HOW MANY ENRICHED TERMNS DO YOU WANT TO PLOT?
show_n_terms <- 30 # how many enriched terms do you want to plot?

#### * 7.2.2 VISUALIZE COLOR-CODING FOR "FDR" OR "PVALUE"  

# Note: visualize with "PValue" is recomended when you set up FDR cutoff to 1 because of heterogeneos data
colorby <- "PValue" 

## *-8. OPTIONAL PARAMETERS: FILL THESE UP DEPENDING ON WHAT YOU CHOOSE IN SECTION *-1. ----

# Note: leave at NULL if not required

### * 8.1. IF "GO" ENRICHMENT WILL BE PERFORMED ----

#### * 8.1.1. WHICH GO ONTOLOGY YOU DO WANT TO EXPLORE? (one of: "MF", "CC" or "BP") ----

ontology = NULL

#### * 8.1.2. DO YOU WANT TO REMOVE REDUNDANT GO TERMS?

simplify <- NULL
cutoff <- NULL # how similar should be two GO terms to be considered redundant?
by = NULL # if two terms are equally similar, which condition you want to use to select between them

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


# EXECUTION OF THE SCRIPT; DON'T MODIFY ANYTHING FROM NOW ON ----  

# Load packages ----

library(dplyr)
#library(stringr)
#library(ggplot2)

# Load MaxQuant data and MSstats input ----
source(here::here("R/MSstatstoRoast.R"))

## Load the MSstats Processed R object if available. If not, load MaxQuant data and MSstats inputs and run MSstats -----
if(file.exists("Data/ProcessedMSstatsData_to_ROAST.Rda") | file.exists("Data/ProcessedMSstatsData_to_ROAST.rds") == TRUE){
            
            if(file.exists("Data/ProcessedMSstatsData_to_ROAST.Rda") == TRUE){
                        loadRData <- function(fileName){
                                    load(fileName)
                                    get(ls()[ls() != "fileName"])
                        }
                        
                        suppressWarnings(
                                    processed_msstats <- loadRData("Data/ProcessedMSstatsData_to_ROAST.Rda"))           
                        
            } else if (file.exists("Data/ProcessedMSstatsData_to_ROAST.rds") == TRUE){
                        
                        suppressWarnings(
                                    processed_msstats <-  readRDS(here::here("Data/ProcessedMSstatsData_to_ROAST.rds")))
                        
                        
            }
            
            # Load annotation data for MSstats and creation of the design output. 
            
            annotation <- read.delim(file = here::here("Data/annotation_to_ROAST.tsv"),
                                     stringsAsFactors = FALSE,
                                     sep = "\t")
            
            # Run MsstatsToRoast function ----
            roast_input <- MSstatsToRoast(data = processed_msstats,
                                          MSstatsAnnotation = annotation,
                                          MSstatsType = MSstatsType,
                                          Paired = Paired,
                                          Conditions = Conditions)
            
} else {
            # Load raw data for processing ----
            
            annotation <- read.delim(file = here::here("Data/annotation.tsv"),
                                     stringsAsFactors = FALSE,
                                     sep = "\t")
            
            proteinGroups <- read.delim(file = here::here("Data/proteinGroups.txt"),
                                        stringsAsFactors = FALSE,
                                        sep = "\t")
            
            evidence <- read.delim(file = here::here("Data/evidence.txt"),
                                   stringsAsFactors = FALSE,
                                   sep = "\t")
            
            if(MSstatsType == "MSstatsTMT"){
                        
            ## Run MSstatsTMT formating protein summarization ----
                        
                        processed_msstats <- MaxQtoMSstatsTMTFormat(evidence = evidence,
                                                              proteinGroups = proteinGroups,
                                                              annotation = annotation,
                                                              which.proteinid = "Leading.razor.protein",
                                                              rmProt_Only.identified.by.site = TRUE)
                        
                        processed_msstats <- proteinSummarization(processed_msstats)
                        
            } else if (MSstatsType == "MSstats"){
            ## Run MSstats formating and normalization ----
                        
                        processed_msstats <- MaxQtoMSstatsFormat(evidence = evidence,
                                                           annotation = annotation,
                                                           proteinGroups = proteinGroups,
                                                           proteinID = "Leading.razor.protein",
                                                           useUniquePeptide = TRUE)
                        
                        processed_msstats <- dataProcess(processed_msstats)
                        
                        
            }
            
            saveRDS(object = processed_msstats,
                    file = here::here("Data/ProcessedMSstatsData.rds"))
            

            stop("Please corroborate that your input only contains two experimental conditions to compare.\n
                 If so, change the names of your input files and run the script again:\n
                 annotation.tsv -> annotation_to_ROAST.tsv\n
                 ProcessedMSstatsData.rds -> ProcessedMSstatsData_to_ROAST.rds")
            
}

## Run roast function ----

if (enrichFunc == "GO"){
            source(file = "R/roastGO.R")
            
            
            roast_result <- roastGO(data = roast_input$tabular_data,
                                    geneIDtype = geneIDtype,
                                    ontology = ontology,
                                    orgDB = orgDB, 
                                    design = roast_input$design,
                                    n_rotations = n_rotations,
                                    minSetSize = minSetSize,
                                    maxSetSize = maxSetSize,
                                    pvalueCutoff = pvalueCutoff)

            if (simplify == TRUE){
                        source("R/simplifyGO.R")
                        roast_result <- simplifyGO(roast_result, cutoff = cutoff, by = by)
            }
            

} else if (enrichFunc == "REACTOME"){
            source(file = "R/roastReactome.R")
            
            roast_result <- roastReactome(data = roast_input$tabular_data,
                                          geneIDtype = geneIDtype,
                                          orgDB = orgDB, 
                                          design = roast_input$design,
                                          n_rotations = n_rotations,
                                          minSetSize = minSetSize,
                                          maxSetSize = maxSetSize,
                                          pvalueCutoff = pvalueCutoff,
                                          exclusionList = exclusionList,
                                          species = species)
} else if (enrichFunc == "KEGG"){
            source(file = "R/roastKEGG.R")
            
            roast_result <- roastKEGG(data = roast_input$tabular_data,
                                      geneIDtype = geneIDtype,orgDB = orgDB,
                                      design = roast_input$design,
                                      n_rotations = n_rotations,
                                      minSetSize = minSetSize,
                                      maxSetSize = maxSetSize,
                                      pvalueCutoff = pvalueCutoff,
                                      exclusionList = exclusionList,
                                      organism = organism)
} else if (enrichFunc == "MSIGDB"){
            source(file = "R/roastMSigDB.R")
            
            roast_result <- roastMSigDB(data = roast_input$tabular_data,
                                        geneIDtype = geneIDtype,orgDB = orgDB,
                                        design = roast_input$design,
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
                                show_n_terms = show_n_terms,
                                colorby = colorby)

# Ridge-line density plots -----

source("R/ridgleplotRoast.R")

ridgelineroast <- ridgeplotRoast(roast_result,
                                 show_n_terms = show_n_terms,
                                 colorby = "PValue")


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
                labs(caption = datasetcode)+
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
                labs(caption = datasetcode)+
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

if (show_n_terms == 25){
        height <- 180
        heightin <- height/25.4
} else if (show_n_terms <= 15){
        height <- 160
        heightin <- height/25.4
} else if (show_n_terms >= 30 & show_n_terms <= 40){
        height <- 190
        heightin <- height/25.4
} else if (show_n_terms >= 41 & show_n_terms <= 50){
        height <- 200
        heightin <- height/25.4
} else if (show_n_terms >= 51 & show_n_terms <= 60){
        height <- 210
        heightin <- height/25.4
} else if (show_n_terms >= 61 & show_n_terms <= 70){
        height <- 220
        heightin <- height/25.4
} else if (show_n_terms >= 71 & show_n_terms <= 100){
        height <- 235
        heightin <- height/25.4
} else if (show_n_terms >= 101){
        height <- 250
        heightin <- height/25.4
}

# Generate figure for prop-change plot ----
ggsave(filename = here::here(paste0("Outputs/Figures/Prop_Change_Plot","_",enrichFunc,ontology,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tiff")),
       plot = prochttofig,
       device = 'tiff',
       width = 297,
       height = height,
       units = 'mm',
       dpi = 300,
       compression = "lzw")

# Generate Ridgeline plot ----
ggsave(filename = here::here(paste0("Outputs/Figures/Ridgeline_plot","_",enrichFunc,ontology,"_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".tiff")),
       plot = ridgelinetofig,
       device = 'tiff',
       width = 297,
       height = height,
       units = 'mm',
       dpi = 300,
       compression = "lzw")

# Generate report ----
excl <- roast_result$ExclusionMessage

rmarkdown::render(input = here::here("R/renderReport.R"),
                  output_file = here::here(paste0("Outputs/Analysis_report","_",enrichFunc,ontology,category,subcategory,specific_category,
                                                  "_","min",minSetSize,"max",maxSetSize,"_","pValueCutoff",pvalueCutoff,".html")))

