# enrichRoast
## Scripts and functions to run limma::roast enrichment analyses  

This repo contains a series of functions and scripts to facilitate the execution of limma::roast (rotation gene set tests) as part of an MSstats workflow and generate tabular outputs and visualizations. 

#### Currently, it is possible to: 

* Run roast against four gene set databases (GO, Reactome, KEGG and MSigDB),
* Select specific subcategories within each database (i.e. `MF` ontology for GO, Hallmark subcategory for MSigDB),
* Setup a p-value cutoff for the tabular output,
* Select the size of the gene sets to be evaluated,
* Filter-out redundant (Semantically similar) GO terms or non-interesting ones (Blacklist) from Reactome or KEGG.

## General workflow (starting from a Max Quant output): 

1. Download this repo into your computer. You should have two scripts and two folders.
2. Initialize the repo in your local computer as an RStudio project.
3. Add the required input data files into the `Data` folder:  
    + `proteinGroups.txt` file from the MaxQuant output.  
    + `evidence.txt` file from the MaxQuant output.
    + `annotation.tsv` file with the MSstats/MSstatsTMT experimental design (You can consult their documentation here: https://msstats.org/). 
    
4. In RStudio, open `prepMSstats_format_script_vx.x.R`
5. Modify the parameters in lines 14 (type of MSstats used) and 19 (the names of the two experimental conditions to compare as they appear in the `annotation.tsv` file.  
6. Click `Source`.
7. The script will run MSstats/MSstatsTMT for summarization and normalization and will produce three new files:
    + `annotation_to_ROAST.tsv`: a reduced version of your initial annotation file with only the selected conditions for comparison.  
    + `ProcessedMSstatsData.rds`: an R object file containing the MSstats/MSstatsTMT processed data. 
    + `ProcessedMSstatsData_to_ROAST.rds`: an R object file containing the MSstats/MSstatsTMT processed data, reduced to only contain information related to the two experimental conditions you want to compare.  
8. Open the main script `enrichRoast_script_vx.x.R`.
9. Modify the parameters between lines 6 and 100, by answering the questions in the comments. 
10. Click `Source`.
11. A new `Outputs` folder should have been generated with figures, tabular outputs and an HTML report. 

### Notes (things to be aware of, implementations required, etc.):

* The `R` folder contains all the functions required to run this analysis in a modular way. The idea is to create a package from them so people in the lab could install it and implement these functions into their own scripts/workflows as they require it. Everyone is encouraged to test the functions individually too. I will soon write individual documentation for each one.
* The processing and summarization with MSstats/MSstatsTMT takes time, this is why the script is generating an R object (`ProcessedMSstatsData.rds`/`ProcessedMSstatsData_to_ROAST.rds`) to avoid running MSstats/MSstatsTMT every time we need to run roast. If you have already processed data through MSstats/MSstatsTMT you can save it using the function `saveRDS()` with the name `ProcessedMSstatsData_to_ROAST.rds`if it already contains only two experimental conditions or with the name `ProcessedMSstatsData.rds` if you would require to use the `prepMSstats_format_script_vx.x.R` for selecting the desired two conditions to compare.
* If you have a `ProcessedMSstatsData_to_ROAST.rds`/`ProcessedMSstatsData.rds` file, you can include it in the `Data` folder as an input to avoid running MSstats/MSstatsTMT every time you want to do a roast enrichment using these scripts/functions.

