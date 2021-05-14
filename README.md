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

## General workflow (starting from log-transformed and normalized data):  

1. You will need a `txt` file named `input_limma.txt`. The first column should be name `ID` and each other column should correspond to different samples/MS runs.
2. Download this repo into your computer. You should have two scripts and two folders.
3. Initialize the repo in your local computer as an RStudio project.
4. Add `input_limma.txt` into the `Data` folder.
5. Open the main script `enrichRoast_script_limma_input_vx.x.R`
6. Modify the parameters between lines 8 and 122, by answering the questions in the comments.
7. Click `Source`. 
11. A new `Outputs` folder should have been generated with figures, tabular outputs and an HTML report. 

## Visualizing your quantitative protein information on KEGG pathways:

The `visualizeKEGG.R` script would help you map your enriched KEGG pathways to the quantified proteins in your data set. This script should only be executed after the ROAST enrichment analysis (after executing any of the previous `enrichRoast_**` scripts).

1. Open the `visualizeKEGG.R` script.
2. Execute the lines between 6 and 21 to check for the proper installation and loading of the required packages.
3. Modify the input in lines 25 and 31 according to the questions in the comments.
4. Execute the rest of the script.
5. A number of `png` files with the protein to KEGG pathway mapping should have been created in your working directory.

## Potential errors:

* `in seq.default(0, maxpval, length = 7) : 'to' must be a finite number`: If this happens when executing either of the plotting functions, you probably need to adapt the `at_least_n_genes*` parameters. Probably some of the enriched terms had < `at_least_n_genes*`.

### Notes (things to be aware of, implementations required, etc.):

* The `R` folder contains all the functions required to run this analysis in a modular way. The idea is to create a package from them so people in the lab could install it and implement these functions into their own scripts/workflows as they require it. Everyone is encouraged to test the functions individually too. I will soon write individual documentation for each one.
* __For advanced users__: If you use the functions in a modular way, you can now use complex `design` matrices as an input. If you do so, please modify the `limma_coef` parameter accordingly. Normally, if you want to include an random effect (i.e. patient-matched design) extracting the `coef = 2` should suffice when doing pairwise comparisons. Including constrast matrices should also be feasible but this is still not implemented, so better stick to 2-group comparisons for now.
* The processing and summarization with MSstats/MSstatsTMT takes time, this is why the script is generating an R object (`ProcessedMSstatsData.rds`/`ProcessedMSstatsData_to_ROAST.rds`) to avoid running MSstats/MSstatsTMT every time we need to run roast. If you have already processed data through MSstats/MSstatsTMT you can save it using the function `saveRDS()` with the name `ProcessedMSstatsData_to_ROAST.rds`if it already contains only two experimental conditions or with the name `ProcessedMSstatsData.rds` if you would require to use the `prepMSstats_format_script_vx.x.R` for selecting the desired two conditions to compare.
* If you have a `ProcessedMSstatsData_to_ROAST.rds`/`ProcessedMSstatsData.rds` file, you can include it in the `Data` folder as an input to avoid running MSstats/MSstatsTMT every time you want to do a roast enrichment using these scripts/functions.

