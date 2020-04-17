#' ---
#' title: "enrichROAST script report"
#' ---

#'## ROAST enrichment report for the analysis against the **`r enrichFunc`** database.  
#' 
#'### Dataset code: `r datasetcode`
#' 
#'##### After analysis, `r dim(roast_result$roastOutput)[1]` terms were enriched under the user-specified parameters.  
#'  
#'### Visualization of results:  
#' 
#'#### Differential proportion plots  
#'  
#'##### Showing the top `r show_n_termsprop` categories with the higher number of proteins identified within them. 
#' 
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE
print(prochttofigngenes)  
#'
#'##### Showing the top `r show_n_termsprop` categories with the higher value of |ProportionUp - ProportionDown|. 
#' 
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE
print(prochangeplotdiff)  
#'
#'#### Ridgeline log2(FC) density plots
#'
#'##### Showing the top `r show_n_termsdens` categories with the higher number of proteins identified within them. 
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE, message = FALSE
print(ridgelinetofigngenes)
#' 
#'##### Showing the top `r show_n_termsdens` categories with the higher value of |ProportionUp - ProportionDown|.
#' 
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE, message = FALSE
print(ridgelinetofigdiff)
#' 
#'### Organism database used: **`r orgDB`**
#'### User parameters: 
#' 
#' * geneIDtype = `r geneIDtype`
#' * minSetSize = `r minSetSize`
#' * maxSetSize = `r maxSetSize`
#' * P-value cut-off = `r pvalueCutoff`
#' * Cut-off by = `r cutoff_by` 
#' 
#'### Experimental design: 
#' 
#' * Conditions = `r Conditions`
#' * Paired = `r Paired`  
#' 
#'### Specific parameters by roastFunction used:  
#'
#'#### roastGO parameters set:
#' 
#' * Ontology tested: `r ontology`
#' * Redundant GO terms removed?: `r simplify`
#' * Similarity cut-off: `r cutoff`
#' * Reduntant terms filtered by: `r by`
#' 
#' #### roastKEGG and roastReactome parameters set: 
#' 
#' * Blacklist: `r exclusionList`
#' * Exclusion message: `r excl`
#' * Organism: `r organism`
#' * Species: `r species`
#' 
#' #### roastMSigDB parameters set: 
#'  
#' * Category: `r category`
#' * Subcategory: `r subcategory`
#' * Specific category: `r specific_category`

