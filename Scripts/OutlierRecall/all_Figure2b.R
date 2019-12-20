#'---
#' title: All Figure 2b plots
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut_binned:  '`sm expand("Output/html/OutlierRecall/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall.html", dataset=config["datasets"], psiType=config["psiTypes"], delta=config["inj_deltas"], outlierType=config["outlier_types"])`'
#'   - allOut_overall:  '`sm expand("Output/html/OutlierRecall/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_overall.html", dataset=config["datasets"], psiType=config["psiTypes"], delta=config["inj_deltas"], outlierType=config["outlier_types"])`'
#' output:
#'  html_document
#'---

if(FALSE){
  snakemake <- readRDS("./tmp/snakemake.RDS")
  
}

#+ source main config
source("./src/r/config.R")

#+ input
allOutFiles       <- snakemake@input$allOut_binned
datasets          <- snakemake@config$datasets
psiTypes          <- snakemake@config$psiTypes
deltas            <- snakemake@config$inj_deltas
outlierTypes      <- snakemake@config$outlier_types

#+ echo=FALSE, results="asis"
cat("<h1>Precision-recall plots per dataset</h1><p>")
devNull <- sapply(datasets, function(name){
  cat(paste0("<h2>Dataset: ", name, "</h1><p>"))
  for(type in psiTypes){
    for(delta in deltas){
      for(outlierType in outlierTypes){
        cat(paste0(
          "</br>", "<a href='OutlierRecall/", name, "/", delta, "/", type, "/", outlierType, "_outlier_recall.html'>", 
          name, " ", type, " ", delta, " ", outlierType, " outliers precision-recall plot</a>"
        ))
        cat(paste0(
          "</br>", "<a href='OutlierRecall/", name, "/", delta, "/", type, "/", outlierType, "_outlier_recall_overall.html'>", 
          name, " ", type, " ", delta, " ", outlierType, " overall precision-recall plot</a>"
        ))
      }
    }
  }
  cat("</br> </p>")
})
