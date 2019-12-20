#'---
#' title: All Precision Recall
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut:  '`sm expand("Output/html/OutlierRecall/{dataset}/{delta}/{psiType}/{method}_precRec_{outlierType}.html", dataset=config["datasets"], psiType=config["psiTypes"], delta=config["inj_deltas"], method=config["methods"], outlierType=config["outlier_types"])`'
#' output:
#'  html_document
#'---

if(FALSE){
  snakemake <- readRDS("./tmp/snakemake.RDS")
  
}

#+ source main config
source("./src/r/config.R")

#+ input
allOutFiles  <- snakemake@input$allOut
datasets     <- snakemake@config$datasets
methods      <- snakemake@config$methods
psiTypes     <- snakemake@config$psiTypes
deltas       <- snakemake@config$inj_deltas
outlierTypes <- snakemake@config$outlier_types

#+ echo=FALSE, results="asis"
cat("<h1>Precision and recall for methods (FraseR, PCA, PEER) per dataset</h1><p>")
devNull <- sapply(datasets, function(name){
  cat(paste0("<h2>Dataset: ", name, "</h1><p>"))
  for(type in psiTypes){
    for(delta in deltas){
      for(method in methods){
        for(outlierType in outlierTypes){
          cat(paste0(
            "</br>", "<a href='OutlierRecall/", name, "/", delta, "/", type, "/", method, "_precRec_", outlierType, ".html'>", 
            name, " ", type, " ", delta, " ", method, " ", outlierType, " outlier recall</a>"
          ))
        }
      }
    }
  }
  cat("</br> </p>")
})
