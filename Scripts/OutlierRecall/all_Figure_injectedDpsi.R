#'---
#' title: All injected dPSI vs fitted dPSI plots
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut_full:  '`sm expand("Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_injectedDpsi.html", dataset=config["datasets"], psiType=config["psiTypes"], delta=config["inj_deltas"], method=config["methods"])`'
#' output:
#'  html_document
#'---

if(FALSE){
  snakemake <- readRDS("./tmp/snakemake.RDS")
  
}

#+ source main config
source("./src/r/config.R")

#+ input
allOutFiles_full  <- snakemake@input$allOut_full
datasets          <- snakemake@config$datasets
psiTypes          <- snakemake@config$psiTypes
deltas            <- snakemake@config$inj_deltas
methods           <- snakemake@config$methods

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  cat(paste0("<h2>Dataset: ", name, "</h1><p>"))
  for(type in psiTypes){
    for(delta in deltas){
      for(method in methods){
        cat(paste0(
          "</br>", "<a href='OutlierInjection/", name, "/", type, "/", delta, "/", method, "_injectedDpsi.html'>", 
          name, " ", type, " ", delta, " ", method, " injected dPSI vs fitted dPSI plot</a>"
        ))
      }
    }
  }
  cat("</br> </p>")
})