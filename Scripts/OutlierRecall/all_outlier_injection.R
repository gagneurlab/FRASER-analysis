#'---
#' title: Outlier injection for all datasets
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut: '`sm expand("Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/outlierInjection.html", dataset=config["datasets"], psiType=config["psiTypes"], delta=config["inj_deltas"])`'
#' output:
#'  html_document
#'---

if(FALSE){
  snakemake <- readRDS("./tmp/snakemake.RDS")
}

#+ source config
source("./src/r/config.R")

#+ input
allOutFiles <- snakemake@input$allOut
datasets    <- snakemake@config$datasets
psiTypes    <- snakemake@config$psiTypes
deltas      <- snakemake@config$inj_deltas

#+ echo=FALSE, results="asis"
cat("<h1>Outlier Injection per dataset</h1><p>")
devNull <- sapply(datasets, function(name){
  cat(paste0("<h2>Dataset: ", name, "</h1><p>"))
  for(type in psiTypes){
    for(delta in deltas){
      cat(paste0(
        "</br>", "<a href='OutlierInjection/", name, "/", type, "/", delta, "/outlierInjection.html'>", name, " ", type, " ", delta, " outlier injection</a>"
      ))
    }
  }
  cat("</br> </p>")
})