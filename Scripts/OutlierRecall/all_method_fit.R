#'---
#' title: Fit PCA, PEER and FraseR on injected data
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut: '`sm expand("Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_fit.html", dataset=config["datasets"], psiType=config["psiTypes"], delta=config["inj_deltas"], method=config["methods"])`'
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
methods     <- snakemake@config$methods
psiTypes    <- snakemake@config$psiTypes
deltas      <- snakemake@config$inj_deltas

#+ echo=FALSE, results="asis"
cat("<h1>Fit of methods (FraseR, PCA, PEER) per dataset</h1><p>")
devNull <- sapply(datasets, function(name){
  cat(paste0("<h2>Dataset: ", name, "</h1><p>"))
  for(type in psiTypes){
    for(delta in deltas){
      for(method in methods){
        cat(paste0(
          "</br>", "<a href='OutlierInjection/", name, "/", type, "/", delta, "/", method, "_fit.html'>", name, " ", type, " ", delta, " ", method, "fit</a>"
        ))
      }
    }
  }
  cat("</br> </p>")
})
