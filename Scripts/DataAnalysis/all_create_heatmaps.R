#'---
#' title: Create heatmaps
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut: '`sm expand("Output/html/DataAnalysis/heatmaps/{dataset}_heatmaps.html", dataset=config["EnrichmentTissues"])`'
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

#+ echo=FALSE, results="asis"
cat("<h1>Heatmaps per dataset</h1><p>")
devNull <- sapply(datasets, function(name){
    cat(paste0(
        "<a href='DataAnalysis/heatmaps/", name, "_heatmaps.html'>",
        name, " heatmaps</a>",
        "</br>"))
})
cat("</p>")
