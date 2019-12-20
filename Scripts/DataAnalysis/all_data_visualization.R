#'---
#' title: Splicing Correlation
#' author: Christian Mertes
#' wb:
#'  input:
#'   - allOut: '`sm expand("Output/html/DataAnalysis/data_viz/{dataset}_data_viz.html", dataset=config["datasets"])`'
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
datasets <- gsub("_data_viz.html$", "", basename(allOutFiles))

#+ echo=FALSE, results="asis"
cat("<h1>Data visualizations per dataset</h1><p>")
devNull <- sapply(datasets, function(name){
    cat(paste0(
        "<a href='DataAnalysis/data_viz/", name, "_data_viz.html'>",
        name, " data vizualization</a>",
        "</br>"))
})
cat("</p>")
