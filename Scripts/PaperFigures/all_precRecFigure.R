#'---
#' title: All paper figure (precRec)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut:  '`sm expand("Output/html/PaperFigures/precRec/{dataset}/{delta}/precRec_{outlierType}_{psiType}.html", dataset=config["datasets"], psiType=config["psiTypes"], delta=config["inj_deltas"], outlierType=config["outlier_types"])`'
#' output:
#'  html_document
#'---

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")

}

#+ source main config
source("./src/r/config.R")

#+ input
allOutFiles       <- snakemake@input$allOut
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
                    "</br>", "<a href='PaperFigures/", name, "/", delta, "/precRec_", outlierType, "_", type, ".html'>",
                    name, " ", type, " ", delta, " ", outlierType, " outliers precision-recall plot</a>"
                ))
            }
        }
    }
    cat("</br> </p>")
})
