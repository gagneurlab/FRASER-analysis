#'---
#' title: All paper figure (heatmaps)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut:  '`sm expand("Output/html/PaperFigures/heatmaps/heatmapFigure_{psiType}.html", psiType=config["psiTypes"])`'
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
psiTypes          <- snakemake@config$psiTypes

#+ echo=FALSE, results="asis"
cat("<h1>Heatmap figures per psiType:</h1><p>")
devNull <- sapply(psiTypes, function(type){
            cat(paste0(
                "</br>", "<a href='PaperFigures/heatmaps/heatmapFigure", "_", type, ".html'>",
                type, " heatmap figure</a>"
            ))
})
