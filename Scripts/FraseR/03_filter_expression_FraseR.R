#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  threads: 1
#'  input:
#'   - psiSS:  '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'   - dPsiSS: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/delta_psiSite.h5"`'
#'  output:
#'   - dPsiSS: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset,...[^-].*}/fds-object.RDS"`'
#'   - wBhtml: 'Output/html/FraseR/{dataset}_filterExpression.html'
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=12, fig.height=8)

#+ example data, echo=FALSE
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer")
    parseWBHeader2("Scripts/FraseR/03_filter_expression_FraseR.R",
            wildcards=wildcards, rerun=TRUE)
}

#+ input
dataset     <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir  <- dirname(dirname(dirname(snakemake@output$dPsiSS)))
register(MulticoreParam(15))

#'
#' # Load count data
#+ echo=TRUE
dataset

#+ load data, echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))


#'
#' Filter FraseR object based on standard values
#'
#+ run filter expression
fds <- filterExpression(fds, filter=FALSE, delayed=TRUE)
table(mcols(fds, type="j")[,"passed"])
plotFilterExpression(fds)


#+ filter bad junctions, echo=FALSE
devNull <- saveFraseRDataSet(fds)
name(fds) <- dataset
fds <- saveFraseRDataSet(fds[mcols(fds, type="j")[,"passed"]])

#'
#' Correlation of counts after filtering
#'
#+ plot correlation heatmap
plots <- lapply(psiTypes, plotCountCorHeatmap, fds=fds, logit=TRUE, topN=100000)

