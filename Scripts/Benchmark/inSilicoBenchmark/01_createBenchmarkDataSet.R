#'---
#' title: Create Benchmarksets - Inject insilico events
#' author: Christian Mertes
#' wb:
#'  input:
#'   - datain: '`sm config["DATADIR"] + "/savedObjects/{dataset}/pvalue_psiSite.h5"`'
#'  output:
#'   - countsS: '`sm config["DATADIR"] + "/savedObjects/inSilico_n{nsamples}_{dataset}/rawCountsSS.h5"`'
#'   - pvalsS:  '`sm config["DATADIR"] + "/savedObjects/inSilico_n{nsamples}_{dataset}/pvalue_psiSite.h5"`'
#'   - wBhtml: 'Output/html/Benchmark/inSilicoBenchmark/inSilico_n{nsamples}_{dataset}.html'
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")

#'
#' # Dataset
#+ echo=TRUE
set.seed(42)
snakemake@wildcards$dataset
ns     <- snakemake@wildcards$nsamples
input  <- dirname(snakemake@input$datain)
outname <- basename(dirname(snakemake@output$pvals))


#'
#' get only expressed junctions
#'
fds <- loadFraseRDataSet(dirname(dirname(input)), basename(input))
fds <- filterExpression(fds, filter=TRUE)


#'
#' take random n samples
#'
fds[sample(1:length(fds), 50000), sample(1:dim(fds)[2], ns)]


#'
#' save first new object
#'
name(fds) <- outname
fds <- saveFraseRDataSet(fds, rewrite=TRUE)



