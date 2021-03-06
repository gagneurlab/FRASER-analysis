#'---
#' title: Calculate P values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - progress: false
#'  input:
#'   - inHtml: 'Output/html/FraseR/{dataset}/{method}_autoencoder_fit.html'
#'  output:
#'   - fdsout: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/padjBetaBinomial_psiSite.h5"`'
#'   - wBhtml: 'Output/html/FraseR/{dataset}/{method}_stat_calculation.html'
#'  threads: 20
#'  type: noindex
#'---


#+ echo=FALSE
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", method="FraseR")
    parseWBHeader2("./Scripts/FraseR/06_calculation_stats_AE_FraseR.R", 
            wildcards=wildcards, rerun=TRUE)
}

#+ input
dataset    <- snakemake@wildcards$dataset
method     <- snakemake@wildcards$method
name       <- paste0(dataset, "__", method)
fdsFile    <- snakemake@output$fdsout
workingDir <- dirname(dirname(dirname(fdsFile)))
bpWorkers  <- bpMaxWorkers(snakemake@threads)
bpProgress <- as.logical(snakemake@params$progress)

#'
#' # Load Zscores data
#+ echo=TRUE
dataset
method
workingDir

#+ echo=FALSE
fds <- loadFraserDataSet(dir=workingDir, name=name)
BPPARAM <- MulticoreParam(bpWorkers, progressbar=bpProgress)
register(BPPARAM)


#'
#' # Calculate stats
#'
#' ## Zscores
#
for(type in psiTypes){
    message(date(), ": running zscore for ", type, " ", method, " ", dataset)
    fds <- calculateZscore(fds, type=type)
}

#'
#' ## Pvalues
for(type in psiTypes){
    message(date(), ": running pvalue for ", type, " ", method, " ", dataset)
    fds <- calculatePvalues(fds, type=type, implementation=method, 
            distributions=c("betabinomial", "normal"))
}

#'
#' ## Adjust Pvalues
for(type in psiTypes){
    message(date(), ": running padjust for ", type, " ", method, " ", dataset)
    fds <- calculatePadjValues(fds, type=type)
}

#'
#' ## Annotate ranges
fds <- annotateRanges(fds)


#'
#' # Save results
#'
fds <- saveFraserDataSet(fds)
fds

