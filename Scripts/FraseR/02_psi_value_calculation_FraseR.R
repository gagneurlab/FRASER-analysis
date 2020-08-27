#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - progress: FALSE
#'  input:
#'   - countsJ:  '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsSS: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - html:     '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_counting.html"`'
#'  output:
#'   - psiSS:  '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'   - dPsiSS: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/delta_psiSite.h5"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_psi_value_calculation.html"`'
#'  threads: 10
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")

# for interactive mode
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer")
    parseWBHeader2("./Scripts/FraseR/02_psi_value_calculation_FraseR.R",
            wildcards=wildcards, rerun=TRUE)
}

#+ input
dataset     <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir  <- dirname(dirname(dirname(snakemake@input$countsJ)))
bpWorkers   <- bpMaxWorkers(snakemake@threads)
bpThreads   <- bpWorkers*4
bpProgress  <- as.logical(snakemake@params$progress)

#'
#' # Load count data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))
BPPARAM <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
register(BPPARAM)

#'
#' Calculating PSI values
#'
fds <- calculatePSIValues(fds)

#'
#' FraseR object after PSI value calculation
#'
fds <- saveFraserDataSet(fds)
fds

