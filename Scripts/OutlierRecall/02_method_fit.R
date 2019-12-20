#'---
#' title: Fit FraseR AE, PCA, PEER or BB on the data with injected outliers
#' author: Ines Scheller
#' wb:
#'  params:
#'  input:
#'   - inHtml: 'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/outlierInjection.html'
#'   - trueOut: '`sm config["DATADIR"] + "/datasets/inject_{delta}/{psiType}/savedObjects/{dataset}/trueOutliers_{psiType}.h5"`'
#'  output:
#'   - wBhtml:  'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_fit.html' 
#'  type: noindex
#'  threads: 20
#'---

#+ source main config
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", psiType="psi5", delta="uniformDistr",
            method="PCA")
    parseWBHeader2("Scripts/OutlierRecall/method_fit.R", 
            rerun=TRUE, wildcards=wildcards)
}

#+ input
fdsFile    <- snakemake@input$trueOut
dataset    <- snakemake@wildcards$dataset
workingDir <- dirname(dirname(dirname(fdsFile)))
psiType    <- snakemake@wildcards$psiType
fitMethod  <- snakemake@wildcards$method
pmethod    <- CONFIG$FDR_METHOD

bpWorkers  <- min(bpworkers(), snakemake@threads)
bpThreads  <- bpWorkers
bpProgress <- FALSE

#' ## Load dataset
#+ load fds
dataset
workingDir
fds     <- loadFraseRDataSet(workingDir, dataset)

#' ## Fit parameters
#' get fit params
fitMethod
psiType

#'
#' ## Use Method (AE, PCA or PEER) to controll for confounders
#'

#+ create output fds
fds_out <- fds
workingDir(fds_out) <- paste0(workingDir(fds), "/", fitMethod)

#+ save output fds to new directory
fds_out <- saveFraseRDataSet(fds_out)

#+ set bpparams
BPPARAM <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
register(BPPARAM)

#+ get dimension of hidden space 
#' q guess for real datasets or use best q from FraseR fit
q <- bestQ(fds_out, psiType)

#+ run fit method
currentType(fds_out) <- psiType
fds_out <- fit(fds_out, correction=fitMethod, q=q, type=psiType, 
        BPPARAM=BPPARAM, verbose=TRUE, iterations=15)

#+ save output fds
fds_out <- saveFraseRDataSet(fds_out)

