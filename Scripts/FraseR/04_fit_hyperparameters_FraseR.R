#'---
#' title: Hyper parameter optimization
#' author: Christian Mertes
#' wb:
#'  params:
#'   - progress: false
#'  input:
#'   - dPsiSS: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'   - html:   '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_filterExpression.html"`'
#'  output:
#'   - dPsiSS: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/delta_psiSite.h5"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_hyper_parameter_optimization.html"`'
#'  threads: 20
#'  type: noindex
#'---
##
## TODO:
##   Add a link to the fraser object to have a proper chain of events in wbuild
##

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", method="FraseR")
    parseWBHeader2("Scripts/FraseR/04_fit_hyperparameters_FraseR.R", 
            wildcards=wildcards, rerun=TRUE)
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
method     <- snakemake@wildcards$method
name       <- paste0(dataset, "__", method)
workingDir <- file.path(snakemake@config$DATADIR, "datasets")
bpWorkers  <- bpMaxWorkers(snakemake@threads)
bpThreads  <- 100
bpProgress <- as.logical(snakemake@params$progress)


#'
#' # Load PSI data
#+ echo=TRUE
dataset
name

#+ echo=FALSE
fds <- loadFraserDataSet(dir=workingDir, name=dataset)

#' Drop unneded columns/assays so we dont carry them along
metadata(fds) <- list()
dropAssays <- c("predictedMeans", "pvalues", "pajd", "zScores")
for(da in dropAssays){
    dropThem <- assayNames(fds)[
            grepl(paste0("^", da, ".*"), assayNames(fds), perl=TRUE)]
    for(i in dropThem){
        assay(fds, i) <- NULL
    }
}
keepMcols <- c("startID", "endID", "maxCount", "quantileValue5",
        "quantileValue3", "maxDPsi3", "maxDPsi5", "passed", "hgnc_symbol",
        "spliceSiteID", "type")
for(type in c("j", "ss")){
    cn <- colnames(mcols(fds, type=type))
    mcols(fds, type=type) <- mcols(fds, type=type)[cn %in% keepMcols]
}
gc()

#' Resave the data to new folder
fds <- saveFraserDataSet(fds, dir=workingDir, name=name, rewrite=TRUE)
BPPARAM <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
gc()

#'
#' # Run hyper parameterization
#'
for(type in psiTypes){
    message(date(), ": ", type)
    q_param     <- c(seq(2, 20, by=2), seq(18, 40, by=5), 50, 70)
    q_param     <- sort(unique(pmin(q_param, ncol(fds))))
    noise_param <- c(0.5) # c(0.5, 1.5)
    minDeltaPsi <- 0.1

    fds <- optimHyperParams(fds, type=type, implementation=method, iterations=4,
            q_param=q_param, noise_param=noise_param, BPPARAM=BPPARAM,
            internalThreads=1)
    fds <- saveFraserDataSet(fds)
}

#'
#' FraseR object
#'
fds <- saveFraserDataSet(fds)
fds

