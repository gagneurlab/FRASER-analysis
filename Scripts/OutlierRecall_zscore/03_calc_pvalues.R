#'---
#' title: Calculate p-values for FraseR AE, PCA, PEER or BB on the data with injected outliers
#' author: Ines Scheller
#' wb:
#'  input:
#'   - inHtml: 'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_fit.html' 
#'   - trueOut: '`sm config["DATADIR"] + "/datasets/inject_{delta}/{psiType}/savedObjects/{dataset}/trueOutliers_{psiType}.h5"`'
#'  output:
#'   - fdsOut: '`sm config["DATADIR"] + "/datasets_zscoreCheck/inject_{delta}/{psiType}/{method}/savedObjects/{dataset}/padjBetaBinomial_{psiType}.h5"`'
#'   - wBhtml:  'Output/html/zscoreCheck/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_pvalues.html' 
#'  type: noindex
#'  threads: 20
#'---

#+ source main config
source("./src/r/config.R")

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    
    fds_bb <- loadFraseRDataSet("Data/paperPipeline/datasets/inject_delta/psi3", "SimulationBB", TRUE)
    fds_dm <- loadFraseRDataSet("Data/paperPipeline/datasets/inject_delta/psi3", "SimulationDM", TRUE)
    
}

#+ input
fdsFile    <- snakemake@input$trueOut
dataset    <- snakemake@wildcards$dataset
psiType    <- snakemake@wildcards$psiType
fitMethod  <- snakemake@wildcards$method
workingDir <- file.path(dirname(dirname(dirname(fdsFile))), fitMethod)
newWorkingdir <- gsub("datasets/", "datasets_zscoreCheck/", workingDir)
pmethod    <- CONFIG$FDR_METHOD

bpWorkers  <- min(bpworkers(), snakemake@threads)
bpThreads  <- bpWorkers
bpProgress <- FALSE

#' ## Load dataset and copy to new dir
#+ load fds
dataset
workingDir
fds     <- loadFraserDataSet(workingDir, dataset)
workingDir(fds) <- newWorkingdir
fds <- saveFraserDataSet(fds, dir=newWorkingdir, name=dataset, rewrite=TRUE)

#' ## Fit parameters
#' get fit params
fitMethod
psiType

#+ set bpparams
bpparam <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
register(bpparam)
# parallel(fds) <- bpparam

#'
#' ## Calculate zscores and p-values 
#'
#+ calculate zscores
fds <- calculateZscore(fds, type=psiType, logit=FALSE)

#+ calculate pvalues
message(date(), " Computing p-values ...")
fds <- calculatePvalues(fds=fds, type=psiType, implementation=fitMethod)

#+ adjust for multiple testing
message(date(), " Computing adjusted p-values ...")
fds <- calculatePadjValues(fds=fds, type=psiType, method=pmethod)

#+ save output fds
fds <- saveFraserDataSet(fds)
