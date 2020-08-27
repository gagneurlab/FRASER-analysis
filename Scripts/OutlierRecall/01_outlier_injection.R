#'---
#' title: Inject artificial outliers into a dataset
#' author: Ines Scheller
#' wb:
#'  input:
#'   - inFile: 'Output/html/FraseR/{dataset}_hyper_parameter_optimization.html'
#'  output:
#'   - fdsOut: '`sm config["DATADIR"] + "/datasets/inject_{delta}/{psiType}/savedObjects/{dataset}/fds-object.RDS"`'
#'   - trueOut: '`sm config["DATADIR"] + "/datasets/inject_{delta}/{psiType}/savedObjects/{dataset}/trueOutliers_{psiType}.h5"`'
#'   - wBhtml:  'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/outlierInjection.html'
#'  type: noindex
#'  threads: 10
#'---

#+ source main config
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", delta="uniformDistr", psiType="psiSite")
    parseWBHeader2("Scripts/OutlierRecall/outlier_injection.R", 
            rerun=TRUE, wildcards=wildcards)
}

#+ input
fdsFile    <- snakemake@output$fdsOut
workingDir <- dirname(dirname(dirname(dirname(dirname(fdsFile)))))
dataset    <- snakemake@wildcards$dataset
psiType    <- snakemake@wildcards$psiType
inj_delta  <- snakemake@wildcards$delta

bpWorkers  <- min(bpworkers(), snakemake@threads)
bpThreads  <- bpWorkers
bpProgress <- FALSE

#' ## Load dataset
#+ load fds
dataset
workingDir
fds     <- loadFraseRDataSet(workingDir, dataset)
dim(fds)

#' ## Copy fds, remove unneccessary assays and save to new workingDir of copy
#+ copy fds
copy_fds <- fds
workingDir(copy_fds) <- dirname(dirname(dirname(fdsFile)))
workingDir(copy_fds)

keep <- c("rawCountsJ", "rawOtherCounts_psi5", "rawOtherCounts_psi3",
        "rawCountsSS", "rawOtherCounts_psiSite", "psi3", "psi5", "psiSite")
if(grepl("Simulation", dataset)){
    keep <- c(keep, "truePSI_psi3", "truePSI_psi5", "truePSI_psiSite")
}
assays(copy_fds) <- assays(copy_fds)[keep]
keep <- c("hyperParams_psi5", "hyperParams_psi3", "hyperParams_psiSite",
        "currentType", "noiseAlpha", "dim")
metadata(copy_fds) <- metadata(copy_fds)[keep]

copy_fds <- saveFraseRDataSet(copy_fds)

#+ get injection params
#' ## Psi type
psiType

#' ## Injection frequency
outlier_freq <- CONFIG$OUTLIER_INJ_FREQ
if(dataset == "example"){
    outlier_freq <- 0.1
}
outlier_freq

#' ## Injection method
injMethod <- ifelse(grepl("Simulation", dataset), "simulatedPSI", "samplePSI")
injMethod

#' ## Injection delta psi
inj_delta

#' ## Min delta psi for injection
minDpsi <- CONFIG$MIN_DPSI
minDpsi

#+ set bpparams
bpparam <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
register(bpparam)
parallel(copy_fds) <- bpparam

#+ inject outliers
#' ## Outlier injection
copy_fds <- injectOutliers(copy_fds, type=psiType, freq=outlier_freq, 
        minCoverage=10, minDpsi=minDpsi, deltaDistr=inj_delta, 
        method=injMethod, verbose=TRUE)
copy_fds <- saveFraseRDataSet(copy_fds)

#+ get nr of injected outliers
#' ## Number of injected primary outliers
nrPrimOut <- sum(abs(getAssayMatrix(copy_fds, "trueOutliers", psiType)) == 1)
nrPrimOut

#' ## Number of injected secondary outliers
nrSecOut <- sum(abs(getAssayMatrix(copy_fds, "trueOutliers", psiType)) == 2)
nrSecOut

