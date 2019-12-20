#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  params:
#'   - progress: FALSE
#'  input:
#'   - dPsiSS: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/delta_psiSite.h5"`'
#'   - html:   '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_hyper_parameter_optimization.html"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_autoencoder_fit.html"`'
#'  threads: 20
#'  type: noindex
#'---
## TODO can not link to predictedMeans since BB method does not have a means object!
##
###'   - fdsout: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/predictedMeans_psiSite.h5"`'
##
## TODO:
##   Add a link to the fraser object to have a proper chain of events in wbuild
##

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", method="BB")
    parseWBHeader2("Scripts/FraseR/05_fit_autoencoder_FraseR.R", wildcards)
    slot(snakemake, "wildcards", check=FALSE) <- wildcards
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input
inputFile  <- snakemake@input$dPsiSS
dataset    <- snakemake@wildcards$dataset
method     <- snakemake@wildcards$method
name       <- basename(dirname(inputFile))
workingDir <- dirname(dirname(dirname(inputFile)))
bpWorkers   <- bpMaxWorkers(snakemake@threads)
bpProgress <- snakemake@params[[1]]$progress


#'
#' # Load PSI data
#+ echo=TRUE
dataset
method
workingDir

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=name)
BPPARAM <- MulticoreParam(bpWorkers, progressbar=bpProgress)
parallel(fds) <- BPPARAM

#'
#' # Fit autoencoder
#'

#'
#' run it for every type
#'
for(type in psiTypes){

    # set current type
    currentType(fds) <- type
    curDims <- dim(K(fds, type))
    q <- bestQ(fds, type)
    noise <- bestNoise(fds)
    probE <- max(0.001, min(1,30000/curDims[1]))

    # run autoencoder
    fds <- fit(fds, q=q, type=type, correction=method, verbose=TRUE,
            iterations=15, nSubset=15000, noiseAlpha=noise, BPPARAM=BPPARAM)

    # save autoencoder fit
    fds <- saveFraseRDataSet(fds)
}


#'
#' # Save results
#'
fds <- saveFraseRDataSet(fds)

