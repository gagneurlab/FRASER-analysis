#'---
#' title: Hyper parameter optimization for different choice of outlier injections
#' author: Ines Scheller
#' wb:
#'  input:
#'   - dPsiSS: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/delta_psiSite.h5"`'
#'   - html:   '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_filterExpression.html"`'
#'  output:
#'   - encDimTable: '`sm config["DATADIR"] + "/processedData/results/{dataset}/{method}_qEstimation/old_qAuroc_{injDistr}_{injStartpoint}_{dpsi}_subset{subset}_{psiType}.tsv"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_qEstimation/old_{injDistr}_{injStartpoint}_{dpsi}_subset{subset}_{psiType}.html"`'
#'  threads: 20
#'  type: noindex
#'---


if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", method="FraseR")
    parseWBHeader2("Scripts/FraseR/qEstimation_outlierDistribution.R", 
                   wildcards=wildcards, rerun=TRUE)
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
method     <- snakemake@wildcards$method
type      <- snakemake@wildcards$psiType
name       <- paste0(dataset, "__", method)
workingDir <- file.path(snakemake@config$DATADIR, "datasets")
bpWorkers   <- bpMaxWorkers(snakemake@threads)
bpThreads  <- 100
bpProgress <- FALSE

injectionDistr <- snakemake@wildcards$injDistr
injectionStartpoint <- snakemake@wildcards$injStartpoint
deltaPsi <- as.numeric(snakemake@wildcards$dpsi)
subset   <- as.logical(snakemake@wildcards$subset)
deltaDistr <- switch(injectionDistr, 
                    "uniform" = "uniformDistr",
                    "fixed" = deltaPsi)


#+ output
outfile <- snakemake@output$encDimTable


#'
#' # Load PSI data
#+ echo=TRUE
dataset
name

#+ echo=FALSE
fds <- loadFraserDataSet(dir=workingDir, name=name)

BPPARAM <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
gc()

#'
#' # Run hyper parameterization
#'
encDimTable <- data.table()
# for(type in psiTypes){
    message(date(), ": ", type)
    
    # define q's to test
    # mp <- 6
    # a <- 2 
    # b <- min(ncol(fds), nrow(fds)) / mp   # N/mp
    # maxSteps <- 15
    # Nsteps <- min(maxSteps, b)
    # q_param <- unique(round(exp(seq(log(a),log(b),length.out = Nsteps))))
    # q_param <- sort(unique(c(q_param, 50, 70, 90)))
    q_param     <- c(seq(2, 20, by=2), seq(18, 40, by=5), 50, 70)
    q_param     <- sort(unique(pmin(q_param, ncol(fds))))
    noise_param <- c(0.5) # c(0.5, 1.5)
    
    # define size of subset for hyper param optimization
    setSubset <- ifelse(isTRUE(subset), 15000, nrow(K(fds, type)))
    
    fds <- optimHyperParams(fds, type=type, implementation=method, iterations=5,
                            q_param=q_param, noise_param=noise_param, 
                            minDpsi=deltaPsi, setSubset = setSubset,
                            BPPARAM=BPPARAM, internalThreads=1,
                            deltaDistr=deltaDistr, method=injectionStartpoint )
    dttmp <- hyperParams(fds, type=type, all=TRUE)
    dttmp[,dataset:=dataset]
    dttmp[,method:=method]
    dttmp[,type:=type]
    dttmp[,bestQ:=q[which.max(aroc)]]
    dttmp[,injectionDistr:=injectionDistr]
    dttmp[,injectionStartpoint:=injectionStartpoint]
    dttmp[,injectedDpsi:=deltaPsi]
    dttmp[,subsetted:=subset]
    encDimTable <- rbind(encDimTable, dttmp)
# }

#'
#' Write table with auroc values for each q
#'
fwrite(encDimTable, outfile)

