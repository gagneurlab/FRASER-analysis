#'---
#' title: Count RNA data with FraseR
#' author: Christian Mertes
#' wb:
#'  params:
#'   - internalThreads: 3
#'   - progress: FALSE
#'  input:
#'   - colData: '`sm config["DATADIR"] + "/annotations/{dataset}.tsv"`'
#'   - html:    '`sm config["htmlOutputPath"] + "/annotations/{dataset}.html"`'
#'  output:
#'   - fdsobj:  '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'   - countsJ: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsS: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - wBhtml:  '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_counting.html"`'
#'  threads: 15
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")

# for interactive mode
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer")
    parseWBHeader2("./Scripts/FraseR/01_countRNA_FraseR.R",
            wildcards=wildcards, rerun=TRUE)
}

#+ input,
dataset     <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir  <- dirname(dirname(dirname(snakemake@output$countsJ)))
bpWorkers   <- bpMaxWorkers(snakemake@threads)
bpProgress  <- as.logical(snakemake@params[[1]]$progress)
iThreads    <- max(1, min(as.integer(c(bpworkers() / 5, snakemake@params[[1]]$internalThreads))))

#'
#' # Dataset
#+ echo=TRUE
dataset

#+ echo=FALSE
colData <- fread(colDataFile)
DT::datatable(colData, options=list(scrollX=TRUE))

genome <- NULL
if("genome" %in% colnames(colData)){
   genome <- colData$genome
   names(genome) <- colData$sampleID
}

#'
#' Counting the dataset
#'
BPPARAM <- MulticoreParam(bpWorkers, nrow(colData), progressbar=bpProgress)
register(BPPARAM)

fds <- FraseRDataSet(colData,
        workingDir = workingDir,
        name       = paste0("raw-", dataset))
fds <- countRNAData(fds, NcpuPerSample=iThreads, recount=TRUE, 
        minAnchor=5, genome=genome)
fds <- saveFraseRDataSet(fds)

#'
#' FraseR object after counting
#'
fds

