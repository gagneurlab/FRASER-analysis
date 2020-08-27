#'---
#' title: Precision Recall per dataset
#' author: Ines Scheller
#' wb:
#'  py:
#'  - |
#'   a = config["N_bins"]
#'   b = config["dPsi_bins"]
#'   indices = range(1, a*b+1)
#'  input:
#'   - fitHtml: 'Output/html/zscoreCheck/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_pvalues.html'
#'  output:
#'   - plotData: '`sm expand(config["DATADIR"] + "/processedData/zscoreCheck/precRec/{{dataset}}/inject_{{delta}}/{{psiType}}/{{method}}_plotData_{{outlierType}}_{index}.tsv.gz", index=indices)`'
#'   - wBhtml: 'Output/html/zscoreCheck/OutlierRecall/{dataset}/{delta}/{psiType}/{method}_precRec_{outlierType}.html'
#'  type: noindex
#'---


#+ source main config
source("./src/r/config.R")
sourceFolder("src/r/precisionRecallHelpers")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", delta="uniformDistr", psiType="psi5", 
            method="PCA-BB-Decoder", outlierType="byJunctionGroup")
    options <- c("--config", "method='[\"PCA\", \"PCA-BB-Decoder\"]'")
    parseWBHeader2("Scripts/OutlierRecall/04_precision_recall.R", 
            wildcards=wildcards, rerun=TRUE, options=options, debug = TRUE)
}


#+ input
dataset     <- snakemake@wildcards$dataset
psiType     <- snakemake@wildcards$psiType
inj_delta   <- snakemake@wildcards$delta
method      <- snakemake@wildcards$method
outlierType <- snakemake@wildcards$outlierType
workingDir  <- file.path(snakemake@config$DATADIR, "datasets_zscoreCheck", 
        paste0("inject_", inj_delta), psiType, method)

#+ output
outTsvs <- snakemake@output$plotData

#+ Load fds
#' ## Load dataset
dataset
workingDir
fds       <- loadFraserDataSet(workingDir, dataset)
q         <- (snakemake@config$Qs)[1]
Nsamples  <- (snakemake@config$N_samples)[1]
pmethod   <- (snakemake@config$FDR_METHOD)[1]
nMeanBins <- as.integer(snakemake@config$N_bins)
dPsiBins  <- as.integer(snakemake@config$dPsi_bins)

# parallel(fds) <- SerialParam()
BPPARAM <- MulticoreParam(10, 10, progress=TRUE)
register(BPPARAM)

#+ Create evaluation table
# debug(createPrecRecTable)
plotDataList <- createPrecRecTable(fds=fds, psiType=psiType, 
        outliers=outlierType, nMeanBins=nMeanBins, dPsiBins=dPsiBins, 
        q=q, correction=method, anTrueCts='originalCounts', 
        Nsamples=Nsamples, inj_value=paste0('deltaPSI=',inj_delta), 
        pmethod=pmethod, dataset=dataset, BPPARAM=BPPARAM)

#' ## Show subset of results 
#' z score is based on deltaPSI
exampleData <- rbindlist(lapply(plotDataList, n=50, function(data, n){
    rbind(head(data, n), tail(data, n)) }))
DT::datatable(exampleData)

#+ Save plot data in tsv
message(date(), " Number of data.tables to write to file: ", length(plotDataList))
message(date(), " Number of output file names: ", length(outTsvs))
dev_null <- lapply(seq_along(plotDataList), FUN=function(i){
    message(date(), " Writing file: ", outTsvs[[i]])
    fwrite(x=plotDataList[[i]], file=outTsvs[[i]], quote=FALSE, 
            sep="\t", buffMB=1024)
})

