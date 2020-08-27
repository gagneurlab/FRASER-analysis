#'---
#' title: Create pr input data for z scores in natural scale
#' author: Ines Scheller
#' wb:
#'   threads: 40
#'   py:
#'   - |
#'    a = config["N_bins"]
#'    b = config["dPsi_bins"]
#'    indices = range(1, a*b+1)
#'   input:
#'     - evalTable:     '`sm expand(config["DATADIR"] + "/processedData/zscoreCheck/precRec/{{dataset}}/inject_{{delta}}/{{psiType}}/{{method}}_plotData_{{outlierType}}_{index}.tsv.gz", index=indices)`'
#'   output:
#'    - outPlotDataPr:  '`sm config["FIGDATADIR"] + "/zscoreCheck/precRec/{dataset}/{delta}/{outlierType}_pr_only_{method}_{psiType}.RDS"`'
#'    - outPlotCutoffs: '`sm config["FIGDATADIR"] + "/zscoreCheck/precRec/{dataset}/{delta}/{outlierType}_pr_cutoffValues_only_{method}_{psiType}.RDS"`'
#'    - wBhtml:    'Output/html/zscoreCheck/OutlierRecall/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_only_{method}.html'
#'   type: noindex
#'---


#+ load config
source('./src/r/config.R')
sourceFolder("src/r/precisionRecallHelpers")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Kremer", delta="uniformDistr", psiType="psi5", 
                      outlierType="byJunctionGroup")
    options <- c("--config", "method='[\"PCA\", \"PCA-BB-Decoder\"]'")
    parseWBHeader2("Scripts/OutlierRecall_zscore/05_prepare_PR_for_plotting.R", 
                   wildcards=wildcards, rerun=TRUE, options=options, debug = TRUE)
}

register(MulticoreParam(20))
max.mc <- 3

#+ input
recallFiles <- snakemake@input$evalTable
dataset     <- snakemake@wildcards$dataset
inj_delta   <- snakemake@wildcards$delta
psiType     <- snakemake@wildcards$psiType
outlierType <- snakemake@wildcards$outlierType
FDR_LIMIT   <- CONFIG$FDR_LIMIT
maxRows     <- CONFIG$maxRows
dPsiCutoff  <- CONFIG$dPsiCutoff
final_methods <- CONFIG$final_methods

#+ output
outData    <- snakemake@output$outPlotDataPr
outCutoffs <- snakemake@output$outPlotCutoffs

#' ## Input and output file
recallFiles

#+ Extract data in parallel
prDatals <- mclapply(recallFiles, mc.allow.recursive=TRUE, mc.preschedule=FALSE, mc.cores=min(max.mc, length(recallFiles)),
                     FUN=readBootstrapData, FDR_LIMI=FDR_LIMIT, dPsiCutoff=dPsiCutoff, rankType=c("zscore"))
prData <- rbindlist(prDatals)

zsDatals <- mclapply(recallFiles, mc.cores=min(max.mc, length(recallFiles)), readZscoreCutoff)
zsData <- rbindlist(zsDatals)

# correct order of facets in the plot
prData$inj_value <- factor(prData$inj_value, levels=sortIntervals(unique(prData$inj_value)))
zsData$inj_value <- factor(zsData$inj_value, levels=sortIntervals(unique(zsData$inj_value)))
prData$junctionMeanBin <- factor(prData$junctionMeanBin, levels=sortIntervals(unique(prData$junctionMeanBin)))
zsData$junctionMeanBin <- factor(zsData$junctionMeanBin, levels=sortIntervals(unique(zsData$junctionMeanBin)))
prData$type <- factor(prData$type, levels=c(paste0('P-value rank & FDR < ', FDR_LIMIT, '\n& deltaPSI > ', dPsiCutoff), paste0("zScore rank\n& deltaPSI > ", dPsiCutoff), "DeltaPSI rank"))
zsData$type <- factor(zsData$type, levels=c(paste0('P-value rank & FDR < ', FDR_LIMIT, '\n& deltaPSI > ', dPsiCutoff), paste0("zScore rank\n& deltaPSI > ", dPsiCutoff), "DeltaPSI rank"))


#+ Save plot data
saveRDS(prData, file=outData)
saveRDS(zsData, file=outCutoffs)
