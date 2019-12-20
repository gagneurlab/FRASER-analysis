#'---
#' title: Figure 2b recall / bootstraped
#' author: Ines Scheller
#' wb:
#'   threads: 40
#'   py:
#'   - |
#'    a = config["N_bins"]
#'    b = config["dPsi_bins"]
#'    indices = range(1, a*b+1)
#'   input:
#'     - evalTable:     '`sm expand(config["DATADIR"] + "/processedData/precRec/{{dataset}}/inject_{{delta}}/{{psiType}}/{method}_plotData_{{outlierType}}_{index}.tsv.gz", method=config["methods"], index=indices)`'
#'   output:
#'    - figurePdf: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall.pdf"`'
#'    - figurePng: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall.png"`'
#'    - ggplot:    '`sm config["DATADIR"] + "/processedData/precRec/{dataset}/inject_{delta}/{psiType}/{outlierType}_outlier_recall_ggplot.RDS"`'
#'    - figurePdf_noConf: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_noConf.pdf"`'
#'    - figurePng_noConf: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_noConf.png"`'
#'    - ggplot_noConf:    '`sm config["DATADIR"] + "/processedData/precRec/{dataset}/inject_{delta}/{psiType}/{outlierType}_outlier_recall_noConf_ggplot.RDS"`'
#'    - figurePdf_sub: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_methodSubset.pdf"`'
#'    - figurePng_sub: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_methodSubset.png"`'
#'    - ggplot_sub:    '`sm config["DATADIR"] + "/processedData/precRec/{dataset}/inject_{delta}/{psiType}/{outlierType}_outlier_recall_methodSubset_ggplot.RDS"`'
#'    - figurePdf_sub_noConf: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_noConf_methodSubset.pdf"`'
#'    - figurePng_sub_noConf: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_noConf_methodSubset.png"`'
#'    - ggplot_sub_noConf:    '`sm config["DATADIR"] + "/processedData/precRec/{dataset}/inject_{delta}/{psiType}/{outlierType}_outlier_recall_noConf_methodSubset_ggplot.RDS"`'
#'    - outPlotDataPr:  '`sm config["FIGDATADIR"] + "/precRec/{dataset}/{delta}/{outlierType}_precision_recall_{psiType}.RDS"`'
#'    - outPlotCutoffs: '`sm config["FIGDATADIR"] + "/precRec/{dataset}/{delta}/{outlierType}_pr_cutoffs_{psiType}.RDS"`'
#'    - wBhtml:    'Output/html/OutlierRecall/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall.html'
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
    parseWBHeader2("Scripts/OutlierRecall/Figure2b_per_psiType.R", 
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
ggplotFile <- snakemake@output$ggplot
outPdf     <- snakemake@output$figurePdf
outPng     <- snakemake@output$figurePng
outData    <- snakemake@output$outPlotData
outCutoffs <- snakemake@output$outPlotCutoffs

ggplotFile_noConf <- snakemake@output$ggplot_noConf
outPdf_noConf     <- snakemake@output$figurePdf_noConf
outPng_noConf     <- snakemake@output$figurePng_noConf

ggplotFileSub <- snakemake@output$ggplot_sub
outPdfSub     <- snakemake@output$figurePdf_sub
outPngSub     <- snakemake@output$figurePng_sub

ggplotFileSub_noConf <- snakemake@output$ggplot_sub_noConf
outPdfSub_noConf     <- snakemake@output$figurePdf_sub_noConf
outPngSub_noConf     <- snakemake@output$figurePng_sub_noConf

#' ## Input and output file
recallFiles
ggplotFile

#+ Extract data in parallel
prDatals <- mclapply(recallFiles, mc.allow.recursive=TRUE, mc.preschedule=FALSE, mc.cores=min(max.mc, length(recallFiles)),
                     FUN=readBootstrapData, pvalue=TRUE, zscoreForAll=TRUE, FDR_LIMI=FDR_LIMIT, dPsiCutoff=dPsiCutoff)
zsDatals <- mclapply(recallFiles, mc.cores=min(max.mc, length(recallFiles)), readZscoreCutoff)
pcDatals <- mclapply(recallFiles, mc.cores=min(max.mc, length(recallFiles)), readPadjCutoff, FDR_LIMIT=FDR_LIMIT, dPsiCutoff=dPsiCutoff)
dpDatals <- mclapply(recallFiles, mc.cores=min(max.mc, length(recallFiles)), readDPSICutoff)

prData <- rbindlist(prDatals)
zsData <- rbindlist(zsDatals)
zsData <- rbind(zsData, rbindlist(pcDatals))
zsData <- rbind(zsData, rbindlist(dpDatals))

q <- unique(prData$q)

# correct order of facets in the plot
prData$inj_value <- factor(prData$inj_value, levels=sortIntervals(unique(prData$inj_value)))
zsData$inj_value <- factor(zsData$inj_value, levels=sortIntervals(unique(zsData$inj_value)))
prData$junctionMeanBin <- factor(prData$junctionMeanBin, levels=sortIntervals(unique(prData$junctionMeanBin)))
zsData$junctionMeanBin <- factor(zsData$junctionMeanBin, levels=sortIntervals(unique(zsData$junctionMeanBin)))
prData$type <- factor(prData$type, levels=c(paste('P-value rank & FDR <', FDR_LIMIT, '\n& deltaPSI >', dPsiCutoff), "zScore rank", "DeltaPSI rank"))
zsData$type <- factor(zsData$type, levels=c(paste('P-value rank & FDR <', FDR_LIMIT, '\n& deltaPSI >', dPsiCutoff), "zScore rank", "DeltaPSI rank"))

#' ## Figure 2b
#+ Create Figure 2b
ggAll <- plotRibbonBenchmark(prData, zscoreData=zsData, linetype=c(1,2,3), maxRows=maxRows, title=paste0(dataset, " outlier recall (", psiType, ", q = ", q ,")"),
                             wrap_function=function() facet_grid(inj_value ~ junctionMeanBin,
                                                                 labeller = label_bquote(rows = "outlier" ~ Delta * Psi == .(as.vector(inj_value)),
                                                                                         cols = "mean coverage" == .(as.vector(junctionMeanBin))) ) ) +
    scale_color_brewer(palette='Dark2') +
    scale_fill_brewer(palette='Dark2')
ggAll

#+ Create Figure 2b (no confidence intervals)
ggNoConfIntervals <- plotRibbonBenchmark(prData, zscoreData=zsData, linetype=c(1,2,3), maxRows=maxRows, title=paste0(dataset, " outlier recall (", psiType, ", q = ", q ,")"),
                                         confidenceIntervals=FALSE, wrap_function=function() facet_grid(inj_value ~ junctionMeanBin,
                                                                                                        labeller = label_bquote(rows = "outlier" ~ Delta * Psi == .(as.vector(inj_value)),
                                                                                                                                cols = "mean coverage" == .(as.vector(junctionMeanBin))) ) ) +
    scale_color_brewer(palette='Dark2') +
    scale_fill_brewer(palette='Dark2')
ggNoConfIntervals

#' ## Figure 2b (method subset)
#+ Create Figure 2b (method subset)
ggSub <- plotRibbonBenchmark(prData[correction %in% final_methods,], zscoreData=zsData[correction %in% final_methods,], linetype=c(1,2,3), maxRows=maxRows,
                             title=paste0(dataset, " outlier recall (", psiType, ", q = ", q ,")"),
                             wrap_function=function() facet_grid(inj_value ~ junctionMeanBin,
                                                                 labeller = label_bquote(rows = "outlier" ~ Delta * Psi == .(as.vector(inj_value)),
                                                                                         cols = "mean coverage" == .(as.vector(junctionMeanBin))) ) ) +
    scale_color_brewer(palette='Dark2') +
    scale_fill_brewer(palette='Dark2')
ggSub

#+ Create Figure 2b (method subset, no confidence intervals)
ggSubNoConfIntervals <- plotRibbonBenchmark(prData[correction %in% final_methods,], zscoreData=zsData[correction %in% final_methods,], linetype=c(1,2,3), maxRows=maxRows,
                                            title=paste0(dataset, " outlier recall (", psiType, ", q = ", q ,")"), confidenceIntervals=FALSE,
                                            wrap_function=function() facet_grid(inj_value ~ junctionMeanBin,
                                                                                labeller = label_bquote(rows = "outlier" ~ Delta * Psi == .(as.vector(inj_value)),
                                                                                                        cols = "mean coverage" == .(as.vector(junctionMeanBin))) ) ) +
    scale_color_brewer(palette='Dark2') +
    scale_fill_brewer(palette='Dark2')
ggSubNoConfIntervals


#+ Save figures as PNG and PDF
ggsave(outPdf, ggAll, 'pdf', width=12, height=8)
ggsave(outPng, ggAll, 'png', width=12, height=8, dpi=900)
saveRDS(ggAll, file=ggplotFile)

ggsave(outPdfSub, ggSub, 'pdf', width=12, height=8)
ggsave(outPngSub, ggSub, 'png', width=12, height=8, dpi=900)
saveRDS(ggSub, file=ggplotFileSub)

ggsave(outPdf_noConf, ggNoConfIntervals, 'pdf', width=12, height=8)
ggsave(outPng_noConf, ggNoConfIntervals, 'png', width=12, height=8, dpi=900)
saveRDS(ggNoConfIntervals, file=ggplotFile_noConf)

ggsave(outPdfSub_noConf, ggSubNoConfIntervals, 'pdf', width=12, height=8)
ggsave(outPngSub_noConf, ggSubNoConfIntervals, 'png', width=12, height=8, dpi=900)
saveRDS(ggSubNoConfIntervals, file=ggplotFileSub_noConf)

#+ Save plot data
saveRDS(prData, file=outData)
saveRDS(zsData, file=outCutoffs)
