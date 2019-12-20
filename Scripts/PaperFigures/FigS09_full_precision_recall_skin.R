#'---
#' title: Plot full Precision recall on skin
#' author: Christian Mertes
#' wb:
#'  input:
#'   - plotData: '`sm expand(config["DATADIR"] + "/processedData/precRec/{dataset}/inject_{delta}/{psiType}/{outlierType}_outlier_recall_overall_ggplot.RDS", dataset=config["reference_dataset_injectedDpsiPlot"], psiType=config["psiTypes"], delta=config["inj_deltas"], outlierType=config["outlier_types"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS9_full_pr_plot.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS9_full_pr_plot.pdf"`'
#'  type: noindex
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    parseWBHeader2("Scripts/PaperFigures/FigS09_full_precision_recall_skin.R",
            rerun=TRUE)
    threads <- 16
}

#+ input
plotDataFiles <- snakemake@input$plotData

#+ output
outPng <- snakemake@output$outPng
outPdf <- snakemake@output$outPdf

plotDataFiles
outPng

#'
#' read data 
#' 
# x <- plotDataFiles[1]
data <- lapply(plotDataFiles, function(x){
    ans <- readRDS(x)
    psiType <- basename(dirname(x))
    ans$dt2p[,psiType:=basename(dirname(x))]
    ans$cutoffdt[,psiType:=basename(dirname(x))]
    ans
})
dt2p     <- rbindlist(lapply(data, "[[", "dt2p"))
cutoffdt <- rbindlist(lapply(data, "[[", "cutoffdt"))

renamePsiType <- function(dt){
    dt[psiType == "psi5",    psiType:="psi[5]"]
    dt[psiType == "psi3",    psiType:="psi[3]"]
    dt[psiType == "psiSite", psiType:="theta"]
    dt[,psiType:=factor(psiType, levels=c("psi[5]", "psi[3]", "theta"))]
    dt
}
renameMethods <- function(dt){
    levels(dt$Method) <- gsub("BB-Decoder", "BB-regression", levels(dt$Method))
    lvOrder <- c("FRASER\n+ P value", "PCA + BB-regression\n+ P value", 
            "naïve BB\n+ P value", "PCA + BB-regression\n+ z score", 
            "PCA\n+ z score")
    dt[,Method:=factor(Method, levels=lvOrder)]
    dt
}

dt2p <- renamePsiType(dt2p)
cutoffdt <- renamePsiType(cutoffdt)
dt2p <- renameMethods(dt2p)
cutoffdt <- renameMethods(cutoffdt)


#' 
#' plot precision recall plot
#' 
ggAll <- ggplot(dt2p, aes(recall, precision, color=Method)) + 
    geom_vline(xintercept=c(0.8, 0.85, 0.9, 0.95), col="gray90") + 
    geom_line() + 
    geom_point(data=cutoffdt, aes(x=recall, y=precision, shape=Cutoff), size=3) + 
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    theme_bw() +
    scale_color_brewer(palette='Dark2') +
    scale_fill_brewer(palette='Dark2') + 
    theme(legend.position = "right") + 
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Method), alpha=0.2, color=NA) + 
    facet_wrap(~psiType, ncol=1, labeller = label_parsed) + 
    ylab("Precision") + 
    xlab("Recall")
ggAll

#'
#' Some statistics for the paper
#' 
statsdt <- rbindlist(lapply(unique(dt2p$Method), function(x){
    rbindlist(lapply(unique(dt2p$psiType), function(xx){
        rbindlist(lapply(c(0.8, 0.85, 0.9), function(y){
            dt2p[Method == x & psiType == xx & recall > y][,.(Method, rank, TP,
                    precision=round(precision, 2), recall=round(recall, 2), cut=y, psiType)][1]
        }))}))}))
statsdt[,fc:=round(precision/precision[grepl("FRASER", Method)], 2),by="cut,psiType"]
statsdt[order(cut, Method)][grepl("(PCA|FRASER)\\n", Method, perl=TRUE)]
cutoffdt[,.(Method, rank, TP, precision=round(precision, 2), recall=round(recall, 2), Cutoff, psiType)][grepl("(PCA|FRASER)\\n", Method, perl=TRUE)]

#+ save plot as PDF and PNG
factor <- 0.6
ggsave(outPng, ggAll, 'png', width=10*factor, height=10*factor)
ggsave(outPdf, ggAll, 'pdf', width=10*factor, height=10*factor)
