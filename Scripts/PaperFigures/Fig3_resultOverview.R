#'---
#' title: Paper figure (Result overview)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - fds:    '`sm config["DATADIR"] + "/datasets/savedObjects/Skin_Not_Sun_Exposed_Suprapubic__" + config["AE_IMPLEMENTATION"] + "/pajdBetaBinomial_psiSite.h5"`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/Figure3_resultOverview.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/Figure3_resultOverview.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    parseWBHeader2("Scripts/PaperFigures/Fig3_resultOverview.R", rerun=TRUE)
}

#+ input
dataset <- "Skin_Not_Sun_Exposed_Suprapubic"
workDir <- file.path(CONFIG$DATADIR, "datasets")
FraseR_implementation <- CONFIG$AE_IMPLEMENTATION
FDR_LIMIT <- 0.1
DELTA_PSI_LIMIT <- 0.1


#+ output
outPdf <- snakemake@output$outPdf
outPng <- snakemake@output$outPng

#+ load fds
fds <- loadFraserDataSet(workDir, paste0(dataset, "__", FraseR_implementation))
fds
outPng

#'
#' Sun not exposed non outlier example
#'
gr1 <- makeGRangesFromDataFrame(keep.extra.columns=TRUE, data.table(
    seqnames = "1",
    start    = 206567031,
    end      = 206574779,
    strand   = "+",
    type     = "psi5"))
gr1

xticks <- c(30, 100, 300)
yticks <- c(50, 100, 300, 500)
g1 <- plotExpression(fds, result=gr1, label=NULL) + 
    scale_y_log10(breaks=yticks, labels=yticks) + 
    scale_x_log10(breaks=xticks, labels=xticks) + 
    ylab("Split reads (K)\nfrom exon 7 to exon 8") + 
    xlab("All split reads (N)\nat exon 7 donor")
g1

g2 <- plotQQ(fds, result=gr1, label=NULL)
g2

#'
#' Sun not exposed outlier example
#' 
gr2 <- makeGRangesFromDataFrame(keep.extra.columns=TRUE, data.table(
    seqnames = "7",
    start    = 100485481,
    end      = 100485674,
    strand   = "+",
    type     = "psi5"))
gr2


xticks <- c(100, 300, 500)
yticks <- c(30, 50, 100)
g3 <- plotExpression(fds, result=gr2, deltaPsiCutoff=0.1, label=NULL) + 
    scale_y_log10(breaks=yticks, labels=yticks) + 
    scale_x_log10(breaks=xticks, labels=xticks) + 
    ylab("Split reads (K)\nfrom exon 17 to exon 18") + 
    xlab("All split reads (N)\nat exon 17 donor")
g3

g4 <- plotQQ(fds, result=gr2, deltaPsiCutoff=0.1, label=NULL)
g4

#'
#' overview number of outliers per sample
#' 
g5tmp <- plotQQ(fds, global=TRUE, aggregate=FALSE) +
    theme(plot.title=element_text(face="bold"))

g5legend <- get_legend(g5tmp + 
        theme(legend.position = "left") + 
        guides(color=guide_legend(element_blank())))
g5yRanges <- layer_scales(g5tmp)$y$range$range

g5 <- g5tmp + 
    theme_bw() + 
    theme(legend.position="none", plot.title=element_text(face="bold")) +
    annotation_custom(g5legend, xmin=0, xmax=2, 
            ymin=g5yRanges[1] + diff(g5yRanges)*0.6, ymax=g5yRanges[2])
g5

#'
#' Assemble figures
#'
#+ assemble figure for psi5
g <- ggarrange(ncol=2, nrow=1, widths = c(2,1), align="hv",
    ggarrange(ncol=2, nrow=2, labels=letters[1:4], align="hv",
        g1 + theme(plot.title=element_blank()),
        g2 + theme(plot.title=element_blank()),
        g3 + theme(plot.title=element_blank()),
        g4 + theme(plot.title=element_blank())),
    ggarrange(labels=letters[5],
        g5 + theme(plot.title=element_blank())))
g


#+ save heatmap figure
factor <- 0.6
ggsave(outPng, g, width = 17*factor, height = 10*factor)
ggsave(outPdf, g, width = 17*factor, height = 10*factor)

