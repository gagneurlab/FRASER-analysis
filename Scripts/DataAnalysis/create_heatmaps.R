#'---
#' title: Heatmaps
#' author: Ines Scheller
#' wb:
#'  params:
#'   - workers: 1
#'  input:
#'   - normFDS:  '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__FraseR-1DecoderBatches/predictedMeans_psiSite.h5"`'
#'  output:
#'   - outPngs: '`sm expand(config["FIGDIR"] + "/heatmaps/{{dataset}}/heatmap{topN}_{psiType}.png", topN=["5k", "10k"], psiType=config["psiTypes"])`'
#'   - outPlotData: '`sm expand(config["FIGDATADIR"] + "/heatmaps/{{dataset}}/{heatmapType}_{norm}_{psiType}.RDS", heatmapType=["sampleCorrelation", "junctionSample5k", "junctionSample10k"], psiType=config["psiTypes"], norm=["raw", "normalized"])`'
#'   - wBhtml: 'Output/html/DataAnalysis/heatmaps/{dataset}_heatmaps.html'
#'  type: noindex
#'---


if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/DataAnalysis/create_heatmaps.R", dataset="Skin_Not_Sun_Exposed_Suprapubic")
}

#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=31, fig.height=28)

#+ input
correctionMethod <- "FraseR-1DecoderBatches"
dataset     <- paste0(snakemake@wildcards$dataset, "__", correctionMethod)
workingDir  <- file.path(CONFIG$DATADIR, "datasets")
bpWorkers   <- min(bpworkers(), as.integer(snakemake@params[[1]]$workers))

#+ output
outPlotData <- snakemake@output$outPlotData
outPngs     <- snakemake@output$outPngs
outPng5k    <- outPngs[which(grepl(pattern="5k", x=outPngs))]
outPng10k   <- outPngs[which(grepl(pattern="10k", x=outPngs))]

#'
#' # Load count data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)

#'
#' Sample correlations of counts and junctions x samples heatmap
#'
require(gridExtra)
require(ggplot2)

# colLevels <- apply(colData(fds), 2, FUN=function(x){
#     length(levels(as.factor(x)))
# })
# annotation_cols_2plot <- names(colLevels[colLevels < 30 & colLevels > 1])
# annotation_cols_2plot
annotation_cols_2plot <- c("DTHHRDY", "AGE", "GENDER", "SMRIN", "SMNABTCHT", "SMCENTER")

#+ create heatmaps for all psi types
plotTypes <- c("sampleCorrelation", "junctionSample")
norm=c("raw", "normalized")
plots <- lapply(psiTypes, function(type){
    plot_list=list()
    for(pt in plotTypes){
        topN <- 100000
        topJ <- 5000
        for(n in norm){
            x=plotCountCorHeatmap(fds=fds, type=type, logit=TRUE, topN=topN,
                                  topJ=topJ, plotType=pt, normalized=(n=="normalized"),
                                  annotation_col=annotation_cols_2plot,
                                  annotation_row=NA, nClust=5 )
            plot_list[[paste0(pt, "_", n)]] = x[[4]]     ##to save each plot into a list. note the [[4]]
            outType <- ifelse(pt == "junctionSample", "junctionSample5k", pt)
            saveRDS(x, file=outPlotData[which(grepl(x=outPlotData, pattern=paste0(outType, ".*", n, ".*", type)))])
        }
    }
    g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
    g
    ggsave(outPng5k[which(grepl(pattern=type, x=outPng5k))], g, width = 31, height = 28)

    # also plot top 10k junctions
    topJ <- 10000
    plot_list=list()
    for (n in norm){
        x=plotCountCorHeatmap(fds=fds, type=type, logit=TRUE, topN=topN, topJ=topJ,
                              plotType="junctionSample", normalized=(n=="normalized"),
                              annotation_col=annotation_cols_2plot,
                              annotation_row=NA, nClust=5 )
        plot_list[[n]] = x[[4]]     ##to save each plot into a list. note the [[4]]
        saveRDS(x, file=outPlotData[which(grepl(x=outPlotData, pattern=paste0("junctionSample10k", ".*", n, ".*", type)))])
    }
    g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
    g
    ggsave(outPng10k[which(grepl(pattern=type, x=outPng10k))], g, width = 31, height = 14)
})

