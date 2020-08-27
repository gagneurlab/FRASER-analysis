#'---
#' title: Figure S14 - Benchmark by swapping
#' author: Christian Mertes
#' wb:
#'   threads: 9
#'   input:
#'     - data: '`sm config["DATADIR"] + "/Benchmark/resBySwap/Skin_Cortex__0.0.RDS"`'
#'   output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS14_BenchmarkBySwapping.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS14_BenchmarkBySwapping.pdf"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")
source("./src/r/precisionRecallHelpers/precRec-helper.R")

# Interactive mode (debuging)
if(FALSE){
    load_wbuild()
    source(".wBuild/wBuildParser2.R")
    options   <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS14_benchmarkBySwapping.R",
            wildcards=wildcards, options=options, rerun=TRUE)
}

#' input
tissue        <- snakemake@wildcards$dataset
deltaPsi      <- as.numeric(snakemake@wildcards$deltaPsi)
data_file     <- snakemake@input$data

#' output
outPdf <- snakemake@output$outPdf
outPng <- snakemake@output$outPng


methods2plot  <- c(
    FRASER="PCA_p",
    `naïve BB`="BB_p",
    `Kremer et al.`="Leafcutter_p",
    LeafcutterMD="LeafcutterMD_p",
    SPOT="SPOT_p")


tissue
deltaPsi

data_file
outPdf

###########################################
#'
#' # Read in data
#'
###########################################
data <- readRDS(data_file)
names(data$data) <- c("dtmerged", "dt2p", "dt2p2")
attach(data$data)

###########################################
#'
#' plot results Precision/Recall
#' 
###########################################

gg1 <- plotRibbonBenchmark(dt2p, title="test", 
        linetype=rep(1, length(methods2plot)), 
        wrap_function=function() {})

gg1 <- gg1 + 
    guides(linetype=FALSE) + 
    scale_color_brewer(palette="Dark2") + 
    scale_fill_brewer(palette="Dark2") +
    theme_bw() + 
    grids() + 
    ylab("Precision") + 
    xlab("Recall") + 
    theme(plot.title=element_blank())

gg1

#'
#' plot numbers (recall)
#' 
dtmerged <- dtmerged[Method %in% methods2plot]
dtmerged[,Method:=factor(Method, levels=methods2plot)]
levels(dtmerged$Method) <- names(methods2plot)

dt2p2 <- rbindlist(lapply(1:10, function(x){
    dtmerged[!is.na(Method) & geneID != "",.(
            cutoff=10^(-x), 
            calls=sum(Score < 10^(-x)),
            total=.N,
            totalOutliers=sum(label),
            recall=sum(label == TRUE & Score < 10^(-x))),
                by="Method"][order(recall)]
}))

gg2 <- ggplot(dt2p2, aes(x=-log10(cutoff), y=recall/nrow(anno), col=Method)) +
    geom_line() +
    theme_bw() + 
    scale_color_brewer(palette="Dark2") + 
    grids() +
    ylim(0,1) +
    xlab(expression(paste(-log[10], "(", italic(P)~"value)"))) + 
    ylab("Recall")
gg2

gg3 <- ggplot(dt2p2, aes(x=-log10(cutoff), y=calls, col=Method)) +
    geom_line() +
    theme_bw() + 
    scale_color_brewer(palette="Dark2") + 
    grids() +
    xlab(expression(paste(-log[10], "(", italic(P)~"value)"))) + 
    ylab("Total number of calls") + 
    scale_y_log10()
gg3


#'
#' Assemble figures
#'
#+ assemble figure for psi5
g <- ggarrange(nrow=2, common.legend=TRUE, legend="bottom", labels=letters[1],
    heights=c(1.5,1),
    gg1, 
    ggarrange(ncol=2, labels=letters[2:3], legend="none",
        gg2,
        gg3))
g <- gg1 + theme(legend.position = "bottom")
g



#+ save figure
factor <- 0.6
outPng
ggsave(outPng, g, width = 11*factor, height = 11*factor)
ggsave(outPdf, g, width = 11*factor, height = 10*factor)

