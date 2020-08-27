#'---
#' title: Paper figure S05 QQ plots
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fds:    '`sm config["DATADIR"] + "/datasets/savedObjects/Skin_Not_Sun_Exposed_Suprapubic__" + config["AE_IMPLEMENTATION"] + "/pajdBetaBinomial_psiSite.h5"`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS5_qq_plots.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS5_qq_plots.pdf"`'
#'  threads: 5
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list()
    opt <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS05_qq_plots.R",
            wildcards=wildcards, options=opt, rerun=TRUE)
}


#+ input
ptype       <- "psi5"
fdsFile     <- snakemake@input$fds
dataset     <- basename(dirname(fdsFile))
tissue      <- gsub("__.*", "", dataset)
workingDir  <- dirname(dirname(dirname(fdsFile)))
BPPARAM     <- MulticoreParam(5, progressbar=TRUE)

#+ output
outPng <- snakemake@output$outPng
outPdf <- snakemake@output$outPdf

#'
#' # Load datasets
#+ echo=TRUE
dataset
tissue
nqqplots <- 12
seed <- 42
    
#+ echo=FALSE
fds <- loadFraserDataSet(dir=workingDir, dataset)
set.seed(seed)

#'
#' get random 24 QQ plots
#' 
rmc <- rowMeans(K(fds, ptype))
rsdp <- rowSds(assay(fds, ptype))
table(rmc > 10 & rsdp > 0.02)

v2p <- sort(sample(which(rmc > 10 & rsdp > 0.02), nqqplots))

qqplot_ls <- bplapply(v2p, plotQQ, fds=fds, type=ptype, 
        BPPARAM=SerialParam(progressbar=TRUE))
qqplot_ls <- lapply(qqplot_ls, function(x){
    x + theme(legend.position="none", axis.title.x=element_blank(),
            axis.title.y=element_blank(), title=element_blank()) })
g1 <- ggarrange(plotlist=qqplot_ls, ncol=2, nrow=ceiling(nqqplots/2), 
        common.legend=FALSE, labels=letters[1:nqqplots],
        hjust=0.2, vjust=0.5)
g1


expplot_ls <- bplapply(v2p, plotExpression, fds=fds, type=ptype, 
        BPPARAM=SerialParam(progressbar=TRUE))
expplot_ls <- lapply(expplot_ls, function(x){
    x + theme(legend.position="none", axis.title.x=element_blank(),
            axis.title.y=element_blank(), title=element_blank()) })
g2 <- ggarrange(plotlist=expplot_ls, ncol=2, nrow=ceiling(nqqplots/2), 
        common.legend=FALSE, labels=letters[1:nqqplots + nqqplots], 
        hjust=0.2, vjust=0.5)
g2

#'
#' Assemble figure
#'
#+ assemble full figure
g <- ggarrange(nrow=3, heights=c(0.4, 23, 1),
    ggplot() + theme_nothing(),
    ggarrange(ncol=4, widths=c(1, 10, 1, 10),
        grid.text(expression(-log[10]~"(observed"~italic(P)~")"), 
                rot=90, gp=gpar(fontsize=14)),
        g1,
        grid.text("Junction count + 1 (K)", rot=90, gp=gpar(fontsize=14)),
        g2),
    ggarrange(ncol=4, widths=c(1, 10, 1, 10),
        ggplot() + theme_nothing(),
        grid.text(expression(-log[10]~"(expected"~italic(P)~")"), 
                gp=gpar(fontsize=14)),
        ggplot() + theme_nothing(),
        grid.text("Total junction coverage + 2 (N)", gp=gpar(fontsize=14))))
g
    
#+ save heatmap figure
factor <- 0.5
outPng
ggsave(outPng, g, width = 17*factor, height = 16*factor)
ggsave(outPdf, g, width = 17*factor, height = 16*factor)

