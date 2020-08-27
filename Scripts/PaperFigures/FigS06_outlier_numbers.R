#'---
#' title: Number of events per method (Figure S6)
#' author: Christian Mertes
#' wb:
#'  threads: 5
#'  input:
#'   - fds:    '`sm expand(config["DATADIR"] + "/datasets/savedObjects/{dataset}__" + config["AE_IMPLEMENTATION"] + "/pajdBetaBinomial_psiSite.h5", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS6_outlier_numbers.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS6_outlier_numbers.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild_small.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS06_outlier_numbers.R",
            options=options, rerun=TRUE)
    threads <- 9
}

#+ input
fdsFiles      <- snakemake@input$fds
datasets      <- gsub("__.*", "", basename(dirname(fdsFiles)))
FRASER_imp    <- unique(gsub(".*__", "", basename(dirname(fdsFiles))))

names(fdsFiles) <- datasets

threads <- snakemake@threads
pvalCutoffs <- 10^-(3:6)
names(pvalCutoffs) <- paste0("italic(P)~value~10^-", 3:6)

length(fdsFiles)
fdsFiles[1:5]

#+ output
outPdf <- snakemake@output$outPdf
outPng <- snakemake@output$outPng


#'
#' ## Load P value FRASER
# f <- fdsFiles[1]
dt2p <- rbindlist(mclapply(fdsFiles, mc.cores=threads, mc.preschedule=FALSE,
    function(f){
        fds    <- loadFraserDataSet(file=f)
        tissue <- gsub("__.*", "", basename(dirname(f)))
        rbindlist(lapply(psiTypes, function(type){
            pv <- as.matrix(pVals(fds, type=type, byGroup=TRUE))
            rbindlist(lapply(pvalCutoffs, function(cut){
                cshit <- colSums(pv < cut)
                ans <- quantile(cshit, probs=c(0.5, 0.9))
                ans <- as.data.table(ans, keep.rownames=TRUE)
                ans[, .(quantile=rn, count=ans, cutoff=cut, type=type,
                        tissue=tissue, totalEvents=prod(dim(pv)), 
                        numHits=sum(cshit))]
            }))
        }))
    }))

# 
# better naming
# 
dt2p[,type:=factor(type, levels=c("psi5", "psi3", "psiSite"))]
levels(dt2p$type) <- c("psi[5]", "psi[3]", "theta")
dt2p[,cutoff:=factor(cutoff)]

names(pvalCutoffs) <- paste0("italic(P) < 10^-", 3:6)
levels(dt2p$cutoff) <- names(sort(pvalCutoffs))
dt2p[,cutoff:=factor(cutoff, levels=rev(levels(cutoff)))]

g1 <- ggplot(dt2p, aes(y=count, x=quantile, fill=type)) + 
    geom_violin() +
    facet_grid(cutoff ~ type, labeller=label_parsed) + 
    scale_fill_brewer(palette="Dark2", labels=function(l) parse(text=l)) + 
    scale_y_log10() + 
    ylab("Number of events") + 
    xlab("Quantile") + 
    labs(fill="Splice metric")
g1

g2 <- ggplot(dt2p, aes(y=numHits/totalEvents, x=type, fill=type)) +
    geom_violin() +
    facet_grid(~cutoff, labeller=label_parsed) + 
    scale_fill_brewer(palette="Dark2", labels=function(l) parse(text=l)) + 
    scale_x_discrete(labels=function(l) parse(text=l)) + 
    scale_y_log10() + 
    labs(fill="Splice metric") + 
    xlab("Splice metric") + 
    ylab("Frequency")
g2


#'
#' Assemble figures
#'
#+ assemble figure for psi5
g <- ggarrange(common.legend=TRUE, legend="right", labels=letters[1:2], ncol=1,
        heights=c(2,1),
        g1 + theme_bw(),
        g2 + theme_bw())
g

#+ save heatmap figure
factor <- 0.4
outPng
ggsave(outPng, g, width = 17*factor, height = 14*factor)
ggsave(outPdf, g, width = 17*factor, height = 14*factor)

