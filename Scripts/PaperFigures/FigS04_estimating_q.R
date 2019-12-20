#'---
#' title: Paper figure S04 Finding Q
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fds:    '`sm expand(config["DATADIR"] + "/datasets/savedObjects/{dataset}__" + config["AE_IMPLEMENTATION"] + "/pajdBetaBinomial_psiSite.h5", dataset=config["heatmap_tissues"])`'
#'   - stats:  '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_stats.RDS", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS4_finding_q.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS4_finding_q.pdf"`'
#'  threads: 5
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list()
    opt <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS04_estimating_q.R",
            wildcards=wildcards, options=opt, rerun=TRUE)
}


#+ input
ptype       <- snakemake@wildcards$psiType
fdsFiles    <- snakemake@input$fds
statsFiles  <- snakemake@input$stats
datasets    <- basename(dirname(fdsFiles))
tissues     <- gsub("__.*", "", datasets)
workingDir  <- dirname(dirname(dirname(fdsFiles[[1]])))
BPPARAM     <- MulticoreParam(5)

#+ output
outPng <- snakemake@output$outPng
outPdf <- snakemake@output$outPdf

#'
#' # Load datasets
#+ echo=TRUE
datasets
tissues

#+ echo=FALSE
fds_ls <- bplapply(datasets, loadFraseRDataSet, dir=workingDir, BPPARAM=BPPARAM)
names(fds_ls) <- basename(dirname(fdsFiles))

stats_ls <- lapply(statsFiles, readRDS)
names(stats_ls) <- basename(dirname(statsFiles))


#'
#' Extract q estimation data
#' 
dt2p <- rbindlist(lapply(names(fds_ls), function(name){
        rbindlist(lapply(psiTypes, function(type){
                dttmp <- hyperParams(fds_ls[[name]], type=type, all=TRUE)
                dttmp[,tissue:=name]
                dttmp[,type:=type]
                dttmp 
        }))}))

# add best q
dt2p[,maxQ:=q[which.max(aroc)], by="tissue,type"]

# better naming
dt2p[,tissue:=dName4plot(tissue, TRUE)]
dt2p[grepl("Skin not Sun Exposed", tissue), tissue:="Suprapubic Skin"]
dt2p[,tissue:=gsub(" ", "~", tissue)]
dt2p[,type:=factor(type, levels=c("psi5", "psi3", "psiSite"))]
levels(dt2p$type) <- c("psi[5]", "psi[3]", "theta")
g1 <- ggplot(dt2p, aes(q, aroc)) + 
    geom_point() + 
    geom_vline(aes(xintercept=maxQ, lty="Optimal q"), ) + 
    scale_linetype_manual(name="", label="Optimal q", values=2) + 
    geom_smooth() + 
    facet_grid(type ~ tissue, labeller=label_parsed) +
    xlab(expression(Latent~space~dimension~italic(q))) + 
    ylab("Area under the curve\nof precision-recall") + 
    theme_bw() + 
    theme(legend.position="bottom")
g1

#'
#' # Get Q across tissues
#'
best_qs <- rbindlist(lapply(psiTypes, function(type){
    q <- sapply(stats_ls, "[[", "ImportantNumbers")[paste0("Q_", type),]
    Nsamples <- sapply(stats_ls, "[[", "ImportantNumbers")["Nsamples",]
    data.table(type=type, q=q, Nsamples=Nsamples, tissues=names(values)) }))
best_qs[,type:=factor(type, levels=c("psi5", "psi3", "psiSite"))]
levels(best_qs$type) <- c("psi[5]", "psi[3]", "theta")
g2 <- ggplot(best_qs, aes(q)) + 
    geom_histogram(stats="identity") + 
    facet_wrap(~type, labeller=label_parsed) +
    theme_bw() + 
    xlab(expression(Latent~space~dimension~italic(q))) + 
    ylab("Frequency")
g2

#' 
#' Best Q correlation with sample size
#' 
g3 <- ggplot(best_qs, aes(Nsamples, q)) +
    geom_point() + 
    geom_smooth(method = "lm") + 
    facet_wrap(~type, labeller=label_parsed) + 
    theme_bw() + 
    xlab("Number of samples") + 
    ylab(expression(Latent~space~dimension~italic(q))) + 
    ylim(c(0, NA)) + 
    xlim(c(0, NA))
g3



#'
#' Assemble figure
#'
#+ assemble full figure
g <- ggarrange(ncol=1, heights=c(2,1), labels=c("A", "B"),
    g1, 
    g3)
g

#+ save heatmap figure
factor <- 0.5
outPng
ggsave(outPng, g, width = 17*factor, height = 14*factor)
ggsave(outPdf, g, width = 17*factor, height = 14*factor)

