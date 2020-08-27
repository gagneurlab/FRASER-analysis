#'---
#' title: Paper figure S04 Finding Q (new for revision)
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fds:    '`sm expand(config["DATADIR"] + "/datasets/savedObjects/{dataset}__" + config["AE_IMPLEMENTATION"] + "/pajdBetaBinomial_psiSite.h5", dataset=config["heatmap_tissues"])`'
#'   - stats:  '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_stats.RDS", dataset=config["EnrichmentTissues"])`'
#'   - encDimTable:  '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + 
#'                               "_qEstimation/old_qAuroc_{injDistr}_{injStartpoint}_{dpsi}_subset{subset}_{psiType}.tsv", 
#'                               dataset=config["heatmap_tissues"], injDistr=["fixed"], 
#'                               injStartpoint=config["injection_startpoint"], 
#'                               dpsi=config["injection_dpsi"], subset=config["injection_subsetting"], psiType=config["psiTypes"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS4_finding_q_new.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS4_finding_q_new.pdf"`'
#'  threads: 5
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list()
    opt <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS04_estimating_q_new.R",
                   wildcards=wildcards, options=opt, rerun=TRUE)
}


#+ input
ptype       <- snakemake@wildcards$psiType
fdsFiles    <- snakemake@input$fds
statsFiles  <- snakemake@input$stats
datasets    <- basename(dirname(fdsFiles))
tissues     <- gsub("__.*", "", datasets)
workingDir  <- dirname(dirname(dirname(fdsFiles[[1]])))
encDimTableFiles  <- snakemake@input$encDimTable
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
fds_ls <- bplapply(datasets, loadFraserDataSet, dir=workingDir, BPPARAM=BPPARAM)
names(fds_ls) <- basename(dirname(fdsFiles))

stats_ls <- lapply(statsFiles, readRDS)
names(stats_ls) <- basename(dirname(statsFiles))

encDimTables <- rbindlist(lapply(encDimTableFiles, fread))
# encDimTables <- encDimTables[injectionDistr == "fixed",]

#'
#' Extract q estimation data
#' 
dt2p <- rbindlist(lapply(names(fds_ls), function(name){
    rbindlist(lapply(psiTypes, function(type){
        dttmp <- hyperParams(fds_ls[[name]], type=type, all=TRUE)
        dttmp[,dataset:=name]
        # dttmp[,tissue:=name]
        dttmp[,method:=CONFIG$AE_IMPLEMENTATION]
        dttmp[,type:=type]
        dttmp[,bestQ:=q[which.max(aroc)]]
        dttmp[,injectionDistr:="uniform"]
        dttmp[,injectionStartpoint:="samplePSI"]
        dttmp[,injectedDpsi:=0.2]
        dttmp[,subsetted:=TRUE]
        dttmp 
    }))}))

dt2p <- rbind(dt2p, encDimTables)

# add best q
# dt2p[,maxQ:=q[which.max(aroc)], by="tissue,type"]

# better naming
dt2p[,dataset:=dName4plot(dataset, TRUE)]
dt2p[grepl("Skin not Sun Exposed", dataset), dataset:="Suprapubic Skin"]
dt2p[,dataset:=gsub(" ", "~", dataset)]
dt2p[,type:=factor(type, levels=c("psi5", "psi3", "psiSite"))]
levels(dt2p$type) <- c("psi[5]", "psi[3]", "theta")
dt2p[,injectionDistr:=factor(injectionDistr, levels=c("uniform", "fixed"))]
dt2p[,injectionStartpoint:=factor(injectionStartpoint, levels=c("samplePSI", "meanPSI"))] 
levels(dt2p$injectionStartpoint) <- c("observed", "mean") 
# levels(dt2p$injectionStartpoint) <- paste0(c("sample~", "mean~"), 
                                           # c(psi5="psi[5]", psi3="psi[3]", psiSite="theta")[ptype])
dt2p[,injectedDpsi:=factor(injectedDpsi)]
dt2p[,subsetted:=factor(subsetted)]
dt2p[q != bestQ, bestQ:=FALSE]
dt2p[q == bestQ, bestQ:=TRUE]


dt2p[, col:=paste0(injectionDistr, "(", injectedDpsi, ")")]
dt2p[injectionDistr == "uniform", col:=paste0("uniform(", injectedDpsi, ", 1)")]

# select to plotting data
dt2p <- dt2p[injectionDistr == "fixed" | 
        (injectionDistr == "uniform" & injectedDpsi == 0.2 & injectionStartpoint == "observed")]
cols <- brewer.pal(length(dt2p[,unique(col)]), "Dark2")
cols[5] <- "black"

# dt2p[injectionDistr == "fixed",]
g1 <- ggplot(dt2p, aes(q, aroc, col=col, linetype=injectionStartpoint)) + 
    geom_point(aes(size=factor(bestQ))) + 
    scale_size_manual(name="Optimal q", labels=c(FALSE, TRUE), 
                      values=c(1, 3)) +
    geom_line() + 
    facet_grid(type ~ dataset, labeller=label_parsed, scales="free_y") +
    # scale_x_log10() +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=c(1,2), 
                          labels=parse(text=levels(dt2p$injectionStartpoint))) +
    xlab(expression(Latent~space~dimension~italic(q))) + 
    ylab("Area under the curve\nof precision-recall") + 
    theme_bw() + 
    theme(legend.position="bottom") + # , legend.box="vertical", legend.margin=margin()
    # annotation_logticks(sides="b")  + 
    # labs(linetype="injection from", col=) +
    guides(color   =guide_legend(order=1, ncol=2, title=parse(text="atop(\"injected\", Delta*psi[5]/Delta*psi[3]/Delta*theta)")),
           linetype=guide_legend(order=2, ncol=1, title="injected\nfrom"),
           size    =guide_legend(order=3, ncol=1))
    

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
# g2 <- ggplot(best_qs, aes(q)) + 
#     geom_histogram(stats="identity") + 
#     facet_wrap(~type, labeller=label_parsed) +
#     theme_bw() + 
#     xlab(expression(Latent~space~dimension~italic(q))) + 
#     ylab("Frequency")
# g2

#' 
#' Best Q correlation with sample size
#' 
g3 <- ggplot(best_qs, aes(Nsamples, q)) +
    geom_point() + 
    geom_smooth(method = "lm") + 
    facet_wrap(~type, labeller=label_parsed) + 
    theme_bw() + 
    xlab("Number of samples") + 
    ylab(bquote(atop("Latent space", "dimension" ~ italic(q)))) + 
    ylim(c(0, NA)) + 
    xlim(c(0, NA))
g3



#'
#' Assemble figure
#'
#+ assemble full figure
g <- ggarrange(ncol=1, heights=c(2.8,1), labels=letters[1:2],
    g1, 
    g3)
g

#+ save heatmap figure
factor <- 0.42
outPng
ggsave(outPng, g, width = 17*factor, height = 19*factor)
ggsave(outPdf, g, width = 17*factor, height = 19*factor)

