#'---
#' title: Paper figure S16 Finding Q for different choices of outlier injection distributions
#' author: Ines Scheller
#' wb:
#'  input:
#'   - encDimTable:  '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + 
#'                               "_qEstimation/old_qAuroc_{injDistr}_{injStartpoint}_{dpsi}_subset{subset}_{{psiType}}.tsv", 
#'                               dataset=config["heatmap_tissues"], injDistr=config["injection_distribution"], 
#'                               injStartpoint=config["injection_startpoint"], 
#'                               dpsi=config["injection_dpsi"], subset=config["injection_subsetting"])`'
#'   - fds:   '`sm expand(config["DATADIR"] + "/datasets/savedObjects/{dataset}__" + config["AE_IMPLEMENTATION"] + "/padjBetaBinomial_psiSite.h5", dataset=config["heatmap_tissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS16_qEst_outlier_distr_{psiType}.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS16_qEst_outlier_distr_{psiType}.pdf"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/PaperFigures/qEstimation/qAucpr_{psiType}.html"`'
#'  type: noindex
#'  threads: 5
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list()
    opt <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigureS16_q_estimation_outlierDistr.R",
                   wildcards=wildcards, options=opt, rerun=TRUE)
}


#+ input
fdsFiles    <- snakemake@input$fds
encDimTableFiles  <- snakemake@input$encDimTable
ptype <- snakemake@wildcards$psiType
datasets    <- basename(dirname(dirname(encDimTableFiles)))
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

#'
#' Extract q estimation data
#' 
message("fds to load: ", paste0(unique(tissues), "__", CONFIG$AE_IMPLEMENTATION))
fds_ls <- bplapply(paste0(unique(tissues), "__", CONFIG$AE_IMPLEMENTATION), loadFraserDataSet, dir=workingDir, BPPARAM=BPPARAM)
names(fds_ls) <- unique(tissues)

dt2p <- rbindlist(lapply(encDimTableFiles, fread))
dt2p <- dt2p[type == ptype,]

dt2p_fds <- rbindlist(lapply(names(fds_ls), function(name){
    dttmp <- hyperParams(fds_ls[[name]], type=ptype, all=TRUE)
    dttmp[,dataset:=name]
    dttmp[,method:=CONFIG$AE_IMPLEMENTATION]
    dttmp[,type:=ptype]
    dttmp[,bestQ:=q[which.max(aroc)]]
    dttmp[,injectionDistr:="uniform"]
    dttmp[,injectionStartpoint:="samplePSI"]
    dttmp[,injectedDpsi:=0.2]
    dttmp[,subsetted:=TRUE]
    dttmp 
}))

dt2p <- rbind(dt2p, dt2p_fds)

# better naming
dt2p[,dataset:=dName4plot(dataset, TRUE)]
dt2p[grepl("Skin not Sun Exposed", dataset), dataset:="Suprapubic Skin"]
dt2p[,dataset:=gsub(" ", "~", dataset)]
# dt2p[,type:=factor(type, levels=c("psi5", "psi3", "psiSite"))]
levels(dt2p$type) <- c("psi[5]", "psi[3]", "theta")
dt2p[,injectionDistr:=factor(injectionDistr, levels=c("uniform", "fixed"))]
dt2p[,injectionStartpoint:=factor(injectionStartpoint, levels=c("meanPSI", "samplePSI"))] 
levels(dt2p$injectionStartpoint) <- paste0(c("mean~", "sample~"), 
                    c(psi5="psi[5]", psi3="psi[3]", psiSite="theta")[ptype])
dt2p[,injectedDpsi:=factor(injectedDpsi)]
dt2p[,subsetted:=factor(subsetted)]
dt2p[q != bestQ, bestQ:=FALSE]
dt2p[q == bestQ, bestQ:=TRUE]

# g <- ggplot(dt2p, aes(q, aroc, col=injectedDpsi, linetype=injectionDistr)) + 
#     geom_point(aes(size=factor(bestQ))) + 
#     scale_size_manual(name="Optimal q", labels=c(FALSE, TRUE), 
#                       values=c(1, 3)) +
#     geom_line() + 
#     facet_grid(type ~ dataset, labeller=label_parsed) +
#     scale_x_log10() +
#     scale_color_brewer(palette="Dark2") +
#     xlab(expression(Latent~space~dimension~italic(q))) + 
#     ylab("Area under the curve\nof precision-recall") + 
#     theme_bw() + 
#     theme(legend.position="bottom")
# g

# g1 <- ggplot(dt2p, aes(q, aroc, col=injectedDpsi, linetype=injectionDistr)) + 
#     geom_point(aes(size=factor(bestQ))) + 
#     scale_size_manual(name="Optimal q", labels=c(FALSE, TRUE), 
#                       values=c(1, 3)) +
#     geom_line() + 
#     facet_grid(type ~ injectionStartpoint, labeller=label_parsed) +
#     scale_x_log10() +
#     scale_color_brewer(palette="Dark2") +
#     xlab(expression(Latent~space~dimension~italic(q))) + 
#     ylab("Area under the curve\nof precision-recall") + 
#     theme_bw() + 
#     theme(legend.position="bottom") +
#     annotation_logticks(sides="b")   
# g1

g2 <- ggplot(dt2p, aes(q, aroc, col=injectedDpsi, linetype=injectionStartpoint)) + 
    geom_point(aes(size=factor(bestQ))) + 
    scale_size_manual(name="Optimal q", labels=c(FALSE, TRUE), 
                        values=c(1, 3)) +
    geom_line() + 
    facet_grid(injectionDistr ~ dataset, labeller=label_parsed) +
    labs(linetype="injection from", col=parse(text="(minimal)~injected~Delta*psi[5]")) +
    # scale_x_log10() +
    scale_color_brewer(palette="Dark2") +
    scale_linetype_manual(values=c(1,2), 
            labels=parse(text=levels(dt2p$injectionStartpoint))) +
    xlab(expression(Latent~space~dimension~italic(q))) + 
    ylab("Area under the curve\nof precision-recall") + 
    theme_bw() + 
    theme(legend.position="bottom") +
    # annotation_logticks(sides="b")  + 
    guides(color=guide_legend(order=1), linetype=guide_legend(order=2), 
            size=guide_legend(order=3))
g2


#+ save figure
factor <- 0.6
outPng
ggsave(outPng, g2, width = 17*factor, height = 14*factor)
ggsave(outPdf, g2, width = 17*factor, height = 14*factor)

