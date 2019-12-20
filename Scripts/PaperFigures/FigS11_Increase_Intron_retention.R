#'---
#' title: Supplemental Increase of Intron retention
#' author: Christian Mertes
#' wb:
#'  threads: 16
#'  input:
#'     - annos:  '`sm expand(config["DATADIR"] + "/annotations/{dataset}.tsv", dataset=config["EnrichmentTissues"])`'
#'     - rlres:  '`sm expand(config["DATADIR"] + "/processedData/leafcutter/{dataset}/results_{dataset}.tsv", dataset=config["EnrichmentTissues"])`'
#'     - fdsres: '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_results.tsv", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS11_IntronRetention.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS11_IntronRetention.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    parseWBHeader2("Scripts/PaperFigures/FigS11_Increase_Intron_retention.R",
            rerun=TRUE)
    threads <- 16
}

#+ input
annoFiles   <- snakemake@input$annos
rlresFiles  <- snakemake@input$rlres
fdsresFiles <- snakemake@input$fdsres

outPng      <- snakemake@output$outPng
outPdf      <- snakemake@output$outPdf

FDR_LIMIT   <- 0.1
MIN_READ_COUNT <- 5
DELTA_PSI_LIMIT <- 0.3
threads     <- snakemake@threads
BPPARAM     <- MulticoreParam(threads, 100, progre=TRUE)
register(BPPARAM)

length(rlresFiles)
rlresFiles[1:5]
fdsresFiles[1:5]
outPdf

#'
#' FUNCTIONS
#'
loadResL <- function(f){
    dt <- fread(f)[gene_p < FDR_LIMIT]
    unique(dt[,.(gene=genes, sampleID=conditionID, type="Leafcutter", 
            tissue=basename(dirname(f)))])
}

loadResF <- function(f){
    dt <- fread(f)[padjust < FDR_LIMIT & totalCounts >= MIN_READ_COUNT &
            abs(deltaPsi) > DELTA_PSI_LIMIT]
    unique(dt[,.(gene=hgncSymbol, tissue=basename(dirname(f)),
            sampleID=indivID, type=ifelse(type != "psiSite", "psi", type))])
}

fdsData <- mclapply(fdsresFiles, loadResF, mc.cores=4)
rlData  <- mclapply(rlresFiles,  loadResL, mc.cores=4)
anno    <- rbindlist(lapply(annoFiles, fread))[,.N,by=SMTSD][,.(N, tissue=SMTSD)]

#'
#' Percentage of intron retention and overlap
#' 
dt2p <- rbindlist(mclapply(fdsData, mc.cores=4, function(x) {
    x[, 
        paste(unique(sort(type)), collapse=","), by="gene,sampleID,tissue"][,
        .(percentage=as.vector(table(V1)/.N), n=.N, type=names(table(V1))), by="tissue"] }))
dt2p[,type:=as.factor(type)]
levels(dt2p$type) <- c("psi~only", "both", "theta~only")
plabels <- setNames(nm=levels(dt2p$type),
        object=lapply(levels(dt2p$type), function(i) parse(text=i))) 

dt2p <- merge(dt2p, anno)
g1 <- ggplot(dt2p, aes(x=reorder(tissue, -N), y=percentage*100, fill=type)) +
    geom_bar(stat="identity") + 
    theme_cowplot() + 
    ylab("Percentage") + 
    xlab("Tissue") + 
    labs(fill="Detection") + 
    scale_fill_brewer(palette="Dark2", labels=plabels) +
    theme(
            axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
            axis.text   = element_text(size=11),
            axis.title  = element_text(face="bold", size=14),
            legend.text = element_text(size=14))
g1


#'
#' Total number of events per psi/theta/overlap
#' 
g2 <- ggplot(dt2p, aes(x=reorder(tissue, -N), y=percentage*n, fill=type)) +
    geom_bar(stat="identity") + 
    theme_cowplot() + 
    ylab("Number of events") + 
    xlab("Tissue") + 
    labs(fill="Detection") + 
    scale_fill_brewer(palette="Dark2", labels=plabels) +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.text   = element_text(size=11),
        axis.title  = element_text(face="bold", size=14),
        legend.text = element_text(size=14))
g2


#'
#' Overlap with leafcutter
#'
tmp_dt2p <- rbindlist(c(rlData, fdsData), use.names=TRUE)
tmp_dt2p[type != "Leafcutter", type:="FRASER"]
dt2p <- rbindlist(mclapply(tmp_dt2p[,unique(tissue)], mc.cores=10, function(x){
    tmp_dt2p[tissue==x, 
        paste(unique(sort(type)), collapse=","), by="gene,sampleID,tissue"][,
        .(percentage=as.vector(table(V1)/.N), n=.N, type=names(table(V1))), by="tissue"] }))

dt2p[,type:=as.factor(type)]
levels(dt2p$type) <- c("FRASER~only", "both", "Leafcutter~only")
plabels <- setNames(nm=levels(dt2p$type),
        object=lapply(levels(dt2p$type), function(i) parse(text=i))) 

dt2p <- merge(dt2p, anno)
g3 <- ggplot(dt2p, aes(x=reorder(tissue, -N), y=percentage, fill=type)) +
    geom_bar(stat="identity") + 
    theme_cowplot() + 
    ylab("Percentage") + 
    xlab("Tissue") + 
    labs(fill="Detection") + 
    scale_fill_brewer(palette="Set2", labels=plabels) +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.text   = element_text(size=11),
        axis.title  = element_text(face="bold", size=14),
        legend.text = element_text(size=14))
g3


#'
#' Overlap with leafcutter
#'
g4 <- ggplot(dt2p, aes(x=reorder(tissue, -N), y=percentage*n, fill=type)) +
    geom_bar(stat="identity") + 
    theme_cowplot() + 
    ylab("Number of events") + 
    xlab("Tissue") + 
    labs(fill="Detection") + 
    scale_fill_brewer(palette="Set2", labels=plabels) +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.text   = element_text(size=11),
        axis.title  = element_text(face="bold", size=14),
        legend.text = element_text(size=14))
g4



#'
#' Arrange the plots
#'
gf1 <- ggarrange(nrow=2, heights=c(20,1),
    ggarrange(ncol=4, labels=LETTERS[1:4], widths=c(2.5,1,1,1), legend="none",
        align="h", label.x=c(0, -0.05, -0.05, -0.05),
        g1 + coord_flip(),
        g2 + theme(axis.text.y=element_blank(), axis.title.y=element_blank()) + coord_flip(),
        g3 + theme(axis.text.y=element_blank(), axis.title.y=element_blank()) + coord_flip(),
        g4 + theme(axis.text.y=element_blank(), axis.title.y=element_blank()) + coord_flip()),
    ggarrange(ncol=3, widths=c(2,5,5),
        ggplot() + theme_void(),
        get_legend(g1 + theme(legend.position="bottom")), 
        get_legend(g3 + theme(legend.position="bottom"))))
gf1

gf2 <- ggarrange(ncol=2, labels=LETTERS[1:2], widths=c(2,1), 
        legend="right", common.legend=TRUE,
        align="h", label.x=c(0, -0.05),
        g1 + coord_flip(),
        g2 + theme(axis.text.y=element_blank(), axis.title.y=element_blank()) + coord_flip())
gf2


#+ save figure
factor <- 0.7
ggsave(outPng, gf2, width = 13*factor, height = 16*factor)
ggsave(outPdf, gf2, width = 13*factor, height = 16*factor)






