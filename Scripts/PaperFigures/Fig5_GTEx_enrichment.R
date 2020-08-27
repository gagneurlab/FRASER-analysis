#'---
#' title: Paper figure GTEx enrichment
#' author: Christian Mertes
#' wb:
#'  threads: 10
#'  input:
#'    - rareSplice:   '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareSplicing__0.0.RDS", dataset=config["EnrichmentTissues"])`'
#'    - rareMMSplice: '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareMMSplice__0.0.RDS", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/Figure5_GTEx_enrichment.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/Figure5_GTEx_enrichment.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/Fig5_GTEx_enrichment.R",
            rerun=TRUE, options=options)
    threads <- 30
}

#+ input
rareSpliceLs   <- snakemake@input$rareSplice
rareMMSpliceLs <- snakemake@input$rareMMSplice
outPng         <- snakemake@output$outPng
outPdf         <- snakemake@output$outPdf
threads        <- snakemake@threads

FraseR_implementation <- snakemake@config$AE_IMPLEMENTATION
METHODS_2_PLOT <- c("BB_p", "Leafcutter_p", "LeafcutterMD_p", "SPOT_p")
PVALUES_2_PLOT <- c(-3, -5)
FDR_LIMIT      <- 0.1
BPPARAM        <- MulticoreParam(threads, 100, progre=TRUE)

length(rareSpliceLs)
rareSpliceLs[1:5]
rareMMSpliceLs[1:5]
outPng


#'
#' FUNCTIONS
#'
readEnrichmentFiles <- function(file){
    x <- readRDS(file)[['enrich_final_ls']]
    rbindlist(lapply(names(x)[1:6], function(name){
        tissue <- strsplit(name, ": ")[[1]][1]
        score <- strsplit(name, ": ")[[1]][2]

        tmp_ans <- x[[name]]
        tmp_ans[,tissue:=tissue]
        tmp_ans[,score:=gsub("Pval", "P <", score)]
        tmp_ans[,score:=gsub("Zscore", "Z score", score)]
        tmp_ans
    }))
}


#' 
#' read enrichment data
#' and only extract the final enrichment scores
dt2plotSplice <- rbindlist(
        bplapply(rareSpliceLs, readEnrichmentFiles, BPPARAM=BPPARAM))
dt2plotMMSplice <- rbindlist(
    bplapply(rareMMSpliceLs, readEnrichmentFiles, BPPARAM=BPPARAM))

dt2plot <- rbind(
        dt2plotSplice[,snptype:="Splice region"],
        dt2plotMMSplice[,snptype:="MMSplice"])
dt2plot[,snptype:=factor(snptype, levels=c("Splice region", "MMSplice"))]

#' 
#' remove multi tissue outlier calls
#'
dt2plot <- dt2plot[multiCall == 1]

#'
#' Create P value enrichment plot
#'
#' FraseR vs Leafcutter p values
#'
AE_METHOD <- paste0(FraseR_implementation, "_p")
frdt   <- dt2plot[cutoff == TRUE & Method == AE_METHOD,      .(
        FraseR=enrichment, FraseR_min=min.ci, FraseR_max=max.ci, tissue, score, snptype)]
otherdt <- dt2plot[cutoff == TRUE & Method %in% METHODS_2_PLOT, .(
        Method, enrich=enrichment, enrich_min=min.ci, 
        enrich_max=max.ci, tissue, score, snptype)]
otherdt[, Method:=mName4Plot(Method, removeTest=TRUE, AE_Name=AE_METHOD)]
otherdt[, Method:=factor(Method, levels=mName4Plot(METHODS_2_PLOT, removeTest=TRUE, AE_Name=AE_METHOD))]

dt <- merge(frdt, otherdt)
dt <- dt[grepl(paste0("(", paste(PVALUES_2_PLOT, collapse="|"), ")$"), score)]
dt[,score:=gsub(" 1e-", " 10^-", score)]
dt[,score:=gsub("P ", "italic(P)", score)]
dt[,snptype:=factor(gsub(" ", "~", snptype), levels=c("Splice~region", "MMSplice"))]
levels(dt$Method) <- gsub(" ", "~", levels(dt$Method))

g1 <- ggplot(dt[snptype == "Splice~region"], aes(enrich, FraseR)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER)") +
    xlab("Enrichment (other method)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
    scale_x_log10() +
    scale_y_log10()
g1

g2 <- ggplot(dt[snptype == "MMSplice"], aes(enrich, FraseR)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER)") +
    xlab("Enrichment (other method)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
    scale_x_log10() +
    scale_y_log10()
g2

#'
#' Enrichment comparison of MMSplice and splice region
#' 
dt <- dt[order(tissue, score, snptype)]
hist(
    dt[Method == "Kremer~et~al." & snptype == "MMSplice",FraseR] / 
    dt[Method == "Kremer~et~al." & snptype == "Splice~region",FraseR])

#'
#' Arrange the plots
#'
g <- ggarrange(labels=letters[1:2], align="hv", ncol=1,
    g1,
    g2)
g


#+ save figure
factor <- 0.55
outPng
ggsave(outPng, g, width = 16*factor, height = 14*factor)
ggsave(outPdf, g, width = 16*factor, height = 10*factor)

