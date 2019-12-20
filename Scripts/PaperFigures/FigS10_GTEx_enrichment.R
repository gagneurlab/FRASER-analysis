#'---
#' title: Supplemental GTEx exnrichment
#' author: Christian Mertes
#' wb:
#'  threads: 16
#'  input:
#'    - rareSplice0:    '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareSplicing__0.0.RDS", dataset=config["EnrichmentTissues"])`'
#'    - rareMMSplice0:  '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareMMSplice__0.0.RDS", dataset=config["EnrichmentTissues"])`'
#'    - rareSplice1:    '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareSplicing__0.1.RDS", dataset=config["EnrichmentTissues"])`'
#'    - rareMMSplice1:  '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareMMSplice__0.1.RDS", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPPng: '`sm config["FIGDIR"] + "/FigureS10_GTEx_enrichment_pval.png"`'
#'   - outPPdf: '`sm config["FIGDIR"] + "/FigureS10_GTEx_enrichment_pval.pdf"`'
#'   - outZPng: '`sm config["FIGDIR"] + "/FigureS11_GTEx_enrichment_zsco.png"`'
#'   - outZPdf: '`sm config["FIGDIR"] + "/FigureS11_GTEx_enrichment_zsco.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS10_GTEx_enrichment.R",
            options=options, rerun=TRUE)
    threads <- 25
}

#+ input
rareSplice0Ls   <- snakemake@input$rareSplice0
rareMMSplice0Ls <- snakemake@input$rareMMSplice0
rareSplice1Ls   <- snakemake@input$rareSplice1
rareMMSplice1Ls <- snakemake@input$rareMMSplice1

outPPng         <- snakemake@output$outPPng
outPPdf         <- snakemake@output$outPPdf
outZPng         <- snakemake@output$outZPng
outZPdf         <- snakemake@output$outZPdf
threads         <- snakemake@threads

FraseR_implementation <- snakemake@config$AE_IMPLEMENTATION
BPPARAM <- MulticoreParam(threads, 100, progre=TRUE)

length(rareSplice0Ls)
rareSplice0Ls[1:5]
rareSplice1Ls[1:5]
outPPdf

#'
#' FUNCTIONS
#'
readEnrichmentFiles <- function(file){
    x <- readRDS(file)[['enrich_final_ls']]
    snptypev <- strsplit(file, "__")[[1]][2]
    dcutv    <- gsub(".RDS$", "", strsplit(file, "__")[[1]][3])
    
    rbindlist(lapply(names(x)[1:6], function(name){
        tissue <- strsplit(name, ": ")[[1]][1]
        score <- strsplit(name, ": ")[[1]][2]
        
        tmp_ans <- x[[name]]
        tmp_ans[,tissue:=tissue]
        tmp_ans[,score:=gsub("Pval", "P <", score)]
        tmp_ans[,score:=gsub("Zscore", "Z score", score)]
        tmp_ans[,dcut:=dcutv]
        tmp_ans[,snptype:=snptypev]
        tmp_ans
    }))
}


#' read enrichment data
#' and only extract the final enrichment scores
dt2plot <- rbindlist(bplapply(
        c(rareSplice0Ls, rareSplice1Ls, rareMMSplice0Ls, rareMMSplice1Ls),
        readEnrichmentFiles, BPPARAM=BPPARAM))
dt2plot[snptype == "rareSplicing", snptype:="Splice region"]
dt2plot[snptype == "rareMMSplice", snptype:="MMSplice"]
dt2plot[,snptype:=factor(snptype, levels=c("Splice region", "MMSplice"))]
dt2plot <- dt2plot[Method != "BB_np"]
dt2plot <- dt2plot[!(Method == "Leafcutter_p" & dcut == "0.0")]
dt2plot

#'
#' Create P value enrichment plot
#'
#' FraseR vs Leafcutter p values
#'
AE_METHOD <- paste0(FraseR_implementation, "_p")
frdt   <- dt2plot[cutoff == TRUE & Method == AE_METHOD, .(
        FRASER=enrichment, FRASER_min=min.ci, FRASER_max=max.ci,
        tissue, score, snptype, dcut)]

tmpdt2p <- merge(dt2plot[cutoff == TRUE & grepl("P < ", score) & 
        Method != AE_METHOD & !grepl("no-weights", Method)], frdt)
tmpdt2p[,pvalType:=ifelse(grepl("_np", Method), "Gaussian P", "BetaBinom P")]
tmpdt2p[,dcut:=paste0("|dPSI| > ", dcut)]
tmpdt2p[grepl("Leafcutter", Method), pvalType:="DM P"]
tmpdt2p[grepl("Leafcutter", Method), dcut:=""]
tmpdt2p[,Method:=mName4Plot(Method, removeTest=TRUE, AE_METHOD)]
tmpdt2p[Method == "BetaBinom", Method:="Naïve"]
tmpdt2p[Method == "BB-AE-weight", Method:="BB-regression"]

tmpdt2p

gp1 <- ggplot(tmpdt2p[snptype == "Splice region" & grepl("> 0.1", dcut)], aes(enrichment, FRASER)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER)") +
    xlab("Enrichment (other method)") +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score ~ Method + pvalType, scales="free") +
    scale_x_log10() +
    scale_y_log10()
gp1

gp2 <- ggplot(tmpdt2p[snptype == "Splice region" & !grepl("> 0.1", dcut)], aes(enrichment, FRASER)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER)") +
    xlab("Enrichment (other method)") +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score ~ Method + pvalType, scales="free") +
    scale_x_log10() +
    scale_y_log10()
gp2

gp3 <- ggplot(tmpdt2p[snptype == "MMSplice" & grepl("> 0.1", dcut)], aes(enrichment, FRASER)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER)") +
    xlab("Enrichment (other method)") +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score ~ Method + pvalType, scales="free") +
    scale_x_log10() +
    scale_y_log10()
gp3

gp4 <- ggplot(tmpdt2p[snptype == "MMSplice" & !grepl("> 0.1", dcut)], aes(enrichment, FRASER)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER)") +
    xlab("Enrichment (other method)") +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score ~ Method + pvalType, scales="free") +
    scale_x_log10() +
    scale_y_log10()
gp4

#'
#' Create Zscore enrichment plot
#'
#' FraseR vs PCA z scores
#'
AE_METHOD <- paste0(FraseR_implementation, "_z")
frdt   <- dt2plot[cutoff == TRUE & Method == AE_METHOD, .(
    FRASER=enrichment, FRASER_min=min.ci, FRASER_max=max.ci,
    tissue, score, snptype, dcut)]

tmpdt2p <- merge(dt2plot[cutoff == TRUE & grepl("Z score ", score) & 
        Method != AE_METHOD & !grepl("no-weights", Method)], frdt)
tmpdt2p[,Method:=mName4Plot(Method, removeTest=TRUE, AE_METHOD)]
tmpdt2p[,dcut:=paste0("|dPSI| > ", dcut)]
tmpdt2p[Method == "BetaBinom", Method:="Naïve"]
tmpdt2p

gz <- ggplot(tmpdt2p, aes(enrichment, FRASER)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("FRASER") +
    xlab("Enrichment") +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score ~ snptype + Method + dcut, scales="free") +
    scale_x_log10() +
    scale_y_log10()
gz


#'
#' Assemble figure
#'
#+ assemble full figure
g <- ggarrange(ncol=2, nrow=2, labels=LETTERS[1:4], align="hv", widths=c(5,4),
    gp2 + ggtitle(expression(paste("Splice variants with no additional ", Delta, psi, " cutoff"))),
    gp1 + ggtitle(expression(paste("Splice variants with a |", Delta, psi, "|" > 0.1~"cutoff"))),
    gp4 + ggtitle(expression(paste("MMSplice variants with no additional ", Delta, psi, " cutoff"))),
    gp3 + ggtitle(expression(paste("MMSplice variants with a |", Delta, psi, "|" > 0.1~"cutoff"))))

g

#+ save p value enrichment figure
factor <- 0.8
outPPng
ggsave(outPPng, g, width = 17*factor, height = 14*factor)
ggsave(outPPdf, g, width = 17*factor, height = 14*factor)


#+ save z score enrichment figure
factor <- 0.8
ggsave(outZPng, gz, width = 16*factor, height = 11*factor)
ggsave(outZPdf, gz, width = 16*factor, height = 11*factor)









