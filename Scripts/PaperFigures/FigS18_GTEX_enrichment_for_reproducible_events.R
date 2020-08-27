#'---
#' title: Supplemental GTEx exnrichment for replicated outliers (S18)
#' author: Christian Mertes
#' wb:
#'  threads: 16
#'  input:
#'    - rareSplice0:    '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareSplicing__0.0.RDS", dataset=config["EnrichmentTissues"])`'
#'    - rareMMSplice0:  '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareMMSplice__0.0.RDS", dataset=config["EnrichmentTissues"])`'
#'    - rareSplice1:    '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareSplicing__0.1.RDS", dataset=config["EnrichmentTissues"])`'
#'    - rareMMSplice1:  '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__rareMMSplice__0.1.RDS", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPPng: '`sm config["FIGDIR"] + "/FigureS18_GTEx_enrichment_for_reproducible_events.png"`'
#'   - outPPdf: '`sm config["FIGDIR"] + "/FigureS18_GTEx_enrichment_for_reproducible_events.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS18_GTEX_enrichment_for_reproducible_events.R",
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
methodsOfInteres <- unique(c(paste0(FraseR_implementation, "_p"), "BB_p",
        # "Leafcutter_p", 
        "LeafcutterMD_p", "PCA_p", "SPOT_p"))
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
dt2plot <- dt2plot[Method != "BB_np" & !grepl("_only_", Method)]
dt2plot <- dt2plot[!(Method == "Leafcutter_p" & dcut == "0.0")]
dt2plot[,multiCall:=factor(multiCall)]
dt2plot

#' 
#' remove methods without reproducible data
#' 
dt2plot <- dt2plot[Method %in% methodsOfInteres]
tibble(dt2plot)
table(dt2plot$Method)

#'
#' Create P value enrichment plot
#'
#' FraseR vs Leafcutter p values
#'
AE_METHOD <- paste0(FraseR_implementation, "_p")
frdt   <- dt2plot[cutoff == TRUE & Method == AE_METHOD, .(
    FRASER=enrichment, FRASER_min=min.ci, FRASER_max=max.ci,
    tissue, score, snptype, dcut, multiCall)]

tmpdt2p <- merge(dt2plot[cutoff == TRUE & grepl("P < ", score) & 
                             Method != AE_METHOD & !grepl("no-weights", Method)], frdt)
tmpdt2p[,pvalType:=ifelse(grepl("_np", Method), "Gaussian P", "BetaBinom P")]
tmpdt2p[,dcut:=paste0("|dPSI| > ", dcut)]
tmpdt2p[grepl("Leafcutter|SPOT", Method), pvalType:="DM P"]
tmpdt2p[grepl("Leafcutter_|SPOT_", Method), dcut:=""]
tmpdt2p[,Method:=mName4Plot(Method, removeTest=TRUE, AE_METHOD)]
tmpdt2p[Method == "BetaBinom", Method:="Naïve"]
tmpdt2p[Method == "BB-AE-weight", Method:="BB-regression"]
tmpdt2p[,Method:=factor(Method, levels=mName4Plot(methodsOfInteres, AE_Name = AE_METHOD))]

tmpdt2p

#' 
#' ggPlotting theming
#' 
numMulticalls <- length(unique(tmpdt2p$multiCall))
figure_theme <- function(){
    list(
        labs(
            y="Enrichment (FRASER)",
            x="Enrichment (other method)",
            color="Minimum of\nreproducibility"),
        geom_abline(slope=1, intercept=0),
        geom_point(aes(x=1,y=1), col="white", alpha=0),
        scale_color_manual(values=colorRampPalette(c("orange", "firebrick"))(numMulticalls)), 
        theme_cowplot(),
        grids(),
        facet_grid(facets=score ~ Method + pvalType, scales="free"),
        scale_x_log10(),
        scale_y_log10()
    )
}

#' 
#' rare splicing with additional cutoff
#' 
tdata <- tmpdt2p[snptype == "Splice region" & grepl("> 0.1", dcut)]
gp1 <- ggplot(tdata, aes(enrichment, FRASER, col=multiCall)) +
    geom_point(alpha=0.7) +
    figure_theme()
gp1


#'
#' rare splicing and no additional cutoff
#' 
tdata <- tmpdt2p[snptype == "Splice region" & !grepl("> 0.1", dcut)]
gp2 <- ggplot(tdata, aes(enrichment, FRASER, col=multiCall)) +
    geom_point(alpha=0.7) + 
    figure_theme()
gp2

#'
#' rare MMSplice and additional cutoff
#' 
tdata <- tmpdt2p[snptype == "MMSplice" & grepl("> 0.1", dcut)]
gp3 <- ggplot(tdata, aes(enrichment, FRASER, col=multiCall)) +
    geom_point(alpha=0.7) +
    figure_theme()
gp3

#'
#' rare MMSplice and no additional cutoff
#' 
tdata <- tmpdt2p[snptype == "MMSplice" & !grepl("> 0.1", dcut)]
gp4 <- ggplot(tdata, aes(enrichment, FRASER, col=multiCall)) +
    geom_point(alpha=0.7) +
    figure_theme()
gp4

#'
#' Assemble figure
#'
#+ assemble full figure
g <- ggarrange(nrow=3, labels=c(letters[1:2], "", ""), 
            common.legend=TRUE, legend="right",
    gp2 + ggtitle(expression(paste("Splice variants with no additional ", Delta, psi, " cutoff"))),
    gp4 + ggtitle(expression(paste("MMSplice variants with no additional ", Delta, psi, " cutoff"))),
    ggarrange(ncol=2, labels=letters[3:4], legend="none", label.x=c(0, 0.05),
        gp1 + ggtitle(expression(paste("Splice variants with a |", Delta, psi, "|" > 0.1~"cutoff"))),
        gp3 + ggtitle(expression(paste("MMSplice variants with a |", Delta, psi, "|" > 0.1~"cutoff")))))
g

#+ save p value enrichment figure
factor <- 0.6
outPPng
ggsave(outPPng, g, width = 16*factor, height = 19*factor)
ggsave(outPPdf, g, width = 16*factor, height = 19*factor)









