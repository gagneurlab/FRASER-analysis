#'---
#' title: Paper figure GTEx enrichment
#' author: Christian Mertes
#' wb:
#'  input:
#'    - fds:     '`sm config["DATADIR"] + "/datasets/savedObjects/Kremer__" + config["AE_IMPLEMENTATION"] + "/padjBetaBinomial_psiSite.h5"`'
#'    - rlds:    '`sm config["DATADIR"] + "/processedData/leafcutter/Kremer/rlds_obj.RDS"`'
#'    - kenrich: '`sm "/s/project/abcd-net/processed/fraser_enrichment_analysis/Kremer__" + config["AE_IMPLEMENTATION"] + "_plot_data.Rds"`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/Figure6_Kremer_results.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/Figure6_Kremer_results.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")
library(ggforce)
library(ggrepel)

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/Fig6_Kremer_results.R",
            rerun=TRUE, options=options)
}

#+ input
FraseR_implementation <- snakemake@config$AE_IMPLEMENTATION
FDR_LIMIT       <- 0.1
DELTA_PSI_LIMIT <- 0.3
MIN_READ_COUNT  <- 5
BPPARAM         <- MulticoreParam(20, 100, progre=TRUE)

fdsFile        <- snakemake@input$fds
rldsFile       <- snakemake@input$rlds
kenrichFile    <- snakemake@input$kenrich       
outPng         <- snakemake@output$outPng
outPdf         <- snakemake@output$outPdf

fdsFile
rldsFile
kenrichFile
outPng

#' 
#' Load Kremer dataset
#' 
fds <- loadFraserDataSet(file=fdsFile)
resFds <- bplapply(psiTypes, aberrant, fds=fds, padjCutoff=FDR_LIMIT,
        deltaPsiCutoff=DELTA_PSI_LIMIT, minCoverage=MIN_READ_COUNT, 
        aggregate=TRUE)
names(resFds) <- psiTypes
# rls <- readRDS(rldsFile)
# clusterGeneMap <- clusterGeneMapping(rls)
# res5 <- collectResults(rls, clusterGeneMap=clusterGeneMap,
#         plot=FALSE, FDR_CUTOFF=0.05, BPPARAM=BPPARAM)
# res1 <- collectResults(rls, clusterGeneMap=clusterGeneMap,
#         plot=FALSE, FDR_CUTOFF=FDR_LIMIT, BPPARAM=BPPARAM)
resP <- fread("https://static-content.springer.com/esm/art%3A10.1038%2Fncomms15824/MediaObjects/41467_2017_BFncomms15824_MOESM395_ESM.txt")


#' 
#' Number of events by type in Kremer
#' 
g1tmp <- plotAberrantPerSample(fds, padjCutoff=FDR_LIMIT, zScoreCutoff=NA, 
        deltaPsiCutoff=DELTA_PSI_LIMIT, minCoverage=MIN_READ_COUNT, 
        aggregated=TRUE)

g1yRanges <- layer_scales(g1tmp)$y$range$range
g1xRanges <- layer_scales(g1tmp)$x$range$range
g1 <- g1tmp + 
    theme_bw() + 
    ylab("Number of\naberrantly spliced genes") + 
    theme(legend.position="none", plot.title=element_text(face="bold")) + 
    theme_cowplot() + 
    theme(plot.title = element_blank(), legend.position="right")
g1


#' 
#' Venn diagram of leafcutter / psi / theta
#' 
# res <- res5
# rl_res <- unique(res[!is.na(genes) & gene_fdr < FDR_LIMIT, 
#        .(sampleID=gsub(",.*", "" , sampleIDs), genes, type="Kremer")])
rl_res <- unique(resP[!is.na(HGNCID) & SPLICE_PADJ < FDR_LIMIT, 
        .(sampleID=RNA_ID, genes=HGNCID, type="Kremer")])
fds_res <- unique(rbindlist(lapply(psiTypes, function(type){
        as.data.table(reshape::melt(resFds[[type]]))[
                !is.na(X1) & value == TRUE][,
                .(sampleID=X2, genes=X1, type=ifelse(
                        type == "psiSite", "theta", "psi"))]})))

# compile venn diagramm data
combData <- rbind(rl_res, fds_res)
vd_data <- table(combData[,
        paste0(sort(unique(type)), collapse=","),by="sampleID,genes"][,V1])

# get circles
r <- 1.6
s <- 0.6
df.circles <- data.table(
    x=c(s, -s,  0),
    y=c(s,  s, -s),
    labels=c("theta", "psi", "Kremer~et~al."))
df.coords <- t(cbind(
    c(order=NA,  x=NA,         y=NA, xstart=NA, xend=NA, ystart=NA, yend=NA),
    c(7,            0,   -s - r*0.6, 0, 0, 0, 0),
    c(4, -s + -r*0.30, -s +  r*0.15, 0, 0, 0, 0),
    c(5,            0,          0.1, 0, 0, 0, 0),
    c(6,  s +  r*0.30, -s +  r*0.15, 0, 0, 0, 0),
    c(1, -s + -r*0.55,  s +   r*0.3, 0, 0, 0, 0),
    c(2,            0,   s +  r*0.5, 0, 0, 0, 0),
    c(3,  s +  r*0.55,   s +  r*0.3, 0, 0, 0, 0)))[-1,]
df.coords

df.anno <- cbind(df.coords, data.table(
    names=names(vd_data),
    count=as.character(vd_data))) 
df.anno

gvd <- ggplot(df.circles) +
    geom_circle(aes(x0=x, y0=y, r=r, fill=labels), 
            alpha=.6, size=0.6, colour='grey') +
    coord_fixed() +
    theme_void() +
    scale_fill_brewer(palette="Dark2", 
            labels=setNames(
                    lapply(df.circles$labels, function(i) parse(text=i)), 
                    df.circles$labels)) +
    scale_colour_brewer(palette="Dark2", guide=FALSE) +
    labs(fill = NULL) + 
    annotate("text", x=df.anno$x, y = df.anno$y, label=df.anno$count, size=5) + 
    # MCOLN1
    annotate("text", label="1x MCOLN1",
            x=df.anno[order == 3,x] + r * 0.4, 
            y=df.anno[order == 3,y] + r*0.8) + 
    annotate("segment", 
            x   =df.anno[order == 3,x] + r * 0.4, 
            xend=df.anno[order == 3,x]*1, 
            y   =df.anno[order == 3,y] + r * 0.6, 
            yend=df.anno[order == 3,y]*1.2) + 
    # TIMMDC1
    annotate("text", label="2x TIMMDC1\n1x CLPP\n1x TAZ", hjust=0,
            x=df.anno[order == 4,x] - r * 1.2, 
            y=df.anno[order == 4,y] - r*0.8) + 
    annotate("segment", 
            x   =df.anno[order == 4,x] - r * 0.7, 
            xend=df.anno[order == 4,x]*1.3, 
            y   =df.anno[order == 4,y] - r * 0.4,
            yend=df.anno[order == 4,y]*1.3) + 
    # extra theme data
    theme(legend.position='bottom', legend.text=element_text(size=14), 
            plot.margin = unit(c(1,3,1,3), "lines")) + 
    coord_cartesian(expand=F, clip = "off")
gvd

  
#' 
#' # Enrichment plot by Liang 
#' 
#' Create variant type enrichment plot
#' -> TODO TODO
#' 
pd <- readRDS(kenrichFile)
pt <- pd[!(var_type %in% c("no_variant", "no_rare_variant", "coding", 
                           "coding_other", "intergenic", "non_coding"))]
pt[,var_type:=gsub("Utr", "UTR", toTitleCase(as.character(var_type)))]
pt[aberrant_splicing == "NA", aberrant_splicing:="Not\nexpressed"]
height <- max(pt$ci_h)
totals_var <- unique(pt[, .(var_type, total_var)])

g4 <- ggplot(pt, aes(var_type, prop, fill = aberrant_splicing)) + 
    geom_bar(stat= 'identity', position = 'dodge') + 
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), position=position_dodge(.9), width=.3) + 
    scale_fill_manual(values = c("royalblue", "firebrick", "gray70")) +
    scale_y_continuous(labels=scales::percent_format(accuracy = 1)) + 
    theme(legend.position="top") + 
    geom_text(aes(var_type, height + 0.005, label = total_var, fill = NULL), data = totals_var) +
    theme_bw() +
    labs(x = "Variant type", y = "Proportion", fill = "Aberrant\nsplicing") +
    theme_cowplot() + 
    grids() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
g4



#'
#' Arrange the plots
#'
#g <- ggarrange(nrow=2, labels=c(NA, letters[3]), heights=c(1,1.3),
#    ggarrange(ncol=2, labels=letters[1:2], widths=c(1,1), 
#        g1, 
#        gvd),
#    g4)
g <- ggarrange(ncol=2, labels=letters[1:2], 
    g1, 
    gvd)
g


#+ save figure
factor <- 0.6
height <- 6.5
ggsave(outPng, g, width = 16*factor, height = height*factor)
ggsave(outPdf, g, width = 16*factor, height = height*factor)

