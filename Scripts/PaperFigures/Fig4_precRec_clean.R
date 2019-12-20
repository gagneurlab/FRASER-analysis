#'---
#' title: Precision-Recall Figure 4
#' author: Ines Scheller, Christian Mertes
#' wb:
#'  input:
#'   - inData: '`sm config["FIGDATADIR"] + "/precRec/{dataset}/{delta}/{outlierType}_precision_recall_{psiType}.RDS"`'
#'   - inCutoffs: '`sm config["FIGDATADIR"] + "/precRec/{dataset}/{delta}/{outlierType}_pr_cutoffs_{psiType}.RDS"`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/Figure4_precRec_clean/{dataset}/{delta}/precRec_{outlierType}_{psiType}_clean.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/Figure4_precRec_clean/{dataset}/{delta}/precRec_{outlierType}_{psiType}_clean.pdf"`'
#'   - wBhtml: 'Output/html/PaperFigures/precRec_clean/{dataset}/{delta}/precRec_{outlierType}_{psiType}_clean.html'
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")
sourceFolder("src/r/precisionRecallHelpers")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic", 
            delta="uniformDistr", outlierType="byJunctionGroup", psiType="psi5")
    parseWBHeader2("Scripts/PaperFigures/Fig4_precRec_clean.R", 
            wildcards=wildcards, rerun=TRUE)
}


#+ input
datain    <- readRDS(snakemake@input$inData)
cutoffsin <- readRDS(snakemake@input$inCutoffs)
dataset   <- snakemake@wildcards$dataset
psiType   <- snakemake@wildcards$psiType
maxRows   <- CONFIG$maxRows
FDR_LIMIT <- CONFIG$FDR_LIMIT
fds_imple <- CONFIG$AE_IMPLEMENTATION
q         <- unique(datain$q)
nPoints4plot <- 30000

#+ output
outPdf <- snakemake@output$outPdf
outPng <- snakemake@output$outPng

#+ filter and reduce points
datain <- datain[correction %in% c("BB", "PCA", fds_imple)]
datain <- datain[type != "DeltaPSI rank"]
datain[, c("precision", "recall", "lower", "upper"):=list(round(precision, 3),
        round(recall, 3), round(lower, 3), round(upper, 3))]
datain <- datain[!duplicated(datain, 
        by=c("type", "precision", "recall", "correction", "junctionMeanBin"))]
datain

#+ correct names
rename_methods <- function(dt){
    levels(dt$inj_value) <- gsub(
        "\\[0\\.000[0-9]+,", "[0.2,", levels(dt[,inj_value]))
    dt[,correction:=factor(correction)]
    levels(dt$correction)[levels(dt[,correction]) == fds_imple] <- "FRASER"
    levels(dt$correction) <- gsub("\\+", "-", levels(dt[,correction]))
    levels(dt$type)[levels(dt$type) == "zScore rank"] <- "\n+ z score"
    levels(dt$type)[levels(dt$type) == 
            "P-value rank & FDR < 0.1 \n& deltaPSI > 0.1"] <- "\n+ P value"
    if(fds_imple == "PCA"){
        dt[correction == "FRASER" & type == "\n+ z score", correction:="PCA"]
    }
    dt
}

data <- rename_methods(datain)
cutoffs <- rename_methods(cutoffsin)

#+ subset data
# Ignore all non matching levels
subset_methods <- function(dt){
    dt[,Method:=factor(paste0(correction, type), levels=c(
        "FRASER\n+ P value", 
        "FRASER\n+ z score", 
        "PCA\n+ z score", 
        "PCA\n+ P value",
        "BB\n+ P value"))]
    dt <- dt[!is.na(Method)]
    if("cutoff" %in% colnames(dt)){
        dt <- dt[order(Method, inj_value, junctionMeanBin, cutoff)]
        dt <- dt[
            cutoff %in% c(0.1) & grepl("P value", Method) | 
                cutoff %in% c(2)      & grepl("z score", Method)]
    } else {
        dt <- dt[order(Method, inj_value, junctionMeanBin, rank)]
    }
    dt
}

data <- subset_methods(data)
cutoffs <- subset_methods(cutoffs)
cutoffs[,Cutoff:=factor(cutoff)]
levels(cutoffs$Cutoff) <- c("FDR < 0.1", "|z score| > 2")

#' 
#' get cutoff straight
#' 
cutoffCorrect <- function(x, data, cutoffs){
    d1 <- cutoffs[x]
    d2 <- data[Method == d1$Method & junctionMeanBin == d1$junctionMeanBin 
               & inj_value == d1$inj_value & recall <= d1$recall]
    d1[,"precision":=d2[nrow(d2), precision]]
    d1[,"recall":=d2[nrow(d2), recall]]
    
    list(d1, d2)
}

res1 <- lapply(seq_row(cutoffs), cutoffCorrect, data=data, cutoffs=cutoffs)
data <- rbindlist(lapply(res1, "[[", 2))
cutoffs <- rbindlist(lapply(res1, "[[", 1))

#'
#' Create figure
#'
#+ create figure
plotdt <- data[data[, rank < max(rank) & (
    rep(c(TRUE, FALSE), c(min(5000, .N), max(0, .N-5000))) |
        sample(c(TRUE, FALSE), .N, replace=TRUE,
               prob=c(nPoints4plot/.N, max(0.1, (.N-nPoints4plot)/.N)))),
    by="Method,inj_value,junctionMeanBin"][,V1]]

g <- ggplot(plotdt, aes(recall, precision, color=Method)) + 
    geom_line() + 
    facet_grid(inj_value ~ junctionMeanBin, labeller=label_bquote(
            rows=atop("Outlier" ~ Delta*psi, .(as.vector(inj_value))),
            cols=atop("Mean coverage", .(as.vector(junctionMeanBin))))) + 
    theme_bw() + 
    theme(legend.position="bottom") + 
    xlab("Recall") + 
    ylab("Precision") + 
    scale_color_brewer(palette="Dark2") + 
    scale_fill_brewer(palette="Dark2") +
    ylim(c(0,1)) + 
    xlim(c(0,1)) + 
    geom_point(data=cutoffs, aes(recall, precision, color=Method, 
            shape=Cutoff), size=3) + 
    theme(strip.text.x=element_text(size=11), 
            strip.text.y=element_text(size=11)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper, 
            fill=Method), alpha=0.2, color=NA)

g


#' 
#' # Save final paper plots
#' 
#+ save plots
fraction <- 0.55
ggsave(outPng, g, width=12*fraction, height=9*fraction)
ggsave(outPdf, g, width=12*fraction, height=9*fraction)


