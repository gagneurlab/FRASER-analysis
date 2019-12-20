#'---
#' title: Paper figure (Precision-Recall)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - inData: '`sm config["FIGDATADIR"] + "/precRec/{dataset}/{delta}/{outlierType}_precision_recall_{psiType}.RDS"`'
#'   - inCutoffs: '`sm config["FIGDATADIR"] + "/precRec/{dataset}/{delta}/{outlierType}_pr_cutoffs_{psiType}.RDS"`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/Figure4_precRec/{dataset}/{delta}/precRec_{outlierType}_{psiType}.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/Figure4_precRec/{dataset}/{delta}/precRec_{outlierType}_{psiType}.pdf"`'
#'   - wBhtml: 'Output/html/PaperFigures/precRec/{dataset}/{delta}/precRec_{outlierType}_{psiType}.html'
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")
sourceFolder("src/r/precisionRecallHelpers")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic", 
            delta="uniformDistr", outlierType="byJunctionGroup", 
            psiType="psi3")
    parseWBHeader2("Scripts/PaperFigures/Fig4_precRec.R", 
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

#' Input
outPng
dataset
snakemake@input$inData

#+ filter
datain <- datain[type != "DeltaPSI rank"]
datain[, c("precision", "recall", "lower", "upper"):=list(round(precision, 3),
        round(recall, 3), round(lower, 3), round(upper, 3))]
datain <- datain[!duplicated(datain, 
        by=c("type", "precision", "recall", "correction", "junctionMeanBin"))]
datain

#+ correct names
rename_methods <- function(dt){
    levels(dt$inj_value)<- gsub(
            "\\[0\\.000[0-9]+,", "[0.2,", levels(dt[,inj_value]))
    dt[correction == "BB", correction:="naïve BB"]
    dt[,correction:=factor(correction)]
    levels(dt$correction)[levels(dt[,correction]) == fds_imple] <- "FRASER"
    levels(dt$correction) <- gsub("\\+", "-", levels(dt[,correction]))
    levels(dt$type)[levels(dt$type) == "zScore rank"] <- "\n+ z score"
    levels(dt$type)[levels(dt$type) == 
            "P-value rank & FDR < 0.1 \n& deltaPSI > 0.1"] <- "\n+ P value"
    if(fds_imple == "PCA"){
        dt[correction == "FRASER" & type == "\n+ z score", correction:="PCA"]
    }
    if("PCA-BB-Decoder" %in% levels(dt$correction)){
        levels(dt$correction)[levels(dt$correction) == "PCA-BB-Decoder"] <- 
                "PCA + BB-Decoder"
    }
    dt
}

data <- rename_methods(datain)
cutoffs <- rename_methods(cutoffsin)

#+ subset data
# Ignore all non matching levels
subset_methods <- function(dt){
    dt[,Method:=factor(paste0(correction, type), levels=c(
            "PCA + BB-Decoder\n+ P value",
            "PCA + BB-Decoder\n+ z score",
            "FRASER\n+ P value",
            "FRASER\n+ z score", 
            "PCA\n+ z score", 
            "PCA\n+ P value",
            "naïve BB\n+ P value"))]
    dt <- dt[!is.na(Method)]
    if("cutoff" %in% colnames(dt)){
        dt <- dt[order(Method, inj_value, junctionMeanBin, cutoff)]
        dt <- dt[
            cutoff %in% c(0.1, 0.05) & grepl("P value", Method) | 
            cutoff %in% c(2, 3)      & grepl("z score", Method)]
    } else {
        dt <- dt[order(Method, inj_value, junctionMeanBin, rank)]
    }
    dt
}

data <- subset_methods(data)
cutoffs <- subset_methods(cutoffs)
cutoffs[,Cutoff:=factor(cutoff)]

#' Adjust lower end
dt <- rbind(data, cutoffs, fill=TRUE)
dt <- dt[,item:=paste(Method, inj_value, junctionMeanBin)]
i <- unique(dt$item)[10]
for(i in unique(dt$item)){
    minVal <- dt[item == i & (grepl("z score", Method) & cutoff == 2 | 
            grepl("P value", Method) & cutoff == 0.1)]
    dt[item == i & recall > minVal$recall, recall:=NA]
    cut <- 0.2
    tmp <- dt[item == i & recall > 0.4 & precision < cut & !is.na(lower)]
    idx <- min(max(1, max(which(diff(tmp[,precision]) < 0))), which(tmp[, precision < 0.02]))
    cut <- tmp[idx,recall]
    if(nrow(tmp) == 0){
        next
    }
    tmp <- tmp[idx]
    getPred <- function(re){
        d <- minVal$recall-tmp$recall
        m <- (minVal$precision-tmp$precision)/d
        tmp$precision + pmin(m*(re - tmp$recall), 0)
    }
    dt[item == i & recall > cut, precision:=getPred(recall)]
    dt[item == i & recall > cut, upper:=getPred(recall)]
    dt[item == i & recall > cut, lower:=getPred(recall)]
}


#'
#' Create figure
#'
#+ create figure
plotdt <- dt[!is.na(recall)]
plotCut <- cutoffs
levels(plotdt$Method) <- gsub("Decoder", "regression", levels(plotdt$Method))
levels(plotCut$Method) <- gsub("Decoder", "regression", levels(plotCut$Method))
lvOrder <- c("FRASER\n+ P value", "PCA + BB-regression\n+ P value", 
        "naïve BB\n+ P value", "PCA + BB-regression\n+ z score", "PCA\n+ z score")
plotdt[ ,Method:=factor(Method, levels=lvOrder)]
plotCut[,Method:=factor(Method, levels=lvOrder)]

g <- ggplot(plotdt, aes(recall, precision, color=Method)) + 
    geom_line() + 
    theme_bw() + 
    theme(legend.position="bottom") + 
    xlab("Recall") + 
    ylab("Precision") + 
    scale_color_brewer(palette="Dark2") + 
    scale_fill_brewer(palette="Dark2") +
    ylim(c(0,1)) + 
    xlim(c(0,1)) + 
    geom_point(data=plotCut, aes(recall, precision, color=Method, 
            shape=Cutoff), size=3) + 
    theme(strip.text.x=element_text(size=11), 
            strip.text.y=element_text(size=11)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper, 
            fill=Method), alpha=0.2, color=NA) + 
    guides(
        color=guide_legend(nrow=2, byrow=TRUE),
        shape=guide_legend(nrow=2, byrow=TRUE))

if(psiType == "psiSite"){
    g <- g + facet_grid(inj_value ~ junctionMeanBin, labeller=label_bquote(
        rows=atop("Outlier" ~ Delta*theta, .(as.vector(inj_value))),
        cols=atop("Mean coverage", .(as.vector(junctionMeanBin)))))
} else { 
    g <- g + facet_grid(inj_value ~ junctionMeanBin, labeller=label_bquote(
        rows=atop("Outlier" ~ Delta*psi, .(as.vector(inj_value))),
        cols=atop("Mean coverage", .(as.vector(junctionMeanBin)))))
}
g


#' 
#' # Save final paper plots
#' 
#+ save plots
factor <- 0.55
outPng
ggsave(outPng, g, width=12*factor, height=10*factor)
ggsave(outPdf, g, width=12*factor, height=8*factor)

