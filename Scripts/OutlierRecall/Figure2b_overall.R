#'---
#' title: Figure 2b recall overall
#' author: Ines Scheller
#' wb:
#'   threads: 12
#'   input:
#'     - evalTable: '`sm expand("Output/html/OutlierInjection/{{dataset}}/{{psiType}}/{{delta}}/{method}_pvalues.html", method=config["methods"])`'
#'   output:
#'     - figurePdf: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_overall.pdf"`'
#'     - figurePng: '`sm config["FIGDIR"] + "/Figure2b/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_overall.png"`'
#'     - ggplot:    '`sm config["DATADIR"] + "/processedData/precRec/{dataset}/inject_{delta}/{psiType}/{outlierType}_outlier_recall_overall_ggplot.RDS"`'
#'     - wBhtml:    'Output/html/OutlierRecall/{dataset}/{delta}/{psiType}/{outlierType}_outlier_recall_overall.html'
#'   type: noindex
#'---

#+ load config
source('./src/r/config.R')
sourceFolder("src/r/precisionRecallHelpers")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic", 
            delta="uniformDistr", psiType="psi5", outlierType="byJunctionGroup")
    parseWBHeader2("Scripts/OutlierRecall/Figure2b_overall.R",
            wildcards=wildcards, rerun=TRUE)
    threads <- 16
}

#+ input
dataset     <- snakemake@wildcards$dataset
inj_delta   <- snakemake@wildcards$delta
psiType     <- snakemake@wildcards$psiType
outlierType <- snakemake@wildcards$outlierType
methods     <- CONFIG$methods
baseDir     <- file.path(CONFIG$DATADIR, "datasets", paste0("inject_", inj_delta))
FDR_LIMIT   <- CONFIG$FDR_LIMIT
maxRows     <- CONFIG$maxRows

FDR_LIMIT       <- 0.1
DELTA_PSI_LIMIT <- 0.0

#+ output
ggplotFile <- snakemake@output$ggplot
outPdf <- snakemake@output$figurePdf
outPng <- snakemake@output$figurePng

dataset
methods
outPdf
ggplotFile


#' 
#' ## Function to read in recall plot
#' 
# correction <- methods[3]
readPrecisionRecallData <- function(correction, psiType, dataset, cores=4){
    
    message(date(), " ", correction, ", ", psiType, ":")
    fds <- loadFraseRDataSet(file.path(baseDir, psiType, correction), dataset)
    
    dtAll <- createPrecRecTable(fds=fds, psiType=psiType, correction=correction,
            outliers="byJunctionGroup", nMeanBins=1, dPsiBins=1, q="best",
            BPPARAM=MulticoreParam(cores))
        
    ans <- readBootstrapData(dtAll[[1]], rankType=c("zscore", "pvalue"), 
            maxRows=1e6, dPsiCutoff=0, FDR_LIMIT=NA, zScore_LIMIT=NA, 
            mc.cores=cores)
    
    return(ans)
}

#' ## Overall precision-recall
#+ get precision recall
dtls <- bplapply(methods, psiType=psiType, dataset=dataset, 
        FUN=readPrecisionRecallData, BPPARAM=MulticoreParam(5))
dt2p <- rbindlist(dtls)

#+ plot overall precision-recall
digits <- 3
p2plot <- !duplicated(dt2p[,.(a=round(precision, digits), 
        b=round(recall, digits), correction, type)])

dt <- dt2p[p2plot]
dt <- dt[,method:=paste(type, correction)]
for(i in unique(dt$method)){
    tmp <- dt[!is.na(cutoff) & method == i & recall > 0.8 & precision < 0.04][min(rank) == rank]
    getPred <- function(re){
        d <- 1-tmp$recall
        m <- -tmp$precision/d
        tmp$precision + m*(re - tmp$recall)
    }
    dt[method == i & recall > tmp$recall, precision:=getPred(recall)]
    tmp[,cutoff:=NA]
    tmp[,c("recall", "precision", "lower", "upper"):=list(1,0,0,0)]
    dt <- rbind(dt, tmp)
}
cutoffdt <- rbindlist(list(
    cbind(cut=0.5, dt[abs(cutoff) < 0.5][grepl("P-value", type)][,
            .SD[max(rank) == rank],by="type,correction"]),
    cbind(cut=0.1, dt[abs(cutoff) < 0.1][grepl("P-value", type)][,
            .SD[max(rank) == rank],by="type,correction"]),
    cbind(cut=0.05, dt[abs(cutoff) < 0.05][grepl("P-value", type)][,
            .SD[max(rank) == rank],by="type,correction"]),
    cbind(cut=2, dt[abs(cutoff) > 2][grepl("zScore", type)][,
            .SD[max(rank) == rank],by="type,correction"]),
    cbind(cut=3, dt[abs(cutoff) > 3][grepl("zScore", type)][,
            .SD[max(rank) == rank],by="type,correction"])))

#+ correct names
rename_methods <- function(dt, fds_imple){
    levels(dt$inj_value)<- gsub(
        "\\[0\\.000[0-9]+,", "[0.2,", levels(dt[,inj_value]))
    dt[correction == "BB", correction:="naïve BB"]
    dt[,correction:=factor(correction)]
    levels(dt$correction)[levels(dt[,correction]) == fds_imple] <- "FRASER"
    levels(dt$correction) <- gsub("\\+", "-", levels(dt[,correction]))
    dt[,type:=factor(type)]
    levels(dt$type) <- gsub("zScore rank", "\n+ z score", levels(dt$type))
    levels(dt$type) <- gsub("P-value rank", "\n+ P value", levels(dt$type))
    if(fds_imple == "PCA"){
        dt[correction == "FRASER" & type == "\n+ z score", correction:="PCA"]
    }
    if("PCA-BB-Decoder" %in% levels(dt$correction)){
        levels(dt$correction)[levels(dt$correction) == "PCA-BB-Decoder"] <- 
            "PCA + BB-Decoder"
    }
    dt
}

data <- rename_methods(dt, "PCA")
cutoffs <- rename_methods(cutoffdt, "PCA")
q <- as.character(unique(data$q))

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

data <- subset_methods(data[,-"cutoff"])
cutoffs[,cutoff:=cut]
cutoffs <- subset_methods(cutoffs)
cutoffs[,Cutoff:=factor(cutoff)]

#' 
#' plot precision recall plot
#' 
ggAll <- ggplot(data, aes(recall, precision, color=Method)) + 
    geom_vline(xintercept=c(0.8, 0.85, 0.9, 0.95), col="gray70") + 
    geom_line() + 
    geom_point(data=cutoffs, aes(x=recall, y=precision, shape=Cutoff), size=3) + 
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    theme_bw() +
    ggtitle(paste0(dataset, " outlier recall (", psiType, ", q=", q, ", TP=", max(data[,TP]), ")")) + 
    scale_color_brewer(palette='Dark2') +
    scale_fill_brewer(palette='Dark2') + 
    theme(legend.position = "bottom", legend.box = "vertical") + 
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Method), alpha=0.2, color=NA)
ggAll

#'
#' Some statistics for the paper
#' 
statsdt <- rbindlist(lapply(unique(data$Method), function(x){
    rbindlist(lapply(c(0.8, 0.85, 0.9), function(y){
        data[Method == x & recall > y][,.(Method, rank, TP,
                precision=round(precision, 2), recall=round(recall, 2), cut=y)][1]
    }))}))
statsdt[,fc:=round(precision/precision[grepl("FRASER", Method)], 2),by=cut]
statsdt[order(cut, Method)][grepl("(PCA|FRASER)\\n", Method, perl=TRUE)]
cutoffs[,.(Method, rank, TP, precision=round(precision, 2), recall=round(recall, 2), Cutoff)][grepl("(PCA|FRASER)\\n", Method, perl=TRUE)]

#+ save plot as PDF and PNG
ggsave(outPng, ggAll, 'png', width=10, height=8)
ggsave(outPdf, ggAll, 'pdf', width=10, height=8)
saveRDS(list(dt2p=data, cutoffdt=cutoffs, ggplot=ggAll), file=ggplotFile)
