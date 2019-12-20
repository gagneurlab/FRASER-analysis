#'---
#' title: "FraseR Benchmark"
#' author: Christian Mertes
#' output:
#'     html_document:
#'         pandoc_args: ["+RTS", "-K64m","-RTS"]
#'         css: ../../js-css-html/fraser-markdown-report-style.css
#'         toc: true
#'         toc_float: true
#'
#'---

#+ load knitr environment, echo=FALSE, cache=FALSE
library(knitr)
library(rmarkdown)
knitrwidth <- 200
ren <- function(){
    try(system(paste(". ~/.bashrc; module load i12g/anaconda;",
            "snakemake --keep-target-files",
            file.path(
                "/s/project/fraser/analysis/datasets/reports",
                "kremer-bader-et-al/benchmark-FraseR.html"))))
    snakemakeLIST <<- readRDS(file.path("tmp/snakemake/src/r/benchmark/",
            "fraser-benchmark.R/datasets:kremer-bader-et-al.RDS"))
    snakemake <<- snakemakeLIST$snakemake
}

#
#+ load config other needed resources for the anaylsis, echo=FALSE, cache=FALSE
suppressPackageStartupMessages({
    source("./src/r/config.R")
    #library(vioplotx)
    library(pROC)
    library(devtools)
    library(parallel)
})

#+ load mito data
load_mitomap()
load_mito_paper_data()

#+ load FraseR package, echo=FALSE, cache=FALSE
load_all(PKG_ROOT)
opts_chunk$set(echo=FALSE, fig.width=12, fig.height=8, cache=TRUE)

#'
#' # Load data
#'
benchmark_results <- list()
rocOptions <- data.table(
    algo         = c(   "FraseR", "leafcutter", "cummings-et-al"),
    displayName  = c(   "FraseR", "Leafcutter", "Cummings et al"),
    psiChangeCut = c(       0.01,           NA,              0.3),
    color        = c("firebrick",  "darkgreen",       "darkblue"),
    aucY         = c(       0.28,         0.18,             0.08)
)

#'
#' * Load benchmark set
#'
benchmark_set <- fread(snakemake@params$benchmarkSet, fill=TRUE)[,1:8,with=FALSE]
benchmark_set <- benchmark_set[hgnc_symbol != "LTBP1"] # is not a true causing gene
benchmark_set

hist(benchmark_set[,.N,by="RNA_ID,hgnc_symbol"][,.N,by="RNA_ID"][,N],
     main="Histogram of #genes with defect per sample",
     xlab="Genes with defect per sample")

#'
#' * Load FraseR data object
#'
#+ load data
name       <- snakemake@params$name
workingDir <- snakemake@params$workingDir
fds        <- loadFraseRDataSet(workingDir, name)
fds

#'
#' extract results
#'
fdsres <- results(fds, sampleIDs=unique(benchmark_set$RNA_ID),
        fdrCut=1, zscoreCut=0)
fdsres <- swapPsiSiteValue(fdsres)
fdsres <- fdsres[!is.na(fdsres$hgnc_symbol)]
res <- fdsres[abs(fdsres$zscore) >= 2.5 & abs(fdsres$deltaPsi) > 0.2]
benchmark_results['FraseR'] <- list(as.data.table(res)[,.(
    RNA_ID=sampleID, hgnc_symbol, pvalue, psiChange=abs(deltaPsi), zscore=zscore)])

#'
#' load cummings results
#'first-kremer-bader-et-al-filteredJunctions.txt
res <- fread("/s/project/fraser/analysis/datasets/cummings-et-al/kremer-bader-et-al-filteredJunctions.txt")
res[,psiChange:=as.numeric(gsub("\\*", "", normValue))]
res[,RNA_ID:=gsub(".bam", "", bamFile)]
benchmark_results['cummings-et-al'] <- list(res[,.(
    RNA_ID, hgnc_symbol, pvalue=1/psiChange, psiChange=psiChange)])

#'
#' load leafcutter results
#'
res <- readRDS("/s/project/fraser/analysis/datasets/leafcutter/kremer-bader-et-al-results.RDS")
res <- res[!is.na(hgncid) & RNA_IDs %in% benchmark_set[,RNA_ID]]
res[is.na(p), p:=1]
res[is.na(loglr), loglr:=min(res[,loglr], na.rm=TRUE)]
benchmark_results['leafcutter'] <- list(res[,.(
    RNA_ID=RNA_IDs, hgnc_symbol=hgncid, pvalue=p, psiChange=loglr)])


algo2Benchmark <- c("FraseR", "leafcutter", "cummings-et-al")
names(algo2Benchmark) <- algo2Benchmark
rankVsHitsObj <- lapply(algo2Benchmark, getRocObj,
        benchmark_results, "pvalue", rocOptions, rankVsHits)
rocPValObjs <- lapply(algo2Benchmark, getRocObj,
        benchmark_results, "pvalue", rocOptions, roc)
rocPsiValObjs <- mclapply(algo2Benchmark, getRocObj,
        benchmark_results, "psiChange", rocOptions, roc)
rocSPValObjs <- lapply(algo2Benchmark, getRocObj,
        benchmark_results, "pvalue", rocOptions, simple_roc)
rocSPsiValObjs <- lapply(algo2Benchmark, getRocObj,
        benchmark_results, "psiChange", rocOptions, simple_roc)
rocSRPValObjs <- lapply(algo2Benchmark, getRocObj,
                       benchmark_results, "pvalue", rocOptions, simple_recall)
rocSRPsiValObjs <- lapply(algo2Benchmark, getRocObj,
                         benchmark_results, "psiChange", rocOptions, simple_recall)
rocPVal100Objs <- mclapply(algo2Benchmark, getRocObj, direction=">", percent = FALSE, algorithm=3,
                        benchmark_results, "pvalue", rocOptions, roc, topN=1000)
rocPVal8000Objs <- rocPVal1000Objs
sapply(algo2Benchmark, plotFraseRRoc, rocPVal100Objs, rocOptions)


sapply(algo2Benchmark, plotFraseRRoc, rocSPValObjs, rocOptions)
sapply(algo2Benchmark, plotFraseRRoc, rocPValObjs, rocOptions)
sapply(algo2Benchmark, plotFraseRRoc, rocPsiValObjs, rocOptions)
sapply(algo2Benchmark, plotFraseRRoc, rocSPsiValObjs, rocOptions)
sapply(algo2Benchmark, plotFraseRRoc, rocPVal1000Objs, rocOptions)


#'
#' * get results for the given samples
#'   * FDR cutoff:    <= 1.0
#'   * delta PSI cutoff:    <= 0.3
#'
unique(benchmarkSet[,RNA_ID])
fdsres   <- results(fds, sampleIDs=unique(benchmarkSet$RNA_ID), fdrCut=1)
fdsres   <- swapPsiSiteValue(fdsres)

numFiltered <- list("unfiltered"=length(fdsres))
fdsres <- fdsres[!is.na(fdsres$hgnc_symbol)]
numFiltered["withGeneName"] <- length(fdsres)
fdsres <- fdsres[fdsres$psiValue    <= 0.6]
numFiltered["afterPSIvalue"] <- length(fdsres)
fdsres <- fdsres[order(fdsres$pvalue)]

#'
#' * Number of filtered results before the ROC curve
numFiltered

#+ prepare benchmark all roc
sgpdt <- as.data.table(fdsres)[,.(pvalue=min(pvalue)), by="sampleID,hgnc_symbol"]
#sgpdt <- merge(sgpdt, as.data.table(fdsres), all.x=TRUE, by=c("sampleID", "hgnc_symbol", "pvalue"))
bsgdt <- unique(benchmarkSet[,.(sampleID=RNA_ID, hgnc_symbol, ROC=TRUE)])
rocdt <- merge(sgpdt, bsgdt, all=TRUE, by=c("sampleID", "hgnc_symbol"))[order(pvalue)][!is.na(hgnc_symbol)]
rocdt[is.na(ROC), ROC:=FALSE]
#rocdt[psiValue <= 0.6]
rocdt[is.na(pvalue), pvalue:=2]

#'
#' # ROC curve
#'
#' * All samples with all data with paper benchmark set
#'
roc2plot <- unique(rocdt[,.(pvalue=min(pvalue), ROC),by=c("sampleID", "hgnc_symbol")])
rocObj <- roc(roc2plot[,ROC], roc2plot[,pvalue])
plot(rocObj)
text(0.5,0.1, paste("AUC: ", round(rocObj$auc, 3)))

#'
#' * additional zscore cutoff: abs(z) >= 3.5
#'
roc2plot <- unique(rocdt[abs(zscore) >=3.5 | ROC,.(pvalue=min(pvalue), ROC),
        by=c("sampleID", "hgnc_symbol")])
rocObj <- roc(roc2plot[,ROC], roc2plot[,pvalue], auc = TRUE)
plot.roc(rocObj, col="firebrick", grid=TRUE, plot.auc=TRUE)
text(0.5,0.1, paste("AUC: ", round(rocObj$auc, 3)))


#'
#' ## ROC curve against splice sites in WES data
#'
ssBench <- mclapply(benchmarkSet$RNA_ID, mc.cores=10, function(id){
    vcfID <- MITOMAP_DATATABLE[RNA_ID == id, EXOME_ID]
    vcfFile <- file.path("/s/project/metaIBD/annotated_r_object_tables/",
            "vep_anno_sample_objects", paste0(vcfID, "-1-annotation_vep.vcf.gz"))
    if(!file.exists(vcfFile)){
        warnings("File does not exist", vcfFile, vcfID, id)
        return(data.table())
    }
    vcf <- readVcf(vcfFile)

    # filter it
    vcfFilt <- filterForSNP(vcf)
    vcfFilt <- filterQC(vcfFilt, minAllel=1)
    vcfFilt <- filterCSQbyRegex(vcfFilt, 'splice_[ad]')

    data.table(sampleID=id, hgnc_symbol=unlist(info(vcfFilt)$DSGENE))
})
ssBench <- unique(rbindlist(ssBench))
ssBench[,ROC:=TRUE]
rocdt <- merge(sgpdt, ssBench, all.x=TRUE, by=c("sampleID", "hgnc_symbol"))[order(pvalue)][!is.na(hgnc_symbol)]
rocdt[is.na(ROC), ROC:=FALSE]
rocdt <- rocdt[psiValue <= 0.6]
rocdt
hist(ssBench[,.N,by="sampleID,hgnc_symbol"][,.N,by="sampleID"][,N],
     main="Histogram of genes with at least a var per sample",
     xlab="Genes with at least one variant")



#'
#' * only psiValue cutoff <= 0.6
#'
roc2plot <- unique(rocdt[,.(pvalue=min(pvalue), ROC),by=c("sampleID", "hgnc_symbol")])
rocObj <- roc(roc2plot[,ROC], roc2plot[,pvalue])
plot(rocObj)
text(0.5,0.1, paste("AUC: ", round(rocObj$auc, 3)))
hist(roc2plot[ROC==TRUE,.N,by=sampleID][,N],
     main="Histogram of genes with at least a var per sample",
     xlab="Genes with at least one variant")

#'
#' * additional zscore cutoff: abs(z) >= 3.5
#'
roc2plot <- unique(rocdt[abs(zscore) >=3.5 | ROC,.(pvalue=min(pvalue), ROC),
        by=c("sampleID", "hgnc_symbol")])
rocObj <- roc(roc2plot[,ROC], roc2plot[,pvalue])
plot(rocObj)
text(0.5,0.1, paste("AUC: ", round(rocObj$auc, 3)))



par(cex=1.4)
plot(NA, xlim=c(1,180000), ylim=c(0,nrow(benchmark_set)), xlab="Rank", ylab="cumSum(TRUE hits)") #log="x",
sapply(algo2Benchmark, rankVsHitsObj, rocOptions, FUN=function(name, data, opt){
    da <- data[[name]]
    currOpt <- opt[algo==name]
    l <- length(da)
    maxL <- min(l, 180000)
    lines(1:maxL, da[1:maxL], col=currOpt[,color])
})
grid()
legend("bottomright", rocOptions[,displayName], col=rocOptions[,color], pch=20, lty=1)

plot(NA, xlim=c(1,100), ylim=c(0,nrow(benchmark_set)), xlab="Rank", ylab="cumSum(TRUE hits)") #log="x",
sapply(algo2Benchmark, rankVsHitsObj, rocOptions, FUN=function(name, data, opt){
    da <- data[[name]]
    currOpt <- opt[algo==name]
    l <- length(da)
    maxL <- min(l, 100)
    lines(1:maxL, da[1:maxL], col=currOpt[,color])
})
grid()
legend("topleft", rocOptions[,displayName], col=rocOptions[,color], pch=20, lty=1)





