#'---
#' title: Validate benchmark
#' author: Christian Mertes
#' wb:
#'  input:
#'   - injections: 'Data/benchmark/insilicoInjection/efreq{efreq}_n{nsamples}_dp{dp}_{dataset}_event_table.tsv'
#'   - pvals: '`sm config["DATADIR"] + "/savedObjects/inSilico_efreq{efreq}_n{nsamples}_dp{dp}_{dataset}/pvalue_psiSite.h5"`'
#'  output:
#'   - res: 'Data/benchmark/insilicoInjection/efreq{efreq}_n{nsamples}_dp{dp}_{dataset}_res_cum_true_hits.rds'
#'   - wBhtml: 'Output/html/Benchmark/inSilicoBenchmark/efreq{efreq}_n{nsamples}_dp{dp}_{dataset}_benchmark.html'
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")

#'
#' # Dataset
#+ echo=TRUE
set.seed(42)
snakemake@wildcards$dataset
inpname <- basename(dirname(snakemake@input$pvals))
oriname <- gsub("_efreq[^_]+|_dp[^_]+", "", inpname)
wddir   <- dirname(dirname(dirname(snakemake@input$pvals)))
oriFds  <- loadFraseRDataSet(wddir, oriname)
fds     <- loadFraseRDataSet(wddir, inpname)
benchset <- fread(snakemake@input$injections)
benchsetgr <- makeGRangesFromDataFrame(benchset, keep.extra.columns = TRUE)

dim(fds)
dim(oriFds)

plotJunctionDistribution(oriFds,  benchsetgr[1], plotLegend=FALSE, plotValVsCounts=character(0), plotRank=character(0))
plotJunctionDistribution(fds,     benchsetgr[1], plotLegend=FALSE, plotValVsCounts=character(0), plotRank=character(0))

plotJunctionDistribution(oriFds,  benchsetgr[2], plotLegend=FALSE, plotValVsCounts=character(0), plotRank=character(0))
plotJunctionDistribution(fds,     benchsetgr[2], plotLegend=FALSE, plotValVsCounts=character(0), plotRank=character(0))

plotJunctionDistribution(oriFds,  benchsetgr[3], plotLegend=FALSE, plotValVsCounts=character(0), plotRank=character(0))
plotJunctionDistribution(fds,     benchsetgr[3], plotLegend=FALSE, plotValVsCounts=character(0), plotRank=character(0))

res <- results(fds, fdrCut=1, dPsiCut=0, zscoreCut=0)
res2 <- results(oriFds, fdrCut=1, dPsiCut=0, zscoreCut=0)
fds <- saveFraseRDataSet(fds)
oriFds <- saveFraseRDataSet(oriFds)

#'
#' prepare data
#'
rocOptions <- data.table(
    algo         = c(   "FraseR", "leafcutter", "cummings-et-al"),
    displayName  = c(   "FraseR", "Leafcutter", "Cummings et al"),
    psiChangeCut = c(        0.0,           NA,              0.3),
    color        = c("firebrick",  "darkgreen",       "darkblue"),
    aucY         = c(       0.28,         0.18,             0.08)
)
benchmark_set <- benchset
benchmark_set[,RNA_ID:=sampleID]
benchmark_results <- list(FraseR=as.data.table(res)[,.(
    RNA_ID=sampleID, hgnc_symbol, pvalue, psiChange=abs(deltaPsi), zscore=zscore)])
rankVsHitsObj <- lapply(c(FraseR="FraseR"), getRocObj, data=benchmark_results,
        predictionName="pvalue", rocOptions=rocOptions, rankVsHits, allSamples=TRUE)


#'
#' # ROC-like curve of injections
#'

par(cex=1.4)
for(n in c(500, 10000)){
plot(NA, xlim=c(1,n), ylim=c(0,1), xlab="Rank",
        ylab=paste0("cumSum(TRUE hits)/hits (n=", nrow(benchmark_set), ")"),
        main=paste("Benchmark FraseR: \ndp:", snakemake@wildcards$dp,
                ", ns=", snakemake@wildcards$nsamples,
                ", efreq=", snakemake@wildcards$efreq))
    name <- "FraseR"
    da <- rankVsHitsObj[[name]]
    l <- length(da)
    maxL <- min(l, n)
    lines(1:n, da[1:n]/nrow(benchmark_set), col="firebrick")
    abline(v=sum(res$p.adj <= 0.1), lty=2)
    abline(0,1/nrow(benchmark_set), lty=3)
    grid()
    legend("bottomright", c("FraseR", "p.adj <= 0.1"), col=c("firebrick", "black"), pch=c(20, NA), lty=c(1,2))
}

#'
#' # sample and gene wide qqplot
plotSampleQQ(fds, maxOutlier = 8)

#'
#' write res
write_tsv(rankVsHitsObj, snakemake@output$res)

if(FALSE){
    parallel(oriFds) <- MulticoreParam(50, 500, progressbar=TRUE)
    register(parallel(oriFds))

    assays(oriFds)[['pvalue_psi5']] <- NULL
    oriFds <- calculatePValues(oriFds)
    assays(oriFds)[['pvalue_psi3']] <- NULL
    oriFds <- calculatePValues(oriFds)
    assays(oriFds)[['pvalue_psiSite']] <- NULL
    oriFds <- calculatePValues(oriFds)
    oriFds <- calculatePValues(oriFds)
    oriFds <- calculatePValues(oriFds)

}
