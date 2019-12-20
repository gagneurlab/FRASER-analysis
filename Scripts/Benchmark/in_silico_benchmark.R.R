

source("./src/r/config.R")
load_all(PKG_ROOT)
library(gplots)
library(LSD)

datadir <- "/s/project/fraser/analysis/datasets/"
if(!file.exists(datadir)){
    datadir <- "~/fraser-datasets/datasets"
}

#'
#' Get smaller dataset
#' 20k junctions and 200 samples
#'
if(FALSE){
    set.seed(42)
    njunc <- 20000
    njunc <- 20000
    nsamp <- 200
    options("FraseR-hdf5-num-chunks"=6)
    fds <- loadFraseRDataSet(datadir, "gtex-skin-all")
    fds <- fds[na2false(mcols(fds, type="j")$psi3_tested)]
    fds <- fds[
        sort(sample(1:dim(fds)[1], njunc)),
        sort(sample(1:dim(fds)[2], nsamp))]

    name(fds) <- "in-silico-bench-small"

    cts <- as(counts(fds, type="psi3") + counts(fds, type="psi3", side="other"), "matrix")
    fds <- fds[apply(cts, 1, quantile, 0.75) > 5,]
    fds <- saveFraseRDataSet(fds, rewrite=TRUE)
}

#'
#' get data points
#'
set.seed(42)
fds <- loadFraseRDataSet(datadir, "in-silico-bench-small")
idx2plot <- idx2plot <- sort(sample(1:dim(fds)[1], 20))

#'
#' get heatmap of junctions
#'
cts  <- as(counts(fds, type="psi3", side="ofInt"), "matrix")
octs <- as(counts(fds, type="psi3", side="other"), "matrix")
ltcorcts <- cor(-log10(1+cts))
heatmap.2(ltcorcts, trace="no")

#'
#' random 16 qq plots
#'
par(mfrow=c(4,4), mar=c(4,4,1,0.5))
for(i in idx2plot){
    plotQQplot(granges(fds)[i], fds, type="psi3")
}

#'
#' random 16 detail plots (delta_psi, counts, qqplot)
#'
for(i in idx2plot){
    plotJunctionDistribution(fds, granges(fds)[i], type="psi3",
            plotValVsCounts=NULL, plotRank="delta_psi3", plotLegend=FALSE)
}


cts <- counts(fds, type="psi3") + counts(fds, type="psi3", side="other")
cts <- as(cts, "matrix")
table(apply(cts, 1, quantile, 0.75) > 5)
par(mfrow=c(1,1))
hist(log10(apply(cts, 1, quantile, 0.95)))

hist(cts[1,])
quantile(cts[1,], 0.75)
i <- 3751
plotJunctionDistribution(fds, granges(fds)[i], type="psi3", plotLegend=FALSE)


hist(cts[i,])
undebug(plotJunctionDistribution)
parallel(fds) <- MulticoreParam(60, 90, progressbar=TRUE)
fds <- calculatePValues(fds, alternative="two.sided")

sort(sample(1:dim(fds)[1], 20))

# assays(fds)[['pvalue_psi3']] <- NULL
