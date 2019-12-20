#'---
#' title: Joined benchmark FraseR
#' author: Christian Mertes
#' wb:
#'  input:
#'   - res: '`sm expand("Data/benchmark/insilicoInjection/efreq{{efreq}}_n{{nsamples}}_dp{dp}_{{dataset}}_res_cum_true_hits.rds", dp=["0.3", "0.4", "0.5"])`'
#'  output:
#'   - wBhtml: 'Output/html/Benchmark/inSilicoBenchmark/efreq{efreq}_n{nsamples}_{dataset}_final_benchmark.html'
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")

#'
#' # Dataset
#+ echo=TRUE
set.seed(42)
snakemake@wildcards$dataset
benchFiles <- snakemake@input$res
benchls <- lapply(benchFiles, fread)

resls <- lapply(benchls, function(x) x[,FraseR])
names(resls) <- gsub("_.*", "", gsub(".*_dp", "dp", benchFiles))
nhits <- max(sapply(resls, max))

par(cex=1.4)
col <- c("firebrick", "forestgreen", "darkblue")
for(n in c(500, 10000)){
    plot(NA, xlim=c(1,n), ylim=c(0,1), xlab="Rank",
         ylab=paste0("cumSum(TRUE hits)/hits (n=", nhits, ")"),
         main=paste("Benchmark FraseR: \ndp:",
                ", ns=", snakemake@wildcards$nsamples,
                ", efreq=", snakemake@wildcards$efreq))
        for(i in 1:length(resls)){
            da <- resls[[i]]
            l <- length(da)
            maxL <- min(l, n)
            lines(1:n, da[1:n]/nhits, col=col[i])
        }
        abline(0,1/nhits, lty=3)
        grid()
        legend("topleft", ncol=2, names(resls), col=col, pch=20, lty=1)
}
