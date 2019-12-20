#'---
#' title: "QQ-plot junction analysis"
#' author: Christian Mertes
#' output:
#'     html_document:
#'         pandoc_args: ["+RTS", "-K64m","-RTS"]
#'         css: ../../js-css-html/fraser-markdown-report-style.css
#'         toc: true
#'         toc_float: true
#' #params:
#' #  workingDir: /tmp/FraseR
#' #  name: test
#' #  threads: 15
#' #  resultDir: /tmp/webspace/fraser/
#' #
#'---
CURR_REPORT_FILE <<- "src/r/benchmark/qq-plot-analysis.R"
if("params" %in% ls()){
    attach(readRDS(file.path("./tmp/snakemake/src/r/analysis/runFraseR.R/",
            "datasets:kremer-bader-et-al.RDS")))
}

#
#+ load FraseR and other needed resources for the anaylsis, cache=TRUE
suppressMessages(source("./src/r/config.R"))
load_all(PKG_ROOT)

#+ load knitr environment1, echo=FALSE, cache=FALSE
opts_chunk$set(echo=FALSE, cache=TRUE)
set.seed(42)

#
#+ functions
qqfitlm <- function(idx, obs, exp, plot.it=FALSE){
    observed_p <- obs[idx,]
    expected_p <- exp[idx,]

    fit       <- lm(observed_p ~ expected_p)
    slope     <- fit$coefficients['expected_p']
    intercept <- fit$coefficients[1]

    if(plot.it){
        lim <- c(0,max(c(observed_p, expected_p)))
        plot(expected_p, observed_p, xlim=lim, ylim=lim)
        abline(0,1, col="gray", lty="dotted")
        abline(intercept,slope, col="firebrick")
    }

    return(slope)
}



#'
#' # Load data
#'
#' The provided RNA_IDs are used
#'
#' Filter data
#+ filter samples, echo=T
fds <- loadFraseRDataSet(params$workingDir, params$name)

#' get pvalues
#+get pvalues
psiType <- "psi3"
se <- as(fds, "RangedSummarizedExperiment")
if(psiType == "psiSite") se <- nonSplicedReads(fds)
tested <- na2false(mcols(fds, type=psiType)[paste0(psiType, "_tested")])
pvalues <- as(assays(fds)[[paste0("pvalue_", psiType)]][tested,], "matrix")

#' remove full NA rows
#+remove nas from pvalues
rangeNames <- paste0(
        seqnames(se[tested]), ":", start(se[tested]), "-", end(se[tested]))
rownames(pvalues) <- rangeNames
pvalues <- pvalues[apply(pvalues, 1, function(x) sum(is.na(x))/length(x) < 0.5),]

#+ plot 1000 random qq plots
sampledRows <- sort(sample(1:dim(pvalues)[1], 1000))
qqp <- fraserQQplotPlotly(pvalues[sampledRows,], sampleWise = FALSE)
qqp

#+ fit lm to qqplots
obs <- do.call(rbind, sapply(qqp$x$visdat, function(x) x()$observ))
exp <- do.call(rbind, sapply(qqp$x$visdat, function(x) x()$expect))
qqfits <- sapply(1:dim(obs)[1], qqfitlm, obs, exp)

#+ plot most outliers
highqqidx <- which(rank(qqfits) > length(qqfits) - 9)
lowqqidx  <- which(rank(qqfits) <= 9)
par(mfrow=c(3,3))
sapply(highqqidx, qqfitlm, obs, exp, TRUE)
sapply(lowqqidx, qqfitlm, obs, exp, TRUE)

#+ idx1high
idx <- 2
high <- TRUE

#+ plot counts
curidx <- ifelse(high, highqqidx[idx], lowqqidx[idx])
seqname  <- strsplit(qqp$x$attrs[[curidx+2]]$name, ":|-")[[1]][1]
startpos <- strsplit(qqp$x$attrs[[curidx+2]]$name, ":|-")[[1]][2]
endpos   <- strsplit(qqp$x$attrs[[curidx+2]]$name, ":|-")[[1]][3]
curgr <- granges(se[seqnames(se) == seqname & start(se) == startpos & end(se) == endpos])
curgr$type = "psi3"
plotJunctionDistribution(fds, curgr, valueVsCounts = TRUE, qqplot = TRUE)
fds[seqnames(fds)]



cfds <- fds[na2false(mcols(fds, type="psi3")$hgnc_symbol == "TIMMDC1")][104]
as.matrix(assays(cfds)[["pvalue_psi3"]][,"MUC1344"])
curgr <- granges(cfds)
curgr$type <- "psi3"
plotJunctionDistribution(fds, curgr, valueVsCounts = TRUE, qqplot = TRUE)
plotCountsAtSite(curgr, fds, curgr$type, plotLog=TRUE)
debug(plotCountsAtSite)
undebug(getPlotDistributionData)
par(mfrow=c(1,1))

