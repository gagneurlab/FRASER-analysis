#'---
#' title: Analyse FraseR's betaBinomial-test runtime
#' author: Christian Mertes
#' output:
#'   html_document:
#'     pandoc_args: ["+RTS", "-K64m","-RTS"]
#'     toc: yes
#'     css: ../../../../gagneurlab_shared/r/knitr_plugins/add_content_table.css
#'---
#+ load packages, echo=FALSE, cache=FALSE
if(FALSE){
    # render it
    library(knitr)
    library(rmarkdown)
    project_webserver_dir <- "/s/public_webshare/project/fraser"
    rscript <- "src/r/runtimeAnalysis/betaBinomRunTime.R"
    outdir <- file.path(project_webserver_dir, dirname(rscript))
    render(rscript, output_format = 'html_document', output_dir = outdir)
}
library(FraseR)
library(LSD)
library(plotly)
library(parallel)
options(width=140)
mcCores <- 20
opts_chunk$set(echo=FALSE, fig.width=12, fig.height=8, cache=TRUE)
#opts_chunk$set(cache=FALSE)

#'
#' # Load the data
#+ load, echo=TRUE
fds <- loadFraseRDataSet("/s/project/fraser/analysis/006-gtex", "gtexfull")
psiType  <- "psi3"
sampledPoints <- 10000
longTimeCutoff <- 55
readType <- FraseR:::whichReadType(fds, psiType)

#'
#+ get raw data
timings <- mcols(fds, type=readType)[[paste0(psiType, "_timing")]]
round0Timings <- round(timings)
rawIC <- counts(fds, type=psiType, side="ofInterest")
rawOC <- counts(fds, type=psiType, side="otherSide")
rawICdt <- as.data.table(rawIC)
rawOCdt <- as.data.table(rawOC)

#'
#+ needed calculated data, cache=FALSE
tested  <- FraseR:::na2false(unlist(
    mcols(fds, type=psiType)[paste0(psiType, "_tested")]
))
testedTime <- timings[tested]
longtimes  <- sapply(timings > longTimeCutoff, isTRUE)
shorttimes <- sapply(timings <  0.5, isTRUE)
sampled <- unique(c(which(longtimes), sample(which(tested), sampledPoints)))

#'
#+ overview, echo=FALSE
#' # Data overview
#' * Number of tested positions
table(tested)
#' * Number of long running tests
table(longtimes)
#' * Histogram of timings
hist(log10(timings[tested]), main="Histogram of test timings", breaks = 50)
#' * Raw counts for long running tests
#'   * rawCounts
rawIC[longtimes,]
#'   * rawOtherCounts
rawOC[longtimes,]

#'
#+ heatscatter
#' # Detailed informations
#' ## Heatscatter of the first 10 long running tests
#' * raw counts versus raw other counts
#+ heatscatter longtime
devNULL <- sapply(which(longtimes)[1:10], function(i){
    rc1 <- as.vector(rawIC[i,])
    rc2 <- as.vector(rawOC[i,])
    cor <- cor(rc1, rc2)
    heatscatter(log10(1+rc1), log10(1+rc2),
        xlab="log10(1+rawCounts)", ylab="log10(1+rawOtherCounts",
        main=paste0("id: ", i,
                "\ttiming: ", round(timings[i],2),
                "\tcor: ", round(cor,2)
        )
    )
})

#'
#' ## Heatscatter of the time versus total counts
FUN <- function(idx, x, f) f(unlist(x[idx]))
medianIC <- unlist(mclapply(sampled, f=median, x=rawICdt, mc.cores=mcCores, FUN=FUN))
medianOC <- unlist(mclapply(sampled, f=median, x=rawOCdt, mc.cores=mcCores, FUN=FUN))
sdIC <- unlist(mclapply(sampled, x=rawICdt, f=sd, mc.cores=mcCores, FUN=FUN))
sdOC <- unlist(mclapply(sampled, x=rawOCdt, f=sd, mc.cores=mcCores, FUN=FUN))

#' * Median of total counts versus time
#+ heat time counts median
heatscatter(timings[sampled], main="Time versus total counts",
    log10(1 + medianIC + medianOC),
    xlab="time (sec)", ylab="log10(1+median(total raw counts))"
)

#+ heat time counts sd
#' * Standard deviation of total counts versus time
heatscatter(timings[sampled], main="Time versus total counts",
    log10(1 + sdIC + sdOC),
    xlab="time (sec)", ylab="log10(1+sd(total raw counts))"
)

#+ summed timing surraounding
#' ## Counts and running time are dependend on the posision
#'
run  <- 50
surrTime <- unlist(mclapply(1:sum(tested), mc.cores=mcCores,
        FUN=function(i,r=run) sum(testedTime[max(0,i-r):min(sum(tested), i+r)])
))
stt <- sample(testedTime)
randSurrTime <- unlist(mclapply(1:sum(tested), mc.cores=mcCores,
    FUN=function(i,r=run) sum(stt[max(0,i-r):min(sum(tested), i+r)])
))
stdf <- data.frame(surrTime=surrTime, time=testedTime, pos=1:length(surrTime),
        randSurrTime=randSurrTime, randTime=stt
)
p <- plot_ly(stdf, type="scattergl", mode="lines")
p <- add_trace(p, x=~pos, y=~surrTime, mode="lines", name="summedTime")
p <- add_trace(p, x=~pos, y=~time, mode="markers", name="time")
p <- add_trace(p, x=~pos, y=~randSurrTime, mode="lines", name="randSummedTime")
p <- add_trace(p, x=~pos, y=~randTime, mode="markers", name="randTime")
p

# /*
#scatterplot3d(timings[sampled],
#        log10(1+rowMedians(as.matrix(rawIC[sampled,]))),
#        log10(1+rowMedians(as.matrix(rawOC[sampled,])))
#)
#
#assays(longfds)$rawOtherCounts_psi3)
#
#dtrrc <- data.table(
#    y=as.vector(assays(sfds)$rawCountsJ[which(longtimes)[1],]),
#    N=as.vector(assays(sfds)$rawOtherCounts_psi3[which(longtimes)[1],])
#)
#testFun <- function(n, mu1, sd1, mu2, sd2){
#    dt <- data.table(
#        y=abs(as.integer(rnorm(n, mu1, sd1))),
#        N=abs(as.integer(rnorm(n, sd2, sd2)))
#    )
#    dt[1,c("y", "N"):=list(y=0,N=0)]
#    betabinVglmTest(countMatrix=as.matrix(dt), y=dt[,y], N=dt[,N])
#}
#mb <- microbenchmark::microbenchmark(
#    times = 100,
#    a = testFun(200, 0, 5,  100,   10),
#    b = testFun(200, 0, 5, 1000,   10),
#    c = testFun(200, 0, 5,  100,  100),
#    d = testFun(200, 0, 5, 1000, 1000),
#    e = testFun(200, 100, 100, 1000, 1000),
#    f = testFun(200, 1000, 1000, 10000, 1000),
#    g = testFun(200, 100, 100, 10000, 10000),
#    h = betabinVglmTest(countMatrix=as.matrix(dtrrc), y=dtrrc[,y], N=dtrrc[,N])
#)
#ggplot2::autoplot(mb)
#
# */

