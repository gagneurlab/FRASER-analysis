#'---
#' title: Plot injected dPSI vs fitted dPSI(=|observed/data PSI - correctedPSI|)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - inHtml: 'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_pvalues.html'
#'   - trueOut: '`sm config["DATADIR"] + "/datasets/inject_{delta}/{psiType}/savedObjects/{dataset}/trueOutliers_{psiType}.h5"`'
#'  output:
#'   - figurePng: '`sm config["FIGDIR"] + "/injectedVsFittedDpsi/{dataset}/{delta}/{psiType}/{method}.png"`'
#'   - ggplot:    '`sm config["DATADIR"] + "/processedData/figData/injectedVsFittedDpsi/{dataset}/{delta}/{psiType}/{method}_ggplot.RDS"`'
#'   - wBhtml:  'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_injectedDpsi.html'
#'  type: noindex
#'---

if(FALSE){
  snakemake <- readRDS("./tmp/snakemake.RDS")

  fds_bb <- loadFraseRDataSet("Data/paperPipeline/datasets/inject_delta/psi3", "SimulationBB", TRUE)
  fds_dm <- loadFraseRDataSet("Data/paperPipeline/datasets/inject_delta/psi3", "SimulationDM", TRUE)

}

#+ source main config
source("./src/r/config.R")

#+ input
fdsFile    <- snakemake@input$trueOut
dataset    <- snakemake@wildcards$dataset
psiType    <- snakemake@wildcards$psiType
fitMethod  <- snakemake@wildcards$method
workingDir <- file.path(dirname(dirname(dirname(fdsFile))), fitMethod)

#+ output
ggplotFile <- snakemake@output$ggplot
outPng <- snakemake@output$figurePng

#' ## Load dataset
#+ load fds
dataset
workingDir
fds     <- loadFraseRDataSet(workingDir, dataset)

#' ## Fit parameters
#' get fit params
fitMethod
psiType

#'
#' ## Plot injected dPSI vs fitted dPSI(=|observed/data PSI - correctedPSI|)
#'

#+ get injected outliers
trueOut <- getAssayMatrix(fds, "trueOutliers", psiType)
trueDpsi <- getAssayMatrix(fds, "trueDeltaPSI", psiType)
ind <- which(trueOut != 0)
#+ get corrected psi
mu <- predictedMeans(fds, psiType)
fitPsi <- mu[ind]
#+ get data psi
k <- K(fds, psiType)
n <- N(fds, psiType)
dataPsi <- (k+1)/(n+2)
#+ get deltaPSI of fit
deltaPSI <- dataPsi-mu
#+ get mean junction coverage
nMeans <- rowMeans2(n)
#+ create data.table
nInd <- matrix(nMeans, nrow=nrow(trueOut), ncol=ncol(trueOut))[ind]
nMeanBins <- 3
quantilesBreaks <- quantile(nInd, p=(1:nMeanBins/nMeanBins))
dt1 <- data.table(outlier=trueOut[ind], tDpsi=trueDpsi[ind], fDpsi=deltaPSI[ind],
                  meanN=nInd, coverageBin=cut(nInd, breaks=c(0, quantilesBreaks), include.lowest = TRUE, ordered_result = TRUE, right=TRUE))
dtc <- dt1[, round(cor(tDpsi, fDpsi), digits=3), by=coverageBin]
dtc2 <- dt1[, round(cor(tDpsi, fDpsi), digits=3), by=outlier]
dtc3 <- dt1[, round(cor(tDpsi, fDpsi), digits=3), by=c('coverageBin', 'outlier')]
# mean squared error function
mse <- function(obs, pred){
  mean((obs - pred)^2)
}
dt_mse <- dt1[, round(mse(tDpsi, fDpsi), digits=5), by=coverageBin]
dt_mse2 <- dt1[, round(mse(tDpsi, fDpsi), digits=5), by=outlier]
dt_mse3 <- dt1[, round(mse(tDpsi, fDpsi), digits=5), by=c('coverageBin', 'outlier')]

regressionLine <- lm(fDpsi~tDpsi,data=dt1)

#+ scatterplot fit of all injected outliers
par(mar=c(5,5,6,2))
heatscatter(dt1$tDpsi, dt1$fDpsi, xlab=bquote("injected" ~ Delta * Psi), ylab=bquote(Delta * Psi ~ "= observed" ~ Psi - "predicted" ~ Psi),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of all injected outliers\n"), cor=TRUE); grid(); abline(0,1)
abline(regressionLine, col="green")
text(x=-0.8,y=0.9, labels=paste0("MSE=",round(mse(dt1$tDpsi, dt1$fDpsi), digits=5)))

require(ggplot2)

g <- ggplot(data=dt1, aes(x=tDpsi, y=fDpsi)) + geom_point(alpha = 0.7) + 
    geom_smooth(method='lm') +
    labs(x=bquote("injected" ~ Delta * Psi), 
         y=bquote(Delta * Psi ~ "= observed" ~ Psi - "predicted" ~ Psi)) +
    theme_bw()
ggsave(outPng, g, 'png')
saveRDS(g, file=ggplotFile)

g1 <- ggplot(data=dt1, aes(x=tDpsi, y=fDpsi)) + geom_point(size=1) + geom_abline(slope=1, intercept=0, color="red") +
  labs(x=bquote("injected" ~ Delta * Psi), y=bquote(Delta * Psi ~ "= observed" ~ Psi - "predicted" ~ Psi),
       title=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of all injected outliers by junction coverage")) +
  facet_grid( . ~coverageBin, labeller=label_both) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))  + coord_fixed(ratio = 1) +
  geom_text(data=dtc, aes(x=-0.8, y=0.94, label=paste("cor =", V1))) + geom_text(data=dt_mse, aes(x=-0.8, y=0.81, label=paste("mse =", V1)))
g1
g2 <- ggplot(data=dt1, aes(x=tDpsi, y=fDpsi)) + geom_point(size=1) + geom_abline(slope=1, intercept=0, color="red") +
  labs(x=bquote("injected" ~ Delta * Psi), y=bquote(Delta * Psi ~ "= observed" ~ Psi - "predicted" ~ Psi),
       title=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of all injected outliers by junction coverage")) +
  facet_grid(. ~ outlier, labeller=label_both) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))  + coord_fixed(ratio = 1) +
  geom_text(data=dtc2, aes(x=-0.8, y=0.94, label=paste("cor =", V1))) + geom_text(data=dt_mse2, aes(x=-0.8, y=0.81, label=paste("mse =", V1)))
g2
g3 <- ggplot(data=dt1, aes(x=tDpsi, y=fDpsi)) + geom_point(size=1) + geom_abline(slope=1, intercept=0, color="red") +
  labs(x=bquote("injected" ~ Delta * Psi), y=bquote(Delta * Psi ~ "= observed" ~ Psi - "predicted" ~ Psi),
       title=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of all injected outliers by junction coverage")) +
  facet_grid(coverageBin ~ outlier, labeller=label_both) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))  + coord_fixed(ratio = 1) +
  geom_text(data=dtc3, aes(x=-0.75, y=0.94, label=paste("cor =", V1))) + geom_text(data=dt_mse3, aes(x=-0.75, y=0.81, label=paste("mse =", V1)))
g3

#' ## Summarize by donor/acceptor site

#+ group by donor/acceptor
trueOut_grouped  <- getTrueOutlierByGroup(fds, psiType)
trueDpsi_grouped <- getAbsMaxByGroup(fds, psiType, trueDpsi)
deltaPSI_grouped <- getAbsMaxByGroup(fds, psiType, deltaPSI)
#+ get pvalues
pvals <- pVals(fds, type=psiType, byGroup = TRUE)
padj <- padjVals(fds, psiType, byGroup = TRUE)

#+ create data table
dt  <- data.table(trueOut=c(trueOut_grouped), pval=c(pvals), padj=c(padj), trueDeltaPSI=c(trueDpsi_grouped), fitDeltaPSI=c(deltaPSI_grouped) )
nMeans_grouped <- getByGroup(fds, psiType, matrix(nMeans, nrow=nrow(n), ncol=ncol(n)))
#+ add to data table
dt[, meanCoverage:=c(nMeans_grouped)]
nMeanBins <- 5
quantilesBreaks <- quantile(nMeans, p=(1:nMeanBins/nMeanBins))
dt[, coverageBin:=cut(nMeans_grouped, breaks=c(0, quantilesBreaks), include.lowest = TRUE, ordered_result = TRUE, right=TRUE)]
#+ subset for TPs
TPs <- dt[trueOut != 0]

par(mar=c(5,5,6,2))
heatscatter(abs(TPs$trueDeltaPSI), abs(TPs$fitDeltaPSI), xlab=bquote("| injected" ~ Delta * Psi ~ "|"), ylab=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of all injected outliers\n"), cor=TRUE); grid(); abline(0,1)
# #+ plot TPs binned by mean coverage
# require(ggplot2)
# g2 <- ggplot(data=TPs, aes(x=abs(trueDeltaPSI), y=abs(fitDeltaPSI))) + geom_point() + stat_density_2d(aes(fill = ..level..), geom="polygon") +
#   scale_fill_gradient(low="blue",high="yellow") + geom_abline(slope=1, intercept=0, color="red") +
#   labs(x=bquote("| injected" ~ Delta * Psi ~ "|"), y=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
#        title=paste0(dataset, ", ", psiType, ": ", fitMethod, " fit of all injected outliers by junction coverage")) +
#   facet_grid(. ~ coverageBin, labeller=label_both) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))  + coord_fixed(ratio = 1)
# g2

#' ## Plot by pvalue low/high

#+ define p-value cut between low/high
cut <- 1e-5
#+ plot pval low/high
par(mar=c(5,5,6,2))
heatscatter(TPs[pval <= cut, abs(trueDeltaPSI)], TPs[pval <= cut, abs(fitDeltaPSI)],
            xlab=bquote("| injected" ~ Delta * Psi ~ "|"), ylab=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of injected outliers with p-value <= ", cut, "\n"), cor=TRUE); grid(); abline(0,1)
par(mar=c(5,5,6,2))
heatscatter(TPs[pval > cut, abs(trueDeltaPSI)],  TPs[pval > cut, abs(fitDeltaPSI)],
            xlab=bquote("| injected" ~ Delta * Psi ~ "|"), ylab=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of injected outliers with p-value > ", cut, "\n"), cor=TRUE); grid(); abline(0,1)

#+ define p-adjust cut between low/high
cutAdj <- 5e-2
#+ plot p-adjust low/high
par(mar=c(5,5,6,2))
heatscatter(TPs[padj <= cutAdj, abs(trueDeltaPSI)], TPs[padj <= cutAdj, abs(fitDeltaPSI)],
            xlab=bquote("| injected" ~ Delta * Psi ~ "|"), ylab=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of injected outliers with adjusted p-value <= ", cutAdj, "\n"), cor=TRUE); grid(); abline(0,1)
par(mar=c(5,5,6,2))
heatscatter(TPs[padj > cutAdj, abs(trueDeltaPSI)],  TPs[padj > cutAdj, abs(fitDeltaPSI)],
            xlab=bquote("| injected" ~ Delta * Psi ~ "|"), ylab=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of injected outliers with adjusted p-value > ", cutAdj, "\n"), cor=TRUE); grid(); abline(0,1)

#+ plot mean junction coverage
par(mar=c(5,5,6,2))
heatscatter(log10(TPs[pval <= cut, meanCoverage]), TPs[pval <= cut, abs(trueDeltaPSI)], ylim=c(0,1),
            xlab="log10(mean junction coverage)", ylab=bquote("| injected" ~ Delta * Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of injected outliers with p-value <= ", cut, "\n"), cor=TRUE); grid()
par(mar=c(5,5,6,2))
heatscatter(log10(TPs[pval > cut, meanCoverage]), TPs[pval > cut, abs(trueDeltaPSI)], ylim=c(0,1),
            xlab="log10(mean junction coverage)", ylab=bquote("| injected" ~ Delta * Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n fit of injected outliers with p-value > ", cut, "\n"), cor=TRUE); grid()

#+ plot coverage vs deltaPSI for "false positives"
FPs <- dt[trueOut == 0]
par(mar=c(5,5,6,2))
heatscatter(log10(FPs[pval <= cut, meanCoverage]), FPs[pval <= cut, abs(fitDeltaPSI)], ylim=c(0,1),
            xlab="log10(mean junction coverage)", ylab=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n potential false positives with p-value <= ", cut, "\n"), cor=TRUE); grid()
par(mar=c(5,5,6,2))
heatscatter(log10(FPs[padj <= cutAdj, meanCoverage]), FPs[padj <= cutAdj, abs(fitDeltaPSI)], ylim=c(0,1),
            xlab="log10(mean junction coverage)", ylab=bquote(Delta * Psi ~ "= | observed" ~ Psi - "predicted" ~ Psi ~ "|"),
            main=paste0(dataset, ", ", psiType, ": ", fitMethod, "\n potential false positives with adjusted p-value <= ", cutAdj, "\n"), cor=TRUE); grid()

