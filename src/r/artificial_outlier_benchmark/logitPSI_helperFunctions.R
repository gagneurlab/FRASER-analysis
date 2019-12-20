
# create data table containing all information needed for ploting the recovery of logit PSI (Fig2a)
getLogitPsiTable <- function(fds, dataset, psiType, correction, nMeanBins=5){
  
  # get true logit PSI 
  logitPsi      <- as.vector(getAssayMatrix(fds, "trueLogitPSI", type=psiType))
  
  # get predicted logit PSI
  if(correction == "PEER"){
    logitPsi_fit <- as.vector(getAssayMatrix(fds, "peerLogitMu", type=psiType))
  } else if(correction == "kerasBBdAE" || correction == "kerasDAE"){
    logitPsi_fit <- as.vector(qlogis(predictedMeans(fds, type=psiType)))
  } else{
    logitPsi_fit  <- as.vector(predictY(fds, type=psiType))
  }
  
  # bin by mean junction coverage (N)
  nMeanBins <- max(1, nMeanBins)
  nMean <- rowMeans2(N(fds, type=psiType))
  quantilesBreaks <- quantile(nMean, p=(1:nMeanBins/nMeanBins))
  bins <- cut(rep(nMean, ncol(fds)), breaks=c(0, quantilesBreaks), right=TRUE, include.lowest=TRUE, ordered = TRUE)
  
  # nr of elements in logitPSI matrix
  nElem <- ifelse(psiType == "psiSite", nrow(mcols(fds, type="psiSite")) * ncol(fds), nrow(fds) * ncol(fds) )
  
  # get information about injected outliers
  if( assayExists(fds, paste0("trueOutliers_", psiType)) ){
    trueOutliers <- as.vector(getAssayMatrix(fds, "trueOutliers", psiType ))
  } else{
    trueOutliers <- double(nElem)
  }
  trueOutliers <- factor(trueOutliers, levels = c(-2,-1,0,1,2))
  
  # create table containing the data to plot
  plotData <- data.table(trueLogitPSI=logitPsi, 
                         predictedLogitPSI=logitPsi_fit,
                         difference=(abs(logitPsi - logitPsi_fit))^2,
                         method=rep(correction, nElem),
                         bin=bins,
                         outlier=trueOutliers,
                         correction=correction,
                         psiType=psiType,
                         dataset=dataset,
                         bestQ=bestQ(fds, psiType),
                         simQ=metadata(fds)[['optimalEncDim']])
  return(plotData)
  
}

# scatterplot true vs predicted logit(psi) (overall)
plotTrueVsPredictedLogitPsi <- function(dt, fds, dataset, psiType, correction){
  
  require(LSD)
  
  plotTitle <- paste0(dataset, ", ", psiType, ", sim_q = ", metadata(fds)[['optimalEncDim']], " best_q = ", unique(dt$bestQ), ", ", correction, ":\n simulated vs predicted logit(PSI)\n")
  
  #+ create plot
  par(mar=c(5,5,6,2))
  heatscatter(dt$trueLogitPSI, plotData$predictedLogitPSI,
    xlab="", ylab="", main=plotTitle, cor=TRUE); grid(); abline(0,1)
  axis(1, font=2)
  axis(2, font=2)
  mtext(side=1, line=3, bquote("simulated logit(" ~ Psi ~ ")"), font=2,cex=1.5)
  mtext(side=2, line=2, bquote("predicted logit(" ~ Psi ~ ")"), font=2,cex=1.5)
  #mtext(side=3, line=-1, plotTitle, font=2,cex=1.5)
}

# scatterplot true vs predicted logit(psi) (binned by outlier and coverage)
plotTrueVsPredictedLogitPsi_binned <- function(dt, fds, dataset, psiType, correction){
  
  dtc <- dt[, round(cor(trueLogitPSI, predictedLogitPSI), digits=3), by=c('bin', 'outlier')]
  
  require(ggplot2)
  
  plotTitle <- paste0(dataset, ", ", psiType, ", sim_q = ", unique(dt$simQ), ", best_q = ",  unique(dt$bestQ), ", ", correction, ":\n simulated vs predicted logit(PSI)")
  
  #+ create plot
  par(mar=c(1,1,1,1))
  g <- ggplot(dt, aes(x=trueLogitPSI, y=predictedLogitPSI)) + geom_point() + geom_abline(slope=1, intercept=0, color="red") +
    labs(x=bquote("simulated logit(" ~ Psi ~ ")"), y=bquote("predicted logit(" ~ Psi ~ ")"), title=plotTitle) + 
    facet_grid( bin ~ outlier, labeller = label_both) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio = 1) + 
    geom_text(data=dtc, aes(x=-5, y=4, label=paste("cor =", V1)))
  
  return(g)
  
}

# boxplot recovery of simulated logit psi
boxplotLogitPsiRecovery <- function(dt){
  
  # read infos from dt
  dataset <- unique(dt$dataset)
  psiType <- unique(dt$psiType)
  simQ <- unique(dt$simQ)
  bestQ <- unique(dt$bestQ)
  
  dt$bin <- factor(dt$bin, levels=sortIntervals(unique(dt$bin)))
  
  require(ggplot2)
  
    # helper function to print nr of data points in each bin
  give.n <- function(x){
    return(c(y = min(x)*1.1, label = length(x))) 
  }
  
  #+ create plot
  plotTitle <- paste0(dataset, ", ", psiType, ", sim_q = ", simQ, " , best_q = ", bestQ, ":\n recovery of simulated logit(PSI)")
  g <- ggplot(dt, aes(x=bin, y=difference, fill=correction)) + geom_boxplot() + 
    stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75)) + 
    scale_y_log10(breaks = function(x){ 10^seq(from=floor(log10(min(x))), to=ceiling(log10(max(x))), by = 2)}) +
    labs(x="mean junction coverage", y=bquote("| logit(" ~ Psi ~ ") - logit(" ~ hat(Psi) ~ ")|" ^2), title=plotTitle) + facet_grid( . ~ outlier, scales="free_x", labeller = label_both) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=15), axis.title = element_text(face = "bold", size=15), axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1), axis.text = element_text(face="bold")) 
  return(g)
    
}




