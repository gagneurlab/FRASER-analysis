
#
# Functions for generating Precision-Recall plots
#

# Full figure plot
fig2b <- function(plotData, dataset, nrJunctions, nrSamples, nrInjOutliers){
  
  require(ggplot2)
  
  title <- paste("Outlier detection benchmark for", dataset, "(", nrJunctions,"x", nrSamples, ") with", nrInjOutliers, "injected outliers")
  
  # remove NAs from data
  plotData$Precision[is.nan(plotData$Precision)] <- 0
  
  # plot all/primary outlier recall x deltaPsi group 
  fullPrecRecPlot <- ggplot(plotData, aes(x=Recall, y=Precision, color=method)) + geom_line() + 
    facet_grid(dPsi_group ~ recall_group, labeller = label_bquote(Delta * Psi == .(as.vector(dPsi_group)))) + 
    labs(title=title) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  return(fullPrecRecPlot)
}

# Full figure plot stratified by N and delta psi
fig2b_byN <- function(plotData, dataset, group=c('all', 'primary'), nrJunctions, nrSamples, nrInjOutliers){
  
  require(ggplot2)
  
  title <- paste(group, "outlier detection benchmark for", dataset, "(", nrJunctions,"x", nrSamples, ") with", nrInjOutliers, "injected outliers")
  
  toPlot <- plotData[recall_group == group]
  
  # remove NAs from data
  toPlot$Precision[is.nan(toPlot$Precision)] <- 0
  
  # plot all/primary outlier recall x deltaPsi group 
  fullPrecRecPlot <- ggplot(toPlot, aes(x=Recall, y=Precision, color=method)) + geom_line() + 
    facet_grid(dPsi_group ~ n_group, 
               labeller = label_bquote(rows = Delta * Psi == .(as.vector(dPsi_group)), cols = "mean junction N" == .(as.vector(n_group))) ) + 
    labs(title=title) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  return(fullPrecRecPlot)
}
