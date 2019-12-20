#
# Get precision and recall for the different methods
#
computePrecRec <- function(fds){

  # Calculate precision recall values for the different methods
  plotData <- data.table()
  q <- metadata(fds)[['optimalEncDim']] # simulated data
  q <- 25 # real data
  
  # get FraseR results
  probE <- 1
  featureExclusionMask(fds, type=type) <- sample(c(TRUE, FALSE),
                                                 nrow(mcols(fds, type=type)), replace=TRUE, prob=c(probE,1-probE))
  plotData <- rbind(plotData, getResFraser(fds, q))
  # get PCA + rho fit results
  plotData <- rbind(plotData, getResPCA(fds, q))
  # get BB results
  plotData <- rbind(plotData, getResBB(fds))
  # get raw data dPsi results
  plotData <- rbind(plotData, getResDeltaPsi(fds))
  
  
  #
  # prepare data for plotting  
  #
  plotData$method <- factor(plotData$method, levels = c("FraseR", "BB", "PCA", "deltaPSI", "deltaPSI_FraseR", "deltaPSI_PCA"))
  plotData$dPsi_group <- factor(plotData$dPsi_group, levels = unique(plotData$dPsi_group))
  plotData$n_group <- factor(plotData$n_group, levels = unique(plotData$n_group))
  # plotData$recall_group[plotData$recall_group == "all"] <- "primary + secondary outliers"
  # plotData$recall_group[plotData$recall_group == "primary"] <- "only primary outliers"
  # plotData$recall_group <- factor(plotData$recall_group, levels = c("primary + secondary outliers", "only primary outliers"))

}

#
# Function to create a data table with precision and recall for an fds based on padj or delta psi values
#
precisionRecall <- function(fds, type, method=c('padj', 'BB', 'deltaPSI', 'predDeltaPSI'), trueOutliers=c("primary", "all", "secondary")){
  
  # get score, i.e either padj or delta psi values
  score <- switch(match.arg(method),
                  padj         = as.vector(padjVals(fds)),
                  BB           = as.vector(getAssayMatrix(fds, "padjBB", type)),
                  deltaPSI     = as.vector(abs(getAssayMatrix(fds, "delta", type))),
                  predDeltaPSI = abs(as.vector((K(fds, type=type) + pseudocount())/(N(fds, type=type) + 2*pseudocount())) - as.vector(t(predictMu(fds, type=type)))) )
  
  # get positions of simulated outliers
  simulatedOutliers <- as.vector(abs(getAssayMatrix(fds, "trueOutliers", type=type)))
  
  # dPsi group (< 0.3, 0.3-0.6, >0.6)
  dPsi <- as.vector(abs(getAssayMatrix(fds, "trueDeltaPSI", type=type)))
  dPsi_group <- cut(dPsi, breaks=c(0, 0.3, 0.6, 1), include.lowest = TRUE, ordered_result = TRUE)
  
  
  # adjust simulated outliers according to outlier type
  simulatedOutliers[simulatedOutliers == 1] <- switch(match.arg(trueOutliers),
                                                      primary   = 1,
                                                      secondary = 0,
                                                      all       = 1 )
  simulatedOutliers[simulatedOutliers == 2] <- switch(match.arg(trueOutliers),
                                                      primary   = 0,
                                                      secondary = 1,
                                                      all       = 1 )
  require(ROCR)
  # for all simulated outliers 
  res <- switch(match.arg(method),
                padj         = performance(prediction(1-score, simulatedOutliers), "prec", "rec"),
                BB           = performance(prediction(1-score, simulatedOutliers), "prec", "rec"),
                deltaPSI     = performance(prediction(score, simulatedOutliers), "prec", "rec"),
                predDeltaPSI = performance(prediction(score, simulatedOutliers), "prec", "rec") )
  cutoffs <- switch(match.arg(method),
                    padj         = 1-unlist(res@alpha.values),
                    BB           = 1-unlist(res@alpha.values), 
                    deltaPSI     = unlist(res@alpha.values),
                    predDeltaPSI = unlist(res@alpha.values) )
  plotData <- data.table(Precision = unlist(res@y.values), Recall = unlist(res@x.values), cutoff=cutoffs, dPsi_group="all", n_group="all", recall_group=trueOutliers)
  
  # for all outliers within each level of dPsi values
  for(group in levels(dPsi_group)){
    
    # group_simOutliers <- simulatedOutliers
    # group_simOutliers[dPsi_group != group] <- 0
    
    # exclude outliers from different psi groups
    group_outliers <- simulatedOutliers == 0 | dPsi_group == group
    group_simOutliers <- simulatedOutliers[group_outliers] 
    group_score <- score[group_outliers]
    
    res <- switch(match.arg(method),
                  padj         = performance(prediction(1-group_score, group_simOutliers), "prec", "rec"),
                  BB           = performance(prediction(1-group_score, group_simOutliers), "prec", "rec"),
                  deltaPSI     = performance(prediction(group_score, group_simOutliers), "prec", "rec"),
                  predDeltaPSI = performance(prediction(group_score, group_simOutliers), "prec", "rec") )
    cutoffs <- switch(match.arg(method),
                      padj         = 1-unlist(res@alpha.values),
                      BB           = 1-unlist(res@alpha.values), 
                      deltaPSI     = unlist(res@alpha.values),
                      predDeltaPSI = unlist(res@alpha.values) )
    plotData <- rbind(plotData, data.table(Precision = unlist(res@y.values), Recall = unlist(res@x.values), cutoff=cutoffs, dPsi_group=group, n_group="all", 
                                           recall_group=trueOutliers))
    
  }
  
  # group by N_outlier: only consider outliers with a<N<b, remove other outliers
  n <- N(fds, type=type)
  # get junction means of N (as mxn matrix to easily access it with outlier indices)
  mean_n <- rep(rowMeans(n), ncol(n))
  n_outlier <- mean_n[simulatedOutliers == 1]
  n_groups <- cut(mean_n, c(0, 100, 1000, max(n_outlier)), include.lowest=TRUE, ordered = TRUE)
  # n_groups <- cut(n, c(0, 10^(2:floor(log10(max(n_outlier)))), max(n_outlier)), include.lowest=TRUE, ordered = TRUE)
  
  for(n_group in levels(n_groups)){
    
    # for all simulated outliers with N in the range of this n_group
    nGroup_outliers <- simulatedOutliers == 0 | n_groups == n_group
    nGroup_simOutliers <- simulatedOutliers[nGroup_outliers]
    nGroup_score <- score[nGroup_outliers]
    res <- switch(match.arg(method),
                  padj         = performance(prediction(1-nGroup_score, nGroup_simOutliers), "prec", "rec"),
                  BB           = performance(prediction(1-nGroup_score, nGroup_simOutliers), "prec", "rec"),
                  deltaPSI     = performance(prediction(nGroup_score, nGroup_simOutliers), "prec", "rec"),
                  predDeltaPSI = performance(prediction(nGroup_score, nGroup_simOutliers), "prec", "rec") )
    cutoffs <- switch(match.arg(method),
                      padj         = 1-unlist(res@alpha.values),
                      BB           = 1-unlist(res@alpha.values),
                      deltaPSI     = unlist(res@alpha.values),
                      predDeltaPSI = unlist(res@alpha.values) )
    plotData <- rbind(plotData, data.table(Precision = unlist(res@y.values), Recall = unlist(res@x.values), cutoff=cutoffs, dPsi_group="all", n_group=n_group,
                                           recall_group=trueOutliers))
    
    # for all outliers within each level of dPsi values and this n_group
    for(group in levels(dPsi_group)){
      
      # group_simOutliers <- simulatedOutliers
      # group_simOutliers[dPsi_group != group] <- 0
      
      # exclude outliers from different psi groups
      group_outliers <- nGroup_simOutliers == 0 | dPsi_group[nGroup_outliers] == group
      group_simOutliers <- nGroup_simOutliers[group_outliers]
      group_score <- nGroup_score[group_outliers]
      
      res <- switch(match.arg(method),
                    padj         = performance(prediction(1-group_score, group_simOutliers), "prec", "rec"),
                    BB           = performance(prediction(1-group_score, group_simOutliers), "prec", "rec"),
                    deltaPSI     = performance(prediction(group_score, group_simOutliers), "prec", "rec"),
                    predDeltaPSI = performance(prediction(group_score, group_simOutliers), "prec", "rec") )
      cutoffs <- switch(match.arg(method),
                        padj         = 1-unlist(res@alpha.values),
                        BB           = 1-unlist(res@alpha.values),
                        deltaPSI     = unlist(res@alpha.values),
                        predDeltaPSI = unlist(res@alpha.values) )
      plotData <- rbind(plotData, data.table(Precision = unlist(res@y.values), Recall = unlist(res@x.values), cutoff=cutoffs, dPsi_group=group, n_group=n_group,
                                             recall_group=trueOutliers))
      
    }
    
    
  }
  
  
  return(plotData)
  
}


#
# Different methods to detect outliers (deltaPSI, PCA, FraseR, BB)
#

# data delta psi = data psi - mean(junction psi) + cutoff
getResDeltaPsi <- function(fds){
  
  type <- currentType(fds)
  
  fds <- calculatePSIValues(fds)
  # fds <- saveFraseRDataSet(fds)
  
  res_deltaPSI <- data.table()
  for(trueOutliers in c("all", "primary")){
    
    res_deltaPSI <- rbind(res_deltaPSI, precisionRecall(fds, type, method = "deltaPSI", trueOutliers=trueOutliers))
    
  }
  res_deltaPSI[, method:="deltaPSI"]
  return(res_deltaPSI)
}

# PCA + rho fitting and delta PSI PCA
getResPCA <- function(fds, q){
  
  type <- currentType(fds)
  
  featureExclusionMask(fds, type=type) <- rep(TRUE, nrow(mcols(fds, type=type)))
  rhoR <- c(1e-5, 1-1e-5)
  currentNoiseAlpha(fds) <- NULL
  fds <- initAutoencoder(fds, q=q, rhoRange=rhoR, type = type) # does E and D with PCA
  fds <- updateRho(fds, type=type, rhoRange=rhoR, BPPARAM=bpparam(), verbose=TRUE) # estimate rho
  predictedMeans(fds, type) <- t(predictMu(fds))
  
  fds <- calculatePvalues(fds)
  fds <- calculatePadjValues(fds)
  fds <- saveFraseRDataSet(fds)
  
  res <- data.table()
  for(trueOutliers in c("all", "primary")){
    
    res_pca <- precisionRecall(fds, type, method="padj", trueOutliers=trueOutliers)
    res_pca[, method:=rep("PCA", nrow(res_pca))]
    
    # delta psi pca = data psi - pca corrected psi + cutoff
    res_deltaPSI_pca <- precisionRecall(fds, type, method = "predDeltaPSI", trueOutliers=trueOutliers)
    res_deltaPSI_pca[, method:=rep("deltaPSI_PCA", nrow(res_deltaPSI_pca))]
    
    res <- rbind(res, res_pca, res_deltaPSI_pca)
    
  }
  return(res)
}

# BB
getResBB <- function(fds){
  
  type <- currentType(fds)
  # fds <- calculatePValues(fds)
  fds <- pvalueByBetaBinomialPerType(fds=fds, aname=paste0("pvalue_", type), psiType=type, pvalFun=betabinVglmTest) 
  fds <- saveFraseRDataSet(fds)
  # multiple testing correction?
  pvals <- replace_na(getAssayMatrix(fds, "pvalue", type), 1)
  padj <- apply(pvals, 2, p.adjust, method="BY")
  setAssayMatrix(fds, "padjBB", type) <- padj
  
  # todo check bb results for bb simulation
  
  res_bb <- data.table()
  for(trueOutliers in c("all", "primary")){
    
    res_bb <- rbind(res_bb, precisionRecall(fds, type, method="BB", trueOutliers=trueOutliers))
    
  }
  
  res_bb[, method:="BB"]
  return(res_bb)
  
}

# FraseR AE and delta PSI FraseR
getResFraser <- function(fds, q){
  
  type <- currentType(fds)
  
  print(table(featureExclusionMask(fds, type=type)))
  
  fds <- fitAutoencoder(fds, q=q, type=type, verbose=TRUE, iterations=15, BPPARAM = parallel(fds))
  
  fds <- calculatePvalues(fds)
  fds <- calculatePadjValues(fds)
  fds <- saveFraseRDataSet(fds)
  
  res <- data.table()
  for(trueOutliers in c("all", "primary")){
    
    res_fraser <- precisionRecall(fds, type, method="padj", trueOutliers=trueOutliers)
    res_fraser[, method:=rep("FraseR", nrow(res_fraser))]
    
    # delta psi pca = data psi - pca corrected psi + cutoff
    res_deltaPSI_fraser <- precisionRecall(fds, type, method = "predDeltaPSI", trueOutliers=trueOutliers)
    res_deltaPSI_fraser[, method:=rep("deltaPSI_FraseR", nrow(res_deltaPSI_fraser))]
    
    res <- rbind(res, res_fraser, res_deltaPSI_fraser)
    
  }
  return(res)
  
}
