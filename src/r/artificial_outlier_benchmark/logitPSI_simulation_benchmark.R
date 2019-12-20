
# Script for FraseR benchmark in recovering logit(psi) with simultated data
devtools::document()
devtools::load_all()


# Simulate data using betaBinomial and Dirichlet-Multinomial distribution
simulation <- "DM"
for(simulation in c("BB", "DM")){
  

  # generate simulated dataset
  m <- 200
  j <- 10000
  q <- 10
  fds <- switch(simulation,
                BB = makeSimulatedFraserDataSet_BetaBinomial(m=m, j=j, q=q),
                DM = makeSimulatedFraserDataSet_Multinomial(m=m, j=j, q=q))
  
  # fds object for pca
  fds_pca <- fds
  
  # fds object for peer
  fds_peer <- fds
  
  q <- metadata(fds)[['optimalEncDim']]
  # once for probE <- 1 and once for 0.1/0.5?
  probE <- 1
  # probE <- 0.5
  probE <- 0.1
  # type <- "psi5"
  
  for(type in c("psi5", "psiSite")){
    
    # set current type
    currentType(fds) <- type
    currentType(fds_pca) <- type
    currentType(fds_peer) <- type
    
    # set featureExclusionMask for pca object
    featureExclusionMask(fds_pca, type=type) <- rep(TRUE, nrow(mcols(fds_pca, type=type)))
    
    # set featureExclusionMask for peer object
    featureExclusionMask(fds_peer, type=type) <- rep(TRUE, nrow(mcols(fds_peer, type=type)))
    
    # subset fitting
    featureExclusionMask(fds, type=type) <- sample(c(TRUE, FALSE),
                                        nrow(mcols(fds, type=type)), replace=TRUE, prob=c(probE,1-probE))
    print(table(featureExclusionMask(fds, type=type)))
    
    # run autoencoder
    fds <- fitAutoencoder(fds, q=q, type=type, BPPARAM=parallel(fds), verbose=TRUE, iterations=15)
    
    # compare AE fit to true values
    K <- K(fds, type=type)
    N <- N(fds, type=type)
    psi <- (K + pseudocount())/(N + 2*pseudocount())
    
    trueMu <- getAssayMatrix(fds, "truePSI", type=type)
    trueRho <- mcols(fds, type=type)[["trueRho"]]
    # trueAlpha <- getAssayMatrix(fds, "trueAlpha", type=type)
    logitPsi <- as.vector(getAssayMatrix(fds, "trueLogitPSI", type=type))
    
    aeMu <- predictedMeans(fds, type=type)
    aeRho <- rho(fds, type=type)
    # aeAlpha <- psi * (1-aeRho)/aeRho
    logitPsi_dAE <- as.vector(predictY(fds))
  
    # par(mfrow=c(2,2))
    # heatscatter(log10(as.vector(trueAlpha)), log10(as.vector(aeAlpha)), xlab="log10(simulated alpha)", ylab="log10(AE fit alpha)", main=paste0(type, ": true vs AE alpha")); grid(); abline(0,1, col="black")
    # heatscatter(as.vector(trueMu), as.vector(aeMu), xlab="simulated mu", ylab="AE fit mu", main=paste0(type, ": true vs AE mu")); grid(); abline(0,1, col="black")
    # # heatscatter(as.vector(trueMu), as.vector(psi), xlab="simulated mu", ylab="data psi", main=paste0(type, ": true mu vs data psi")); grid(); abline(0,1, col="black")
    # # heatscatter(as.vector(psi), as.vector(aeMu), xlab=type, ylab="AE fit mu", main=paste(type, "vs AE mu")); grid(); abline(0,1, col="black")
    # # heatscatter(trueRho, aeRho, xlab="simulated rho", ylab="AE fit rho", main=paste0(type, ": true vs AE rho")); grid(); abline(0,1, col="black")
    # heatscatter(log10(trueRho), log10(aeRho), xlab="log10(simulated rho)", ylab="log10(AE fit rho)", main=paste0(type, ": log10(true vs AE rho)")); grid(); abline(0,1, col="black")
    # heatscatter(logitPsi, logitPsi_dAE, xlab="true logit PSI", ylab="dAE fit logit PSI", main = "true logit PSI vs dAE logit PSI"); grid(); abline(0,1)
    
    # do a PCA as comparison
    pca <- pca(as.matrix(x(fds_pca)), nPcs=q)
    pc <- pcaMethods::loadings(pca)
    D(fds_pca) <- pc
    E(fds_pca) <- pc
    b(fds_pca) <- double(nrow(mcols(fds, type=type)))
    # rho(fds_pca) <- methodOfMomemtsRho(K(fds_pca))
    pcaMu <- t(predictMu(fds_pca, type=type))
    logitPsi_pca <- as.vector(predictY(fds_pca))
    
    # par(mfrow=c(1,2))
    # heatscatter(as.vector(trueMu), as.vector(pcaMu), xlab="simulated mu", ylab="PCA fit mu", main=paste0(type, ": true vs PCA mu")); grid(); abline(0,1, col="black")
    # # heatscatter(as.vector(psi), as.vector(pcaMu), xlab=type, ylab="PCA fit mu", main=paste(type, "vs PcA mu")); grid(); abline(0,1, col="black")
    # heatscatter(logitPsi, logitPsi_pca, xlab="true logit PSI", ylab="PCA fit logit PSI", main = "true logit PSI vs PCA logit PSI"); grid(); abline(0,1)
  
    #  
    # Correct using PEER
    #
    require(peer)
    
    # q is the number of known hidden factors
    maxFactors <- q
    # default and recommendation by PEER: min(0.25*n, 100)
    # maxFactors <- min(as.integer(0.25* ncol(fds_peer)), 100)
    
    # prepare PEER model
    model <- PEER()
    PEER_setPhenoMean(model, as.matrix(x(fds_peer)))
    PEER_setNk(model, maxFactors)          # nr of hidden confounders
    PEER_setNmax_iterations(model, 1000)   # 1000 iterations is default
    # PEER_setAdd_mean(model, TRUE)        # should mean expression be added as additional factor?
    
    # run fullpeer pipeline
    PEER_update(model)
    
    # extract PEER data
    peerResiduals <- PEER_getResiduals(model)
    peerLogitMu <- t(as.matrix(x(fds_peer)) - peerResiduals)
    peerMu <- plogis(peerLogitMu)
    
    # save peer model in fds_peer object
    setAssayMatrix(fds_peer, "peerLogitMu", type=type) <- peerLogitMu
    metadata(fds_peer)[[paste0("PEERmodel_", type)]] <- list(
      alpha     = PEER_getAlpha(model),
      residuals = PEER_getResiduals(model),
      W         = PEER_getW(model),
      hiddenSpace = PEER_getX(model))
    
    # par(mfrow=c(1,2))
    # heatscatter(logitPsi, as.vector(peerLogitMu), xlab="true logit PSI", ylab="PEER fit logit PSI", main = "true logit PSI vs PEER logit PSI"); grid(); abline(0,1)
    # heatscatter(as.vector(trueMu), as.vector(peerMu), xlab="simulated mu", ylab="PEER fit mu", main = paste0(type, ": true vs PEER mu")); grid(); abline(0,1)
    
    # # SVA
    # require(sva)
    # 
    # mod <- model.matrix(~1, data=as.data.frame(x(fds_pca)))
    # dat <- as.matrix(t(x(fds_pca)))
    # sva_fit <- sva(dat, mod, method="two-step")
    # 
    # ## Regress out SVs to get logit(psi) prediction?
    # 
    
  }
  
  psi_nElem <- nrow(fds) * ncol(fds)
  se_nElem  <- nrow(mcols(fds, type="psiSite")) * ncol(fds)
  n_elem <- psi_nElem + se_nElem
  
  #logit PSI (true and predicted (AE, PCA)) for PSI
  logitPsi      <- as.vector(getAssayMatrix(fds, "trueLogitPSI", type="psi5"))
  logitPsi_dAE  <- as.vector(predictY(fds, type="psi5"))
  logitPsi_pca  <- as.vector(predictY(fds_pca, type="psi5"))
  logitPsi_peer <- as.vector(getAssayMatrix(fds_peer, "peerLogitMu", type="psi5"))
  
  # bin by N
  # N_psi <- as.vector(N(fds, type="psi5"))
  # bins <- cut(rep(N_psi, 3), breaks=c(0, 10^(0:(round(log10(max(N_psi)))-1)), max(N_psi)), include.lowest=TRUE, ordered = TRUE)
  N_psi <- N(fds, type="psi5")
  N_mean_psi <- rowMeans(N_psi)
  bins <- cut(rep(rep(N_mean_psi, ncol(fds)), 3), breaks=c(0, 10^(0:(round(log10(max(N_mean_psi)))-1)), max(N_mean_psi)), include.lowest=TRUE, ordered = TRUE)
  # nBins <- 5
  # bins <- cut(rep(rep(N_mean_psi, ncol(fds)), 3), breaks=c(0, quantile(N_mean_psi, p=(1:nBins)/nBins)), include.lowest=TRUE, ordered = TRUE)
  # bins <- cut(rep(rep(N_mean_psi, ncol(fds)), 3), breaks=nBins, include.lowest=TRUE, ordered = TRUE)
  
  logitPSI_diffs <- data.table(difference=c((abs(logitPsi - logitPsi_dAE))^2,(abs(logitPsi - logitPsi_pca))^2, (abs(logitPsi - logitPsi_peer))^2),
                                    method=c(rep('FraseR', psi_nElem), rep('PCA', psi_nElem), rep('PEER', psi_nElem)),
                                    bin=bins)
  
  # trueOutliers <- as.vector(getAssayMatrix(fds, "trueOutliers", type ))
  # logitPSI_diffs[, outlier:=rep(trueOutliers, 3)]
  # logitPSI_diffs$outlier <- factor(logitPSI_diffs$outlier, levels = c(-2,-1,0,1,2))
  
  #logit PSI (true and predicted (AE, PCA)) for SE
  logitSE       <- as.vector(getAssayMatrix(fds, "trueLogitPSI", type="psiSite"))
  logitSE_dAE   <- as.vector(predictY(fds, type="psiSite"))
  logitSE_pca   <- as.vector(predictY(fds_pca, type="psiSite"))
  logitSE_peer  <- as.vector(getAssayMatrix(fds_peer, "peerLogitMu", type="psiSite"))
  
  # bin by N
  # N_se <- as.vector(N(fds, type="psiSite"))
  # bins=cut(rep(N_se, 3), breaks=c(0, 10^(0:(round(log10(max(N_se)))-1)), max(N_se)), include.lowest=TRUE, ordered = TRUE)
  N_se <- N(fds, type="psiSite")
  N_mean_se <- rowMeans(N_se)
  bins=cut(rep(rep(N_mean_se, ncol(fds)), 3), breaks=c(0, 10^(0:(round(log10(max(N_mean_se)))-1)), max(N_mean_se)), include.lowest=TRUE, ordered = TRUE)
  # bins=cut(rep(rep(N_mean_se, ncol(fds)), 3), breaks=6, include.lowest=TRUE, ordered = TRUE)
    
  logitSE_diffs <- data.table(difference=c((abs(logitSE - logitSE_dAE))^2,(abs(logitSE - logitSE_pca))^2, (abs(logitSE - logitSE_peer))^2),
                                   method=c(rep('FraseR', se_nElem), rep('PCA', se_nElem), rep('PEER', se_nElem)),
                                   bin=bins)
  
  logit_diffs_full_binN <- rbind(logitPSI_diffs, logitSE_diffs)
  logit_diffs_full_binN[,psi_type:=c(rep("PSI_5", 3*psi_nElem), rep("SE", 3*se_nElem))]
  
  #
  # boxplot difference between between true and predicted logit PSI of AE, PCA and PEER, binned by N
  #
  require(ggplot2)
  plotTitle <- switch(simulation,
                      BB = "BetaBinomial simulation",
                      DM = "Dirichlet Multinomial simulation")
  plotTitle <- paste0(plotTitle, " (q = ", q, ", samples = ", m, ", junctions = ", j, ")")
  
  give.n <- function(x){
    return(c(y = min(x)*1.1, label = length(x))) 
    # experiment with the multiplier to find the perfect position
  }
  
  fig2a_binN <- ggplot(logit_diffs_full_binN, aes(x=bin, y=difference, fill=method)) + geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75)) + scale_y_log10() +
    labs(x="junction mean N", y="|logit(PSI) - logit(PSI_hat)|^2", title=plotTitle) + facet_grid( . ~ psi_type, scales="free_x") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename=paste0("fig2a_", simulation, "_binN.png"), plot = fig2a_binN, width=14, height=7)
  
  fig2a_outStrat <- ggplot(logitPSI_diffs, aes(x=bin, y=difference, fill=method)) + geom_boxplot() + 
    stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75)) + scale_y_log10() +
    labs(x="junction mean N", y="|logit(PSI) - logit(PSI_hat)|^2", title=plotTitle) + facet_grid( . ~ outlier, scales="free_x", labeller = label_both) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  fig2a_outStrat
  ggsave(filename=paste0("fig2a_", simulation, "_outlierStratified.png"), plot = fig2a_outStrat, width=25, height=10)
  
  # # 
  # # bin by logit psi
  # #
  # lowerBin_psi <- if(min(logitPsi) < -5){ c(min(logitPsi), -5) } else{ min(logitPsi) }
  # bins <- cut(rep(logitPsi, 3), c(lowerBin_psi, -1, -0.1, 0.1, 1, 5, max(logitPsi)), include.lowest = TRUE, ordered = TRUE)
  # logitPSI_diffs <- data.table(difference=c((abs(logitPsi - logitPsi_dAE))^2,(abs(logitPsi - logitPsi_pca))^2, (abs(logitPsi - logitPsi_peer))^2),
  #                              method=c(rep('FraseR', psi_nElem), rep('PCA', psi_nElem), rep('PEER', psi_nElem)),
  #                              bin=bins)
  # 
  # lowerBin_se <- if(min(logitSE) < -5){ c(min(logitSE), -5) } else{ min(logitSE) }
  # bins <- cut(rep(logitSE, 3), c(lowerBin_se, -1, -0.1, 0.1, 1, 5, max(logitSE)), include.lowest = TRUE, ordered = TRUE)
  # logitSE_diffs <- data.table(difference=c((abs(logitSE - logitSE_dAE))^2,(abs(logitSE - logitSE_pca))^2, (abs(logitSE - logitSE_peer))^2),
  #                             method=c(rep('FraseR', se_nElem), rep('PCA', se_nElem), rep('PEER', se_nElem)),
  #                             bin=bins)
  # 
  # logit_diffs_full_binLPsi <- rbind(logitPSI_diffs, logitSE_diffs)
  # logit_diffs_full_binLPsi[,psi_type:=c(rep("PSI_5", 3*psi_nElem), rep("SE", 3*se_nElem))]
  # 
  # # boxplot difference between between true and predicted logit PSI of AE, PCA and PEER, binned by logitPsi
  # fig2a_binLPsi <- ggplot(logit_diffs_full_binLPsi, aes(x=bin, y=difference, fill=method)) + geom_boxplot() + 
  #   stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75)) + scale_y_log10() +
  #   labs(x="logit(PSI)", y="|logit(PSI) - logit(PSI_hat)|^2", title=plotTitle) + facet_grid( . ~ psi_type, scales="free_x") + 
  #   theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  # # fig2a_binLPsi
  # ggsave(filename=paste0("fig2a_", simulation, "_binLogitPsi.png"), plot = fig2a_binLPsi, width=14, height=7)
  

  # # plot fit of logit PSI for AE, PCA, PEER
  # par(mfrow=c(2,3))
  # # plot logit psi for PSI5 and AE, PCA, PEER
  # heatscatter(logitPsi, logitPsi_dAE, xlab="true logit PSI", ylab="AE fitted logit PSI", main="PSI5: true vs AE logit PSI"); grid(); abline(0,1, col="black")
  # heatscatter(logitPsi, logitPsi_pca, xlab="true logit PSI", ylab="PCA fitted logit PSI", main="PSI5: true vs PCA logit PSI"); grid(); abline(0,1, col="black")
  # heatscatter(logitPsi, logitPsi_peer, xlab="true logit PSI", ylab="PEER fitted logit PSI", main="PSI5: true vs PEER logit PSI"); grid(); abline(0,1, col="black")
  # # plot logit psi for SE  and AE, PCA, PEER
  # heatscatter(logitSE, logitSE_dAE, xlab="true logit SE", ylab="AE fitted logit SE", main="SE: true vs AE logit SE"); grid(); abline(0,1, col="black")
  # heatscatter(logitSE, logitSE_pca, xlab="true logit SE", ylab="PCA fitted logit SE", main="SE: true vs PCA logit SE"); grid(); abline(0,1, col="black")
  # heatscatter(logitSE, logitSE_peer, xlab="true logit SE", ylab="PEER fitted logit SE", main="SE: true vs PEER logit SE"); grid(); abline(0,1, col="black")

}

