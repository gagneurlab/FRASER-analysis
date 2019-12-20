
#
# Inject artificial outliers into simulated data
#
simulateNInject <- function(simulation=c("BB", "DM"), nrJunctions, nrSamples, q, nrOutliers, type, verbose=FALSE){

  # Get simulated data
  fds <- switch(match.arg(simulation),
                DM = makeSimulatedFraserDataSet_Multinomial(m=nrSamples, j=nrJunctions, q=q),
                BB = makeSimulatedFraserDataSet_BetaBinomial(m=nrSamples, j=nrJunctions, q=q) )
  
  currentType(fds) <- type
  
  # Define deltaPSI to use in outlier injection
  # deltaPsi <- c(rep(0.5, nrOutliers/2), rep(-0.5, nrOutliers/2))
  # deltaPsi <- pmin(0.9, pmax(-0.9, rnorm(nrOutliers, 0, 0.8)))
  deltaPsi <- runif(nrOutliers, -0.9, 0.9)
  
  # Inject artificial outliers
  fds <- switch(match.arg(simulation),
                DM = injectOutliers(fds, type=type, nrOutliers=nrOutliers, deltaPSI=deltaPsi, method="simulatedPSI", verbose=verbose), 
                BB = injectOutliers(fds, type=type, nrOutliers=nrOutliers, deltaPSI=deltaPsi, method="simulatedPSI", swap = FALSE, verbose=verbose) ) 
  
  
  # save fds with injected outliers for later
  name(fds) <- paste0("outliers_", simulation, "_simulation")
  fds <- saveFraseRDataSet(fds)

  return(fds)
}

#
# Inject artificial outliers into existing data
#
injectInFDS <- function(fds, nrOutliers, type, injectionMethod=c("meanPSI", "samplePSI","simulatedPSI"), delta=c("normalDistr", "uniformDistr"), verbose=FALSE){
  
  currentType(fds) <- type
  
  # Define deltaPSI to use in outlier injection
  # deltaPsi <- c(rep(0.5, nrOutliers/2), rep(-0.5, nrOutliers/2))
  # deltaPsi <- pmin(0.9, pmax(-0.9, rnorm(nrOutliers, 0, 0.8)))
  deltaPsi <- switch(match.arg(delta),
                     normalDistr  = pmax(0,pmin(1, rnorm(nrOutliers, 0.5, 0.2))) * sample(c(-1,1), nrOutliers, replace=TRUE),
                     uniformDistr = runif(nrOutliers, -0.9, 0.9) )
  
  # Inject artificial outliers
  fds <- injectOutliers(fds, type=type, nrOutliers=nrOutliers, deltaPSI=deltaPsi, method=injectionMethod, verbose=verbose) 
  
  # # save fds with injected outliers for later
  # name(fds) <- paste0(name(fds), "_", injectionMethod)
  # fds <- saveFraseRDataSet(fds)
  
  return(fds)
}