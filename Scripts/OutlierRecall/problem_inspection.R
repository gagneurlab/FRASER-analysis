#'---
#' title: Inspect p-value problems in FraseR fit
#' author: Ines Scheller
#' wb:
#'  input:
#'   - inHtml: 'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_pvalues.html' 
#'  output:
#'   - wBhtml:  'Output/html/OutlierInjection/{dataset}/{psiType}/{delta}/{method}_fitProblems.html' 
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
dataset    <- snakemake@wildcards$dataset
psiType    <- snakemake@wildcards$psiType
delta      <- snakemake@wildcards$delta
fitMethod  <- snakemake@wildcards$method
workingDir <- file.path(CONFIG$DATADIR, "datasets", paste0("inject_",delta), psiType, fitMethod)

#' ## Load dataset
#+ load fds
dataset
workingDir
fds     <- loadFraseRDataSet(workingDir, dataset)

#' ## Fit parameters
#+ get fit params
fitMethod
psiType


#' ## Get pvalues, mu, ...
#+ get data
pvals <- pVals(fds, psiType)
trueOut <- getAssayMatrix(fds, "trueOutliers", psiType)
mu <- predictedMeans(fds, psiType)
k  <- K(fds, psiType)
n  <- N(fds, psiType)
# x  <- x(fds, psiType)
E  <- E(fds, psiType)
H  <- H(fds, psiType)
D  <- D(fds, psiType)
b  <- b(fds, psiType)
rho <- rho(fds, psiType)

#' ## Get junctions where mu==1
#+ get problematic junctions
mu1 <- which(mu == 1, arr.ind=TRUE)
mu1_junctions <- unique(mu1[,1])
mu1_junctions

#' ## Get params for those junctions
#' b:
b[mu1_junctions]
#' rho:
rho[mu1_junctions]

#' ## K vs N for those junctions (green: injected outliers)
#+ plot problematic junctions
dev_null <- sapply(mu1_junctions, function(i){
  print(paste("mu fit of junction", i))
  print(summary(mu[i,]))
  print("mu == 1?:")
  print(table(mu[i,] == 1))
  print("K==N where mu==1?")
  print(table(k[i,which(mu[i,] == 1)] == n[i,which(mu[i,] == 1)]))
  print("H * D[i,]:")
  HD <- as.vector(H %*% D[i,])
  print(summary(HD))
  print("y[i,]:")
  print(summary(HD + b[i]))
  print(paste("b =", b[i]))
  print(paste("rho =", rho[i]))
  
  outlier <- trueOut[i,]
  pos <- which(outlier != 0)
  plotData(idx=i, k, n, pvals=NULL)
  points(n[i,pos] + 2*pseudocount(), k[i, pos] + pseudocount(), pch=2, col="green")
    
  plot(n[i,]-k[i,]+1, log = "y", xlab="sample", ylab="log10(N-K+1)")
  plot(1:ncol(fds), mu[i,], xlab = "sample", ylab="mu fit")
  # plot(1:ncol(fds), (k[i,]+1)/(n[i,]+2), xlab = "sample", ylab="data psi")
  heatscatter((k[i,]+1)/(n[i,]+2), mu[i,], xlab = "data psi", ylab="mu fit", main=paste("mu (raw data vs fit), junction =", i))
})


