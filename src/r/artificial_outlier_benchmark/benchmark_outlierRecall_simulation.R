
devtools::load_all("~/Projects/FraseR/FraseR/")
source("Projects/FraseR-analysis/src/r/artificial_outlier_benchmark/injectOutliers.R")
source("Projects/FraseR-analysis/src/r/artificial_outlier_benchmark/computePrecisionRecall.R")
source("Projects/FraseR-analysis/src/r/artificial_outlier_benchmark/figure2b.R")


#
# Benchmark of recall of artificially injected outliers into simulated data
#

# Define params
simulation <- "BB"
nrOut <- 500
# simulation <- "DM"
# nrOut <- 250
type <- "psi5"

# Get simulated data with injected outliers and save fds
fds <- simulateNInject(simulation=simulation, nrJunctions = 10000, nrSamples = 100, q = 10, nrOutliers = nrOut, type=type, verbose = TRUE)
dataset <- paste0("outliers_", simulation, "_simulation")
# or load existiing simulated data
# fds <- loadFraseRDataSet("~/FraseR/", dataset)

plotData <- computePrecRec(fds, dataset)

# save plotData 
saveRDS(plotData, paste0(name(fds), '_PrecRecData.RDS'))

# # restore old counts
# removeInjectedOutliers(fds, currentType(fds))

#
# Plot figure 2B (precision-recall already computed)
#

# read in data for plotting
plotData <- readRDS(paste0(dataset, '_PrecRecData.RDS'))
nrInjOut <- sum(abs(getAssayMatrix(fds, "trueOutliers", currentType(fds))) == 1)

# get name of data for plot title (simulated or real data respectively)
datasetInfos <- strsplit(name(fds), "_")[[1]]
datasetName <- paste0(datasetInfos[2], " (", datasetInfos[3], ")")

# create and save plot
# fullPrecRecPlot <- fig2b(plotData, simulation, nrow(fds), ncol(fds), nrInjOut)
for(trueOutliers in c('all', 'primary')){
  
  fullPrecRecPlot <- fig2b_byN(plotData, datasetName, trueOutliers, nrow(fds), ncol(fds), nrInjOut)
  # fullPrecRecPlot
  # save plot
  ggsave(filename=paste0("fig2b_", datasetName, "_", trueOutliers, ".png"), plot = fullPrecRecPlot, width=15, height=12)
}


