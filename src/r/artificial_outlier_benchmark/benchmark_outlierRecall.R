
# devtools::load_all()
# source("injectOutliers.R")
# source("getPrecisionRecall.R")
devtools::load_all("~/Projects/FraseR/FraseR/")
source("Projects/FraseR-analysis/src/r/artificial_outlier_benchmark/injectOutliers.R")
# source("Projects/FraseR-analysis/src/r/artificial_outlier_benchmark/getPrecisionRecall.R")


#
# Benchmark of recall of artificially injected outliers into real data (Kremer dataset + GTEx skin)
#

# load fds and copy to another dir
# fds <- loadFraseRDataSet("/s/project/fraser/analysis/datasets", "kremer-bader-et-al", upgrade=TRUE)
fds <- loadFraseRDataSet("~/Projects", "kremer-bader-et-al", upgrade=TRUE)
workingDir(fds) <- "~/FraseR"
name(fds) <- "outliers_kremer-bader-et-al"
fds <- saveFraseRDataSet(fds)

# define params for outlier injection
nrOut <- 500
type <- "psi5"
method <- "meanPSI"
method <- "samplePSI"

# inject outliers and save fds
fds <- injectRealData(fds=fds, nrOutliers=nrOut, type=type, injectionMethod=method, verbose=TRUE)


dataset <- "outliers_kremer-bader-et-al_meanPSI"
