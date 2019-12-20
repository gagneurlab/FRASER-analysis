#
# FraseR shiny example
#
# @author: Christian Mertes
#

# load FraseR package
library(devtools)
load_all("../fraser")
useRealWorldData <- TRUE

# create example data set
if(useRealWorldData){
    # load precomputed real world data
    kerberos.init()
    path <- "/s/project/fraser/analysis/prokisch-full"
    name <- "Full_Prokisch_Analysis"
    fds <- loadFraseRDataSet(path, name)
    res <- results(fds)

} else {
    # extract toy data
    fds <- createTestFraseRSettings()
    parallel(fds) <- MulticoreParam(4)

    # run the stats
    fds <- countRNAData(fds)
    fds <- calculatePSIValues(fds)
    fds <- calculateZScores(fds)
    fds <- calculatePValues(fds)

    # annotate data with ENSEMBL
    fds <- annotateRanges(fds)

    # extract results
    res <- results(fds, fdrCut = 0.5, zscoreCut = 0.5)
}

# start shiny server
FraseRShinyApp(fds, res)

