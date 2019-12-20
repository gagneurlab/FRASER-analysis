#'---
#' title: Define Simulation dataset (BB & DM)
#' author: Christian Mertes
#' wb:
#'  output:
#'  - fdsobjBB:  '`sm config["DATADIR"] + "/datasets/savedObjects/raw-SimulationBB/fds-object.RDS"`'
#'  - countsJBB: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-SimulationBB/rawCountsJ.h5"`'
#'  - countsSBB: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-SimulationBB/rawCountsSS.h5"`'
#'  - fdsobjDM:  '`sm config["DATADIR"] + "/datasets/savedObjects/raw-SimulationDM/fds-object.RDS"`'
#'  - countsJDM: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-SimulationDM/rawCountsJ.h5"`'
#'  - countsSDM: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-SimulationDM/rawCountsSS.h5"`'
#' output:
#'  html_document
#'---

outFile <- "Data/paperPipeline/datasets/savedObjects/raw-SimulationBB/fds-object.RDS"
outFileBB <- snakemake@output$fdsobjBB
outFileDM <- snakemake@output$fdsobjDM

#+ load main config, echo=FALSE
source("./src/r/config.R", echo=FALSE)

# dataset name
workingDir <- dirname(dirname(dirname(outFileBB)))
nameBB <- basename(dirname(outFileBB))
nameDM <- basename(dirname(outFileDM))

#'
#' # Create Simulation dataset based on BetaBinomial distributio
#'
nSamples   <- 200
nJunctions <- 5000
simuQ      <- 10

#'
#' Dataset Name
#' and working directory
#+ echo=TRUE
nameBB
nameDM
workingDir

#'
#' ## Get simulations
#'
fdsBB <- makeSimulatedFraserDataSet_BetaBinomial(m=nSamples, j=nJunctions,
        q=simuQ, name=nameBB, workingDir=workingDir)

fdsDM <- makeSimulatedFraserDataSet_Multinomial(m=nSamples, j=nJunctions,
        q=simuQ, name=nameDM, workingDir=workingDir)

#'
#' ## save objects
#'
fdsBB <- saveFraseRDataSet(fdsBB)
fdsDM <- saveFraseRDataSet(fdsDM)


# pcaobj <- pca(t(xmat_rc[biggest_sds,]), nPcs=20)
# pcaobj
# plot(pcaMethods::scores(pcaobj)[,1], pcaMethods::scores(pcaobj)[,2])
