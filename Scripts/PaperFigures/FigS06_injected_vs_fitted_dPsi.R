#'---
#' title: Plot injected dPSI vs fitted dPSI(=|observed/data PSI - correctedPSI|)
#' author: Ines Scheller
#' wb:
#'  params:
#'   - geom_hex: FALSE
#'   - binwidth: 0.02
#'  input:
#'   - in_weighted:  '`sm expand(config["DATADIR"] + "/datasets/inject_{delta}/{psiType}/{AE}/savedObjects/{dataset}/predictedMeans_{psiType}.h5", dataset=config["reference_dataset_injectedDpsiPlot"], AE=config["AE_IMPLEMENTATION"], psiType=config["psiTypes"], delta=config["inj_deltas"])`'
#'   - in_noWeights: '`sm expand(config["DATADIR"] + "/datasets/inject_{delta}/{psiType}/{AE}/savedObjects/{dataset}/predictedMeans_{psiType}.h5", dataset=config["reference_dataset_injectedDpsiPlot"], AE=config["nonWeighted_AE_IMPLEMENTATION"], psiType=config["psiTypes"], delta=config["inj_deltas"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS6_injectedVsFittedDpsi.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS6_injectedVsFittedDpsi.pdf"`'
#' output:
#'  html_document
#'---


#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=10, fig.height=10)

#+ input
fdsFiles_weighted    <- snakemake@input$in_weighted
fdsFiles_notWeighted <- snakemake@input$in_noWeights
psiTypes      <- CONFIG$psiTypes
dataset       <- CONFIG$reference_dataset_injectedDpsiPlot
fitMethods    <- c(CONFIG$AE_IMPLEMENTATION, CONFIG$nonWeighted_AE_IMPLEMENTATION)
use_geom_hex  <- as.logical(snakemake@params[[1]]$geom_hex)
binwidth      <- as.numeric(snakemake@params[[1]]$binwidth) 

#+ output
outPng <- snakemake@output$outPng
outPdf <- snakemake@output$outPdf

#'
#' ## Plot injected dPSI vs fitted dPSI(=|observed/data PSI - correctedPSI|)
#'
#+ get data for plotting (injected vs fitted dpsi)
getPlotData <- function(fds, psiType){
    
    message(date(), ": getting plot data for ", psiType, " from ", workingDir(fds), "/savedObjects/", name(fds))
    
    # get  injected outliers
    trueOut <- getAssayMatrix(fds, "trueOutliers", psiType)
    trueDpsi <- getAssayMatrix(fds, "trueDeltaPSI", psiType)
    # all injected outliers
    ind <- which(trueOut != 0)
    # # only injected primary outliers
    # ind <- which(abs(trueOut) > 0)
    
    # get corrected psi
    mu <- predictedMeans(fds, psiType)
    
    # get data psi
    k <- K(fds, psiType)
    n <- N(fds, psiType)
    dataPsi <- (k+pseudocount())/(n+2*pseudocount())
    
    # get deltaPSI of fit
    deltaPSI <- dataPsi-mu
    
    # get n==0 frequency at junctions with injections
    nrNzero <- matrix(nrow=nrow(trueOut), ncol=ncol(trueOut), rowSums(n == 0))
    
    # create data.table
    dt <- data.table(outlier=trueOut[ind], tDpsi=trueDpsi[ind], 
                     fDpsi=deltaPSI[ind], nrNzeroSamples=nrNzero[ind])
    
    dt
}

getPlotPerFile <- function(fdsFile, weighted){
    
    message(date(), ": reading in ", fdsFile)
    fds    <- loadFraseRDataSet(file=fdsFile)
    psiType <- unlist(strsplit(unlist(
                                strsplit(basename(fdsFile), "_"))[2], "[.]")
                                )[1]
    if(isTRUE(weighted)){
        method <- "weighted"
    } else {
        method <- "not weighted"
    }
    
    # get data for weighted implementation
    dt <- getPlotData(fds, psiType=psiType)
    dt <- dt[, psiType:=psiType]
    dt <- dt[, method:=method]
    
    # estimate point densities
    # found here: https://slowkow.com/notes/ggplot2-color-by-density/
    
    # Get density of points in 2 dimensions.
    # @param x A numeric vector.
    # @param y A numeric vector.
    # @param n Create a square n by n grid to compute density.
    # @return The density within each square.
    get_density <- function(x, y, ...) {
        dens <- MASS::kde2d(x, y, ...)
        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ii <- cbind(ix, iy)
        return(dens$z[ii])
    }
    dt[,pointdens:=get_density(tDpsi, fDpsi, n = 100)]
    
    return(dt)
}

dtls_w <- lapply(fdsFiles_weighted, getPlotPerFile, weighted=TRUE)
dt_w <- rbindlist(dtls_w)

dtls_nonW <- lapply(fdsFiles_notWeighted, getPlotPerFile, weighted=FALSE)
dt_nonW <- rbindlist(dtls_nonW)

dt <- rbind(dt_w, dt_nonW)

# # for the plot: irgnore injections at junctions with >50 samples with n==0
# dt <- dt[nrNzeroSamples < 50]

#'
#' Create and save figure
#'
#+ create and save full figure
g <- ggplot(data=dt, aes(x=tDpsi, y=fDpsi)) + facet_grid(psiType ~ method)

if(isFALSE(use_geom_hex)){
    g <- g + geom_point(aes(color=pointdens), size=1, show.legend=FALSE) + 
        scale_color_gradientn(colors = colorpalette('heat', 30))
} else {
    # alternative using geom_hex
    g <- g + geom_hex(aes(fill=stat(log(count))), binwidth = binwidth, 
                  show.legend = FALSE) + 
    scale_fill_gradientn(colors = colorpalette('heat', 30))
    
}
    
g <- g + geom_smooth(method='lm') + 
    geom_abline(intercept=0, slope=1, col="firebrick") + 
    labs(x=bquote("injected" ~ Delta * Psi),
         y=bquote(Delta * Psi ~ "= observed" ~ Psi - "predicted" ~ Psi)) +
    theme_bw()
g

factor <- 0.5
ggsave(outPng, g, width = 17*factor, height = 14*factor)
ggsave(outPdf, g, width = 17*factor, height = 14*factor)
