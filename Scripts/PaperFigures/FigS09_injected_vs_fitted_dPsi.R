#'---
#' title: Plot injected dPSI vs fitted dPSI(=|observed/data PSI - correctedPSI|)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - in_weighted:  '`sm expand(config["htmlOutputPath"] + "/OutlierInjection/{dataset}/{psiType}/{delta}/{{ae_impl}}_fit.html",            dataset=config["reference_dataset_injectedDpsiPlot"], psiType=config["psiTypes"], delta=config["inj_deltas"])`'
#'   - in_noWeights: '`sm expand(config["htmlOutputPath"] + "/OutlierInjection/{dataset}/{psiType}/{delta}/{{ae_impl}}-no-weights_fit.html", dataset=config["reference_dataset_injectedDpsiPlot"], psiType=config["psiTypes"], delta=config["inj_deltas"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS9_injectedVsFittedDpsi_{ae_impl}.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS9_injectedVsFittedDpsi_{ae_impl}.pdf"`'
#'   - wBhtml: 'Output/html/PaperFigures/injected_vs_fitted_dPsi/{ae_impl}.html'
#'  type: noindex
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(ae_impl="PCA-BB-Decoder")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS09_injected_vs_fitted_dPsi.R",
            wildcards=wildcards, options=options, rerun=TRUE)
    threads <- 25
}

#+ input
psiTypes   <- CONFIG$psiTypes
dataset    <- CONFIG$reference_dataset_injectedDpsiPlot
fitMethods <- paste0(snakemake@wildcards$ae_impl, c("", "-no-weights"))
register(MulticoreParam(5))

#+ output
outPng <- snakemake@output$outPng
outPdf <- snakemake@output$outPdf

#'
#' ## Plot injected dPSI vs fitted dPSI(=|observed/data PSI - correctedPSI|)
#'
#+ get data for plotting (injected vs fitted dpsi)
getPlotData <- function(fds, psiType){
    
    message(date(), ": getting plot data for ", psiType, " from ", workingDir(fds))
    
    # get injected outliers
    trueOut <- getAssayMatrix(fds, "trueOutliers", psiType)
    trueDpsi <- getAssayMatrix(fds, "trueDeltaPSI", psiType)
    ind <- which(trueOut != 0)
    
    # get corrected psi
    mu <- predictedMeans(fds, psiType)
    fitPsi <- mu[ind]
    
    # get data psi
    k <- K(fds, psiType)
    n <- N(fds, psiType)
    dataPsi <- (k+pseudocount())/(n+2*pseudocount())
    
    # get deltaPSI of fit
    deltaPSI <- dataPsi-mu
    
    # create data.table
    dt <- data.table(outlier=trueOut[ind], tDpsi=trueDpsi[ind], 
                     fDpsi=deltaPSI[ind])
    dt
}

getPlotPerPsiType <- function(psiType, dataset, fitMethods, datadir){
    message(date(), ": getting full plot data for type ", psiType)
    fdsFiles <- paste0(datadir, "/datasets/inject_uniformDistr/",
                       psiType, "/", fitMethods, "/savedObjects/", dataset,
                       "/fds-object.RDS")
    fds_weighted    <- loadFraseRDataSet(file=fdsFiles[1])
    fds_notWeighted <- loadFraseRDataSet(file=fdsFiles[2])
    
    # get data for weighted implementation
    dt_w <- getPlotData(fds_weighted, psiType=psiType)
    dt_w <- dt_w[,`:=`(psiType=psiType, method="weighted")]
    
    # get data for non-weighted implementation
    dt_nonW <- getPlotData(fds_notWeighted, psiType=psiType)
    dt_nonW <- dt_nonW[,`:=`(psiType=psiType, method="not weighted")]
    
    # combine and return
    dt <- rbind(dt_w, dt_nonW)
    return(dt)
}

dt_ls <- bplapply(psiTypes, getPlotPerPsiType, dataset=dataset, 
        fitMethods=fitMethods, datadir=CONFIG$DATADIR)
dt <- rbindlist(dt_ls)

#' 
#' Add density
#' 
#' * estimate point densities
#' * adapted from http://auguga.blogspot.com/2015/10/r-heat-scatter-plot.html

getDensity <- function(x,y, ...){
    # create density map
    dens <- kde2d(x, y, ...)
    gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model and apply model to get estimates for original data points
    mod <- loess(zgr~xgr*ygr, data=gr)
    ans <- predict(mod, newdata=data.frame(xgr=x, ygr=y))
    ans
}

dt[,density:=getDensity(tDpsi, fDpsi, h=0.1),by="method,psiType"]


#'
#' Create and save figure
#'
#+ create and save full figure
dt[,method:=gsub(" ",  "~", method)]
dt[,psiType:=factor(psiType, levels=c("psi5", "psi3", "psiSite"))]
levels(dt$psiType) <- c("psi[5]", "psi[3]", "theta")

g <- ggplot(data=dt, aes(x=fDpsi, y=tDpsi, color=density)) + 
    geom_point(alpha = 0.7, size=1) + 
    facet_grid(psiType ~ method, labeller=label_parsed) + 
    geom_smooth(method='lm') + 
    geom_abline(intercept=0, slope=1, col="firebrick") +
    scale_color_gradientn(colours=colorpalette('heat', 5)) + 
    ylab(expression("Injected" ~ Delta ~ "splice metric")) +
    xlab(expression(atop("Predicted" ~ Delta ~ "splice metric",
            "observed - predicted splice metric"))) + 
    theme_bw() + 
    theme(legend.position="none")
g


#+ save p value enrichment figure
factor <- 0.5
outPng
ggsave(outPng, g, width = 14*factor, height = 16*factor)
ggsave(outPdf, g, width = 14*factor, height = 16*factor)

