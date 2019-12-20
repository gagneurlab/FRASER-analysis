#'---
#' title: Create tissue sample correlation boxplots
#' author: Ines Scheller
#' wb:
#'  params:
#'   - workers: 10
#'  input:
#'   - inputTissues: '`sm expand("Output/html/FraseR/{tissue}/FraseR-1DecoderBatches_autoencoder_fit.html", tissue=config["EnrichmentTissues"])`'
#'  output:
#'   - outData: '`sm expand(config["FIGDATADIR"] + "/heatmaps/boxplot_{psiType}.RDS", psiType=config["psiTypes"])`'
#' output:
#'  html_document
#'---

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/DataAnalysis/sampleCorrelations.R", dataset="Skin_Not_Sun_Exposed_Suprapubic")
}

#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=30, fig.height=10)

#+ input
tissues     <- CONFIG$EnrichmentTissues
workingDir  <- file.path(CONFIG$DATADIR, "datasets")
correctionMethod <- "FraseR-1DecoderBatches"
bpWorkers   <- min(bpworkers(), as.integer(snakemake@params[[1]]$workers))

#+ output
outData <- snakemake@output$outData

#'
#' Boxplot of sample correlations of raw and normalized counts for different GTEx tissues
#'
plots <- lapply(psiTypes, function(type){

    message(date(), ": ", type)

    #+ create plot data for one psiType
    dtls <- bplapply(tissues, function(tissue, type, minMedian=1, logit=TRUE, topN=100000){

        # load fds
        fds <- loadFraseRDataSet(workingDir, paste0(tissue, "__", correctionMethod))

        # get data counts
        kmat <- as.matrix(counts(fds, type=type, side="ofIn"))
        nmat <- kmat + as.matrix(counts(fds, type=type, side="other"))

        expRowsMedian <- rowMedians(kmat) >= minMedian
        expRowsMax    <- rowMax(kmat)     >= 10
        table(expRowsMax & expRowsMedian)

        skmat <- kmat[expRowsMax & expRowsMedian,]
        snmat <- nmat[expRowsMax & expRowsMedian,]

        xmat <- (skmat + 1)/(snmat + 2)
        if(isTRUE(logit)){
            xmat <- qlogis(xmat)
        }
        xmat_rc    <- xmat - rowMeans(xmat)

        # get raw sample correlations
        xmat_rc_sd <- rowSds(xmat_rc)
        plotIdx <- rank(xmat_rc_sd) >= length(xmat_rc_sd) - topN
        xmat_rc_2_plot <- xmat_rc[plotIdx,]
        cormatRaw <- cor(xmat_rc_2_plot)

        # get sample correlation after normalization
        pred_mu <- as.matrix(predictedMeans(fds, type=type)[
            expRowsMax & expRowsMedian,][plotIdx,])
        if(isTRUE(logit)){
            pred_mu <- qlogis(pred_mu)
        }
        lpred_mu_rc <- pred_mu - rowMeans(pred_mu)
        xmat_rc_2_plot <- xmat_rc_2_plot - lpred_mu_rc
        cormatNorm <- cor(xmat_rc_2_plot)

        # create data table
        dtRaw  <- data.table(cor=cormatRaw[upper.tri(cormatRaw)], normalized="raw", tissue=tissue)
        dtNorm <- data.table(cor=cormatNorm[upper.tri(cormatNorm)], normalized="normalized", tissue=tissue)
        dt <- rbind(dtRaw, dtNorm)

        dt

    }, type=type, BPPARAM=MulticoreParam(bpWorkers))

    dt2p <- rbindlist(dtls)
    dt2p$normalized <- factor(dt2p$normalized, levels=c("raw", "normalized"))

    g <- ggplot(dt2p, aes(x=tissue, y=abs(cor), fill=normalized)) + geom_boxplot() +
        theme_bw() + xlab("GTEx tissue") + ylab("|sample correlation|") +
        ggtitle(paste0("Sample correlations in GTEx tissues (", type, ")")) +
        labs(fill="")
        # + theme(legend.position = c(0.925, 0.08))
    saveRDS(g, file=outData[which(grepl(x=outData, pattern=type))])

    g

})

require(gridExtra)
gall <- grid.arrange(arrangeGrob(grobs= plots,ncol=1))
gall

