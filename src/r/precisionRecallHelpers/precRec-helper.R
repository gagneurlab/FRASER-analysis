
# Create table with precision recall results
createPrecRecTable <- function(fds, psiType, BPPARAM=bpparam(),
                    outliers=c('byJunctionGroup', 'all','primary'), nMeanBins=3,
                    dPsiBins=3, q, correction, anTrueCts='originalCounts',
                    Nsamples='All', inj_value='deltaPSI=uniformDistr', 
                    pmethod='BY', dataset='dataset'){
    # check if we go for grouping of results
    groupRes <- outliers == 'byJunctionGroup' 
    useIndex    <- TRUE
    if(isTRUE(groupRes)){
        useIndex <- !duplicated(getSiteIndex(fds, psiType))
    }
    
    # start the cluster if appropriate
    if(!bpisup(BPPARAM)){
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }
    
    # get counts 
    message(date(), " Extracting raw counts (k and n) ...")
    k <- as.matrix(K(fds, psiType)[useIndex,])
    n <- as.matrix(N(fds, psiType)[useIndex,])
    nAlternatives <- data.table(index=getSiteIndex(fds, psiType))[,
            .(rowID=.I, .N), by=index][order(rowID)][useIndex]
    
    message(date(), " Computing deltaPSI and rank ...")
    # neg (-) for descending order]
    deltaPSI <- getDeltaPsi(fds, psiType, byGroup=groupRes, BPPARAM=BPPARAM)
    deltaPSI_rank <- frank(c(-abs(deltaPSI)))
    
    message(date(), " Getting z scores ...")
    zscores <- zScores(fds, type=psiType, byGroup=groupRes, BPPARAM=BPPARAM)
    
    message(date(), " Getting pvalues and true outliers ...")
    trueOutlier <- getTrueOutliers(fds, psiType, byGroup=groupRes, BPPARAM=BPPARAM)
    trueDPsi    <- getTrueDeltaPsi(fds, psiType, byGroup=groupRes, BPPARAM=BPPARAM)
    pvals       <- pVals(fds, psiType, byGroup=groupRes)
    padjust     <- padjVals(fds, psiType, byGroup=groupRes)
    
    if(outliers == "primary"){
        trueDPsi[abs(trueOutlier) == 2] <- 0
        trueOutlier[abs(trueOutlier) == 2] <- 0
    } else if(outliers == 'all'){
        trueOutlier[trueOutlier == 2] <- 1
        trueOutlier[trueOutlier == -2] <- -1
    }
    totalTrueOut <- sum(trueOutlier != 0)
    
    message(date(), " Computing deltaPSI bins ...")
    dPsiBins <- max(1, dPsiBins)
    quantilesBreaks <- quantile(abs(trueDPsi[trueDPsi != 0]), 
            p=(1:dPsiBins/dPsiBins))
    deltaPsiBins <- cut(as.vector(abs(trueDPsi)), 
            breaks=c(min(abs(trueDPsi[trueDPsi != 0])), 
                    quantilesBreaks), 
            right=TRUE, include.lowest=TRUE, ordered_result=TRUE, dig.lab=2)
    
    
    # junction coverage bin
    message(date(), " Computing junction mean bins ...")
    nMeanBins <- max(1, nMeanBins)
    rMeans <- rowMeans2(n)
    
    quantilesBreaks <- quantile(rMeans, p=(1:nMeanBins/nMeanBins))
    quantilesBreaks[3] <- pretty(quantilesBreaks[3])[2]
    meanBins <- cut(rMeans, breaks=c(0, quantilesBreaks), right=TRUE, 
            include.lowest=TRUE, ordered_result=TRUE, dig.lab=2)
    
    nZeros <- rowSums(n == 0)
    
    if(q == "best"){
        q <- bestQ(fds, psiType)
    }
    
    message(date(), " Creating data table ...")
    ans <- data.table(
        rowInAssay=seq_along(getSiteIndex(fds, psiType))[useIndex],
        colInAssay=factor(rep(samples(fds), each=nrow(n))),
        numInjected = totalTrueOut,
        trueOutlier = c(trueOutlier),
        pvalues = c(pvals),
        padjust = c(padjust),
        zScore = c(zscores),
        dPsi = c(deltaPSI),
        dPsi_rank = c(deltaPSI_rank),
        k = c(k),
        n = c(n),
        Nsamples=factor(Nsamples),
        inj=factor(outliers),
        inj_value=factor("all"),
        dPsiBin=factor(deltaPsiBins),
        injDpsi=c(trueDPsi),
        correction=factor(correction),
        pmethod=factor(pmethod),
        dataset=factor(dataset),
        q=factor(q),
        geneMean=rMeans,
        geneMeanBin=factor(meanBins),
        junctionMeanBin=factor("all"),
        nAlternatives=nAlternatives[,N],
        nZeros=nZeros)
    
    # bin by N and dPsi
    if(dPsiBins > 1 || nMeanBins > 1){
        message(date(), " Binning by deltaPSI and junction mean ...")
        bpstop(BPPARAM)
        gc()
        ansls <- lapply(levels(ans$geneMeanBin), d=ans, function(nBin, d){
            ls <- lapply(levels(ans$dPsiBin), d=d, nBin=nBin, function(dpsiBin, d, nBin){
                message(date(), ": working on: ", dpsiBin, " and ", nBin)
                
                d1 <- d[(dPsiBin == dpsiBin & geneMeanBin == nBin) | trueOutlier == 0,]
                d1[,numInjected:=sum(abs(trueOutlier))]
                d1[,dPsi_rank:=frank(c(-abs(dPsi)))]
                d1[,inj_value:= dpsiBin]
                d1[,junctionMeanBin:= nBin]
                setorder(d1, dPsi_rank)
                d1[,recall := cumsum(abs(trueOutlier))/numInjected]
                d1
            })
        })
        
        # reset result list since we do have it now mixed 
        resultls <- list()
        gc()
        for(i in seq_along(ansls)){
            resultls <- append(resultls, ansls[[i]])
        }
    } else {
        resultls <- list(ans)
    }
    
    return(resultls)
    
}

pvalue_score_sign <- function(pvalue, dPsi, cutoff=0.1){
    score <- pvalue
    score[abs(dPsi) < cutoff] <- 1
    score
}

readDPSICutoff <- function(f, cutoffs=c(0.2,0.5,0.8)){
    dt <- fread(f)
    dt[,score:=score_sign(dPsi)]
    ans <- lapply(cutoffs, function(x){
        ans <- dt[score > x]
        if(nrow(ans) == 0){
            return(dt[,.(recall=0, precision=0, Nsamples, inj, inj_value, junctionMeanBin,
                         correction, pmethod, type = "DeltaPSI rank", cutoff=x)][1])
        }
        setorder(ans, -score)
        ans[, .(
            recall=sum(trueOutlier != 0)/unique(numInjected),
            precision=sum(trueOutlier != 0)/.N,
            Nsamples = unique(Nsamples),
            inj = unique(inj),
            inj_value = unique(inj_value),
            junctionMeanBin = unique(junctionMeanBin),
            correction = unique(correction),
            pmethod = unique(pmethod),
            type = "DeltaPSI rank",
            cutoff=x)]
    })
    rbindlist(ans)
}

readZscoreCutoff <- function(f, cutoffs=c(2,3,5), dPsiCut=0.1){
    if(isScalarCharacter(f)){
        dt <- fread(f)
    } else {
        dt <- f
    }
    dt[,score:=score_sign(zScore, abs(dPsi), dPsiCut, 0)]
    ans <- lapply(cutoffs, function(x){
        ans <- dt[score > x]
        if(nrow(ans) == 0){
            return(dt[,.(recall=0, precision=0, Nsamples, inj, inj_value, 
                    junctionMeanBin, correction, pmethod, 
                    type = "zScore rank", cutoff=x)][1])
        }
        setorder(ans, -score)
        ans[, .(
            recall=sum(trueOutlier != 0)/unique(numInjected),
            precision=sum(trueOutlier != 0)/.N,
            Nsamples = unique(Nsamples),
            inj = unique(inj),
            inj_value = unique(inj_value),
            junctionMeanBin = unique(junctionMeanBin),
            correction = unique(correction),
            pmethod = unique(pmethod),
            type = "zScore rank",
            cutoff=x)]
    })
    rbindlist(ans)
}

readPadjCutoff <- function(f, cutoffs=c(0.001, 0.01, 0.05, 0.1), FDR_LIMIT=0.1, dPsiCutoff=0.1){
    dt <- fread(f)
    dt[, score:=pvalue_score_sign(pvalues, dPsi, dPsiCutoff)]
    dt[, score_adjust:=pvalue_score_sign(padjust, dPsi, dPsiCutoff)]
    dt[, label:=trueOutlier != 0]
    dt <- dt[score_adjust < FDR_LIMIT,]
    
    dt <- dt[order(score, decreasing=FALSE)]
    dt[,rank      := 1:.N]
    dt[,TP        := cumsum(label)]
    dt[,c('TP', 'rank'):=list(max(TP), max(rank)), by=score]
    dt[,precision := TP/rank]
    dt[,recall    := TP/max(numInjected)]
    
    ans <- lapply(cutoffs, function(x){
        ans <- dt[score_adjust < x]
        if(nrow(ans) == 0){
            return(dt[,.(recall=0, precision=0, Nsamples, inj, inj_value, junctionMeanBin,
                         correction, pmethod, 
                         type = paste('P-value rank & deltaPSI >', 
                                      dPsiCutoff, '\n& FDR <', FDR_LIMIT), cutoff=x)][1])
        }
        ans[, .(
            recall=ans[nrow(ans),recall],
            precision=ans[nrow(ans),precision],
            Nsamples = unique(Nsamples),
            inj = unique(inj),
            inj_value = unique(inj_value),
            junctionMeanBin = unique(junctionMeanBin),
            correction = unique(correction),
            pmethod = unique(pmethod),
            type = paste('P-value rank & FDR <', FDR_LIMIT, '\n& deltaPSI >', dPsiCutoff),
            cutoff=x)]
    })
    rbindlist(ans)
}

getBootstrapScore <- function(dt, pval, dpsi, maxRows, dPsiCutoff=0.1){
    if(isTRUE(dpsi)){
        dt[, score:=score_sign(dPsi, values4Cut=dPsi, cutoff=dPsiCutoff, set2z=0)]
    } else if(isTRUE(pval)){
        dt[, score:= -pvalue_score_sign(pvalues, dPsi=dPsi, cutoff=dPsiCutoff)]
    } else  {
        dt[, score:=score_sign(zScore, values4Cut=dPsi, cutoff=dPsiCutoff, set2z=0)]
    }
    dt <- dt[order(score, decreasing=TRUE)][label == TRUE | 1:.N <= maxRows]
    dt
}


#
# gather data
#
readBootstrapData <- function(f, rankType=c("pvalue", "zscore"), 
                    maxRows=1e6, FDR_LIMIT=0.1, dPsiCutoff=0.1, mc.cores=10,
                    zScore_LIMIT=NA){
    rankType <- match.arg(rankType, c("pvalue", "zscore", "deltaPsi"), 
            several.ok=TRUE)
    
    if(is.data.table(f)){
        dt <- f
        f <- paste0('.../', unique(dt$correction), '/...')
    } else {
        dt <- fread(f)
    }
    dt[, label:=trueOutlier != 0]
    ans <- list()
    
    # 
    # zscore based curve
    # 
    if("zscore" %in% rankType){
        dtz <- getBootstrapScore(dt, FALSE, FALSE, maxRows, dPsiCutoff)
        currType <- 'zScore rank'
        if(!is.na(zScore_LIMIT)){
            dtz <- dtz[abs(score) > zScore_LIMIT]
            currType <- paste0("zScore rank & abs(zScore) > ", zScore_LIMIT)
        }
        dtz <- dtz[,.(bootstrapPR(score, label, total=max(numInjected), 
                n=200, mc.cores=mc.cores), cutoff=zScore)]
        dtz[,type:=currType]
        ans <- append(ans, list(dtz))
    }
    
    # 
    # pvalue based curve
    # 
    if("pvalue" %in% rankType){
        dtp <- getBootstrapScore(dt, TRUE, FALSE, maxRows, dPsiCutoff)
        currType <- 'P-value rank'
        if(!is.na(FDR_LIMIT)){
            dtp <- dtp[padjust < FDR_LIMIT]
            currType <- paste('P-value rank & FDR <', FDR_LIMIT)
        }
        dtp <- dtp[,.(bootstrapPR(score, label, total=max(numInjected),
            n=200, mc.cores=mc.cores), cutoff=padjust)]
        dtp[,type:=currType]
        ans <- append(ans, list(dtp))
    }
    
    # 
    # deltaPsi based curve
    # 
    if("deltaPsi" %in% rankType){
        stop("not implemented yet. Please check precRec-helper for more details")
        # dt3 <- getBootstrapScore(dt, !zscoreForAll & pvalue,  TRUE, maxRows, dPsicutoff)
        # dt3 <- dt3[,bootstrapPR(score, label, total=max(numInjected), n=200, mc.cores=mc.cores)]
        
        # dt3 <- rbind(dt3, data.table(score=-Inf, label=FALSE, rank=Inf,
        #                              TP=max(dt3$TP), precision=0, recall=1, lower=0, upper=0))
        # dt3 <- rbind(dt3, dt3[min(rank)==rank][1][,.(
        #     score, label, rank=0, TP, precision, recall=0, lower, upper)])
        # dt3[, type:='DeltaPSI rank']
    }
    
    # 
    # combine results and add additional informations
    # 
    ans <- rbindlist(lapply(ans, function(x){
        x[,.(type, rank, TP, precision, recall, lower, upper, cutoff)]}))
    
    ans$correction      <- unique(dt$correction)
    ans$q               <- unique(dt$q)
    ans$inj_value       <- unique(dt$inj_value)
    ans$inj             <- unique(dt$inj)
    ans$junctionMeanBin <- unique(dt$junctionMeanBin )
    if(!is.na(dPsiCutoff) & dPsiCutoff > 0){
        ans[,type:=paste0(type, '\n& deltaPSI >', dPsiCutoff)]
    }
    
    # return it
    ans[!is.na(cutoff)]
}

plotRibbonBenchmark <- function(data, title, linetype=c(1,2,3), maxRows=1e5,
                                FDR_LIMIT=0.1,
                                wrap_function=function() facet_grid(inj_value ~ junctionMeanBin),
                                zscoreData=NULL, confidenceIntervals=TRUE){
    #' Reduce plotting points
    prProbs <- c(10000/maxRows, 1-10000/maxRows)
    data <- data[rank < 10000 | sample(x=c(TRUE, FALSE), size=.N, replace=TRUE, prob=prProbs)]
    #data <- data[order(type, correction, inj_value, inj, precision, recall)]
    #dup <- duplicated(data[,.(type, correction, inj_value, inj, precision, recall)])
    #data <- data[!dup]
    # print(head(data))
    print(DT::datatable(rbind( head(data[type == "DeltaPSI rank",], n = 1000), head(data[type == paste0("P-value rank \n& FDR < ", FDR_LIMIT),], n = 1000) )))
    
    
    #' Plot it
    gg <- ggplot(data, aes(recall, precision, color=correction, linetype=type)) +
        geom_line()
    if(isTRUE(confidenceIntervals)){
        gg <- gg + geom_ribbon(alpha=0.2, col=NA, aes(x=recall, fill=correction,
                                                      linetype=type,
                                                      ymin = pmax(0, pmin(1, lower)),
                                                      ymax = pmax(0, pmin(1, upper))))
    }
    gg <- gg + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
        scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
        scale_linetype_manual(values=linetype, labels = function(lab){
            sapply(lab, function(lab){
                switch(lab,
                       'DeltaPSI rank' = bquote(Delta*Psi ~ "rank"),
                       lab )} )
        }) +
        wrap_function() +
        grids(linetype="dashed") +
        theme(strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm")),
              strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
        theme(legend.key.height=unit(0.8, "cm")) +
        labs(x='Recall', y='Precision', color="Correction",
             fill='Correction', linetype="Method")
    if(!missing(title)){
        gg <- gg + ggtitle(title)
    }
    if(!is.null(zscoreData)){
        gg <- gg +
            geom_point(data=zscoreData,
                       aes(recall, precision, shape=factor(cutoff)), size=3) +
            scale_shape_manual(values=seq_len(length(levels(factor(zscoreData$cutoff))))) +
            labs(shape=bquote("padjust/zScore/" ~ Delta * Psi ~ " cutoff"))
    }
    gg
}

score_sign <- function(zScore, values4Cut, cutoff, set2z=0){
    score <- zScore
    score <- abs(score)
    if(!missing(values4Cut)){
        if(!(is.null(cutoff) | is.na(cutoff) | cutoff == 0)){
            score[values4Cut < cutoff] <- set2z
        }
    }
    return(score)
}

renameCorrectionMethods <- function(dt, column){
    naming <- c(
        PCA               = 'pca',
        PEER              = 'peer',
        "Cook's"          = 'baseCooks',
        Pearson           = 'basePearsonRes',
        OUTRIDER          = CONFIG_YAML$FIGURE_IMPLEMENTATION,
        'robust OUTRIDER' = CONFIG_YAML$FIGURE_ROB_IMPLEMENTATION)
    
    x <- dt[,as.character(get(column))]
    for(i in seq_along(naming)){
        x <- gsub(paste0('^', naming[i], '$'), names(naming)[i], x)
    }
    dt[,c(column):=list(x)]
    
    return(dt)
}

correctPrecisionRankPlotNames <- function(dt){
    # Injection name
    dt[inj == 'all',     inj:='Primary + secondary outliers']
    dt[inj == 'primary', inj:='Only primary outliers']
    
    # Injection value
    dt[,inj_value:=gsub('deltaPSI=', 'Simulated |dPSI| = ', inj_value)]
    
    dt[type == 'DeltaPSI rank', type:='|dPSI| ranked']
    dt[type == 'DeltaPSI rank', type:='|dPSI| ranked']
    dt[type == 'FDR filtered', type:=paste('Ranked & FDR <', FDR_LIMIT)]
    
    return(dt)
}

sortIntervals <- function(intervals, decr=FALSE){
    if(!is.character(intervals)){
        intervals <- as.character(intervals)
    }
    split <- strsplit(intervals, ",")
    intervalStart <- as.numeric(substring(lapply(split , `[`, 1), 2))
    order <- order(intervalStart, decreasing = decr)
    sorted <- intervals[order]
    return(sorted)
}

