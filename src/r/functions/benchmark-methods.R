#'
#' get ROC objects
#'
simple_roc <- function(labels, scores, decreasing=TRUE){
    labels <- labels[order(scores, decreasing=decreasing)]
    data.frame(labels,
               TPR=cumsum(labels)/sum(labels),
               FPR=cumsum(!labels)/sum(!labels))
}

simple_recall <- function(labels, scores, topN=100, decreasing=TRUE){
    labels <- labels[order(scores, decreasing=decreasing)]
    labels <- labels[1:min(length(labels), topN)]
    sum(labels)/length(labels)
}

rankVsHits <- function(labels, scores, decreasing=FALSE){
    labelOrder <- order(scores, decreasing = decreasing)
    labels <- labels[labelOrder]
    scores <- scores[labelOrder]
    labelsum <- sum(labels)
    if(labels[length(labels)]){
        maxVal <- rev(unique(scores))[1:2]
        if(diff(maxVal) <= -1){
            if(all(labels[which(scores == maxVal[1])])){
                labels <- labels[scores != maxVal[1]]
            }
        }
    }
    cumsum(labels)
}

getRocObj <- function(dataName, data, predictionName, rocOptions, FUN=roc,
                allSamples=FALSE, topN=NA, ...){
    currOpt <- rocOptions[algo==dataName]

    # get real data
    currRes <- data[[dataName]]
    if(isFALSE(allSamples)){
        currRes <- currRes[RNA_ID %in% benchmark_set[,RNA_ID]]
    }
    if(!is.na(currOpt[,psiChangeCut])){
        currRes <- currRes[psiChange>currOpt[,psiChangeCut]]
    }
    currRes[,prediction:=get(predictionName)]
    if(predictionName != "pvalue"){
        currRes[,prediction:=1/(prediction + 1e-10)]
        currRes <- currRes[,.(prediction=min(prediction, na.rm=TRUE)),
                           by="RNA_ID,hgnc_symbol"]
    } else {
        currRes <- currRes[,.(prediction=min(p.adjust(prediction, "hochberg"),
                                             na.rm=TRUE)), by="RNA_ID,hgnc_symbol"]
    }


    if(is.na(topN) | topN == 0){
        topN <- length(currRes$prediction)
    }
    currRes <- currRes[order(prediction)][1:min(nrow(currRes), topN)]

    # merge real data with benchmark set
    benchmark2merge <- benchmark_set[,.(RNA_ID, hgnc_symbol, event=TRUE)]
    currRes <- merge(currRes, benchmark2merge, all=TRUE)
    currRes[is.na(event),event:=FALSE]
    maxVal <- max(currRes[,prediction], na.rm=TRUE)
    currRes[is.na(prediction), prediction:=maxVal+1]

    if(nrow(currRes) > topN){
        currRes <- currRes[order(prediction)]
        callOrders <- currRes[,c(which(event), which(!event))]
        currRes <- currRes[sort(callOrders[1:topN])]
    }

    return(FUN(currRes$event, currRes$prediction, ...))
}


#'
#' plot ROC curves
#'
plotFraseRRoc <- function(dataName, rocObjs, rocOptions){
    currOpt <- rocOptions[algo==dataName]
    currRoc <- rocObjs[[dataName]]

    # plot simple way
    if(is(currRoc, "data.frame")){
        if(dataName == "FraseR"){
            plot(NA, ylim=c(0,1), xlim=c(1,0), main="ROC curve",
                 xlab="FPR / 1 - Specificity", ylab="TPR / Sensitivity")
            grid()
            abline(1,-1,col="gray")
        }
        lines(x=c(1,1-currRoc[,"FPR"]), y=c(0,currRoc[,"TPR"]),
              col=currOpt[,color])

        #currRoc <- rocObjs[[3]]
        #lines(col="red",
        #    x=currRoc[,2],
        #    y=rev(currRoc[,1]))
        #a<-1

        # plot pROC way
    } else {
        plot.roc(currRoc, add=dataName != "FraseR", grid=TRUE,
                 plot.auc=FALSE, col=currOpt[,color])
        aucText <- paste0(currOpt[,displayName],
                          ": AUC: ", round(currRoc$auc, 3),
                          "; N = ", length(currRoc$predictor))
        text(0.4, currOpt[,aucY]-0.05, aucText, col=currOpt[,color])
    }
}
