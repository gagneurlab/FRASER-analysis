#'
#' # DROP bug analysis 
#' 
devtools::load_all("../fraser")

#' input
fdsFile <- "/s/project/drop-analysis/processed_data/aberrant_splicing/datasets/savedObjects/all/fds-object.RDS"

#' load object
fds <- loadFraseRDataSet(file=fdsFile)

#' show object
fds
colData(fds)

#' extract results
res <- results(fds)
res
plotAberrantPerSample(fds)

plotCountCorHeatmap(fds, type="psi5", normalized=FALSE, logit=TRUE)
plotCountCorHeatmap(fds, type="psi5", normalized=TRUE, logit=TRUE)

bestQ(fds, type="psi3")
bestQ(fds, type="psi3")
bestQ(fds, type="psiSite")

#' problematic junction
jid <- 42850
type <- "psi3"
plotExpectedVsObservedPsi(fds, type=type, idx=jid)
plotQQ(fds, type=type, idx=jid)


res
i <- 1
pdf("drop_splicing.pdf")
    for(i in 1:10){
        message(date(), ": work on res: ", i)
        print(ggarrange(nrow=2, ncol=2,
            plotVolcano(fds, as.character(res[i]$sampleID), 
                    as.character(res[i]$type), aggregate=TRUE),
            plotExpression(fds, result=res[i]),
            plotExpectedVsObservedPsi(fds, result=res[i]),
            plotQQ(fds, result=res[i])))
    }
dev.off()
