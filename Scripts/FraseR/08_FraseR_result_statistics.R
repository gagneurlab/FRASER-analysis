#'---
#' title: Get Result Statistics
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fdsin:       '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/pajdBetaBinomial_psiSite.h5"`'
#'   - rawfdsin:    '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/delta_psiSite.h5"`'
#'   - resultTable: '`sm config["DATADIR"] + "/processedData/results/{dataset}/{method}_results.tsv"`'
#'   - gtf:         '`sm config["GTF_FILE"]`'
#'  output:
#'   - stats:  '`sm config["DATADIR"] + "/processedData/results/{dataset}/{method}_stats.RDS"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_final_stats.html"`'
#'  type: noindex
#'  threads: 3
#'---

#+ load config and setup, echo=FALSE
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Minor_Salivary_Gland", method="PCA")
    parseWBHeader2("./Scripts/FraseR/08_FraseR_result_statistics.R",
            wildcards, rerun=TRUE)
}

#+ input
dataset    <- snakemake@wildcards$dataset
method     <- snakemake@wildcards$method
gtfFile    <- snakemake@input$gtf
fdsFile    <- snakemake@input$fdsin
outFile    <- snakemake@output$stats

workingDir <- dirname(dirname(dirname(fdsFile)))
name       <- basename(dirname(fdsFile))

FDR_LIMIT          <- 0.1
DELTA_PSI_LIMIT    <- 0.3
Z_SCORE_LIMIT      <- 2
MIN_COVERAGE       <- 10
CORRELATION_DIGITS <- 2
#'
#' # Load data
#'
#' Dataset:
#+ echo=TRUE
dataset
method
workingDir

#+ echo=FALSE
fds <- loadFraserDataSet(dir=workingDir, name=name)
fds_raw <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))
BPPARAM <- MulticoreParam(10, 10)
register(BPPARAM)

ans <- list()

#' 
#' Annotate ranges first
#' 
txdb <- makeTxDbFromGFF(gtfFile)
fds <- findAnnotatedJunction(fds, annotation=txdb, annotateNames=FALSE)

knownStarts <- mcols(fds, type="psi5")[
        !duplicated(getSiteIndex(fds, "psi5")), "known_start"]
knownEnds <- mcols(fds, type="psi3")[
    !duplicated(getSiteIndex(fds, "psi3")), "known_end"]

#' ## Important numbers
ans['ImportantNumbers'] <- list(c(
    NStarts       = length(knownStarts),
    NEnds         = length(knownEnds),
    knownStarts   = sum(knownStarts),
    knownEnds     = sum(knownEnds),
    NJunctions    = nrow(assay(fds, "psi5")),
    NSites        = nrow(assay(fds, "psiSite")),
    NRawJunctions = nrow(assay(fds_raw, "psi5")),
    NRawSites     = nrow(assay(fds_raw, "psiSite")),
    Nsamples      = ncol(fds)))
for(type in psiTypes){
    tmp <- c(bestQ(fds, type), bestNoise(fds, type))
    names(tmp) <- paste0(c("Q_", "Noise_"), type)
    ans[['ImportantNumbers']] <- c(ans[['ImportantNumbers']], tmp)
}
ans['ImportantNumbers']

#' 
#' Percentage of novel splice sites (donor and acceptor)
#' 
round((1 - ans[['ImportantNumbers']]['knownStarts'] / 
        ans[['ImportantNumbers']]['NStarts']) * 100, 2)
round((1 - ans[['ImportantNumbers']]['knownEnds']   / 
        ans[['ImportantNumbers']]['NEnds'])   * 100, 2)

#'
#' ## Added value of intron retention detection
#'
addedValIR_ls <- lapply(psiTypes, function(type){
        pv  <- pVals(fds, type=type)
        dpv <- deltaPsiValue(fds, type=type)
        n   <- N(fds, type=type)
        
        hgnc_symbols <- mcols(fds, type=type)[,"hgnc_symbol"]
        pv[abs(dpv) < DELTA_PSI_LIMIT] <- 1
        pv[n < MIN_COVERAGE] <- 1
    
        pvpg <- getPvalsPerGene(fds, type=type, pvals=pv)
        pvpg })

asplice <- addedValIR_ls[[1]] < FDR_LIMIT | addedValIR_ls[[2]] < FDR_LIMIT
irsplice <- addedValIR_ls[[3]] < FDR_LIMIT

goodGeneIDs <- intersect(rownames(asplice), rownames(irsplice))

as_num <- sum(asplice[goodGeneIDs,], na.rm=TRUE)
ir_num <- sum(irsplice[goodGeneIDs,], na.rm=TRUE)
only_ir_num <- sum(irsplice[goodGeneIDs,] & !asplice[goodGeneIDs,], na.rm=TRUE)
ir_increase <- (only_ir_num + as_num) / as_num * 100

as_num
ir_num
ir_increase

ans[['ImportantNumbers']]['AltSplicing'] <- as_num
ans[['ImportantNumbers']]['IRSplicing']  <- ir_num
ans[['ImportantNumbers']]['IR_increase'] <- ir_increase

#'
#' # computing result numbers for different cutoff approaches
#' 
getAberrantNumbers <- function(fds, pac, dpc, zsc, mcc){
    rbindlist(mclapply(psiTypes, function(x){
        tmp1 <- aberrant(fds, type=x, padjCutoff=pac, deltaPsiCutoff=dpc, 
                zScoreCutoff=zsc, minCoverage=mcc, aggregate=TRUE)
        tmp1 <- quantile(colSums(tmp1), probs=c(0.1, 0.5, 0.9, 1))
        as.data.table(tmp1, keep.rownames=TRUE)[,.(quantile=rn, outlier=tmp1, 
                type=x, padjCut=pac, dPsiCut=dpc, zsCut=zsc, minCov=mcc)]
    }))
}

resNumbers <- list()
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=NA,  dpc=NA,  zsc=2, mcc=10)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=NA,  dpc=0.1, zsc=2, mcc=10)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=NA,  dpc=0.3, zsc=2, mcc=10)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=NA,  dpc=NA,  zsc=3, mcc=10)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=0.1, dpc=NA,  zsc=0, mcc=10)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=0.1, dpc=0.1, zsc=0, mcc=10)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=0.1, dpc=0.3, zsc=0, mcc=10)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=0.1, dpc=NA,  zsc=0, mcc=5)))
resNumbers <- append(resNumbers, list(
        getAberrantNumbers(fds, pac=NA,  dpc=NA,  zsc=2, mcc=5)))

allResNumbers <- rbindlist(resNumbers)
allResNumbers[,tissue:=dataset]
allResNumbers

ans[["allResNumbers"]] <- allResNumbers

#' 
#' ## Sample correlation
#'
correlations <- list()
corDT <- rbindlist(lapply(psiTypes, function(type){
    minMedian=1
    logit=TRUE
    topN=100000

    # get data counts
    kmat <- as.matrix(K(fds, type=type))
    nmat <- as.matrix(N(fds, type=type))

    # subset to most variable junctions
    expRowsMedian <- rowMedians(kmat) >= minMedian
    expRowsMax    <- rowMax(kmat)     >= 20

    skmat <- kmat[expRowsMax & expRowsMedian,]
    snmat <- nmat[expRowsMax & expRowsMedian,]
    xmat  <- (skmat + 1)/(snmat + 2)
    xmat  <- pmin(pmax(round(xmat, CORRELATION_DIGITS), 
              10^-CORRELATION_DIGITS),
            1-10^-CORRELATION_DIGITS)
    xmat_sd <- rowSds(xmat)
    plotIdx <- rank(xmat_sd) >= length(xmat_sd) - topN
    xmat <- xmat[plotIdx,]

    xmat <- qlogis(xmat)
    xmat_rc <- xmat - rowMeans(xmat)

    # get raw sample correlations
    cormatRaw <- cor(xmat_rc)

    # get sample correlation after normalization
    pred_mu <- qlogis(as.matrix(predictedMeans(fds, type=type)[
            expRowsMax & expRowsMedian,][plotIdx,]))
    lpred_mu_rc <- pred_mu - rowMeans(pred_mu)
    norm_xmat_rc <- xmat_rc - lpred_mu_rc
    cormatNorm <- cor(norm_xmat_rc)

    # create data table
    rbind(
        data.table(cor=cormatRaw[upper.tri(cormatRaw)],
                normalized="raw", dataset=dataset, method=method, type=type),
        data.table(cor=cormatNorm[upper.tri(cormatNorm)],
                normalized="normalized", dataset=dataset, method=method, type=type))
}))


ggplot(corDT, aes(y=abs(cor), x=normalized, fill=type)) + geom_boxplot()
corDT[,.(mean=mean(abs(cor)), sd=sd(abs(cor))),by="normalized,type"]

ans["SampleCors"] <- list(corDT)
ans


#' # Save numbers to RDS file
#'
saveRDS(ans, outFile)

