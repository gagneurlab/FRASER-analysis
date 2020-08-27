#'---
#' title: Get Result Statistics
#' author: Ines Scheller
#' wb:
#'  input:
#'   - fdsin:       '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/pajdBetaBinomial_psiSite.h5"`'
#'   - resultTable: '`sm config["DATADIR"] + "/processedData/results/{dataset}/{method}_results.tsv"`'
#'   - gtf:         '`sm config["GTF_FILE"]`'
#'  output:
#'   - stats:  '`sm config["DATADIR"] + "/processedData/results/{dataset}/{method}_stats_minN.RDS"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_stats_minN.html"`'
#'  type: noindex
#'  threads: 3
#'---

#+ load config and setup, echo=FALSE
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Minor_Salivary_Gland", method="PCA")
    parseWBHeader2("./Scripts/FraseR/09_FraseR_sampleCorrelations_minN.R",
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

FDR_LIMIT       <- 0.1
DELTA_PSI_LIMIT <- 0.3
Z_SCORE_LIMIT   <- 2
MIN_COVERAGE    <- 10

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
BPPARAM <- MulticoreParam(10, 10)
register(BPPARAM)

ans <- list()

#' 
#' ## Sample correlation with min N
#'
minN <- c(0, 1, 10, 50, 100)
combinations <- CJ(psiTypes, minN, unique=TRUE)
combinations
corDT_minN <- rbindlist(bpmapply(combinations$psiTypes, combinations$minN, 
                        SIMPLIFY = FALSE, FUN=function(type, minN){
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
                              
    xmat_sd <- rowSds(xmat)
    plotIdx <- rank(xmat_sd) >= length(xmat_sd) - topN
    xmat <- xmat[plotIdx,]
                              
    # round to limit effect of pseudocounts
    xmat <- round(xmat, digits=2)
    xmat[xmat == 0] <- 0.01
    xmat[xmat == 1] <- 0.99
                              
    xmat <- qlogis(xmat)
    xmat_rc <- xmat - rowMeans(xmat)
                              
    # set NAs if minRowN filter not met
    xmat_rc[snmat[plotIdx,] < minN] <- NA
                              
    # get raw sample correlations
    cormatRaw <- cor(xmat_rc, use="pairwise.complete.obs")
                          
    # get sample correlation after normalization
    pred_mu <- as.matrix(predictedMeans(fds, type=type)[
    expRowsMax & expRowsMedian,][plotIdx,])
    pred_mu <- round(pred_mu, digits=2)
    pred_mu[pred_mu == 0] <- 0.01
    pred_mu[pred_mu == 1] <- 0.99
    pred_mu <- qlogis(pred_mu)
    lpred_mu_rc <- pred_mu - rowMeans(pred_mu)
    norm_xmat_rc <- xmat_rc - lpred_mu_rc
    norm_xmat_rc[snmat[plotIdx,] < minN] <- NA # set NAs if minRowN filter not met
    cormatNorm <- cor(norm_xmat_rc, use="pairwise.complete.obs")
                              
    # calculate mean number of pairwise junctions 
    getPairwiseNrJunctionsFun <- function(i,j,data){ 
        sum(!(is.na(data[,i]) | is.na(data[,j])))
    }
    getPairwiseNrJunctions <- Vectorize(getPairwiseNrJunctionsFun, 
                                        vectorize.args=list("i","j"))
    pairwise_nrJunctions_raw <- outer(1:ncol(xmat_rc), 1:ncol(xmat_rc), 
                                    getPairwiseNrJunctions, data=xmat_rc)
    pairwise_nrJunctions_norm <- outer(1:ncol(norm_xmat_rc), 1:ncol(norm_xmat_rc), 
                                    getPairwiseNrJunctions, data=norm_xmat_rc)
    
    # create data table
    rbind(
          data.table(cor=cormatRaw[upper.tri(cormatRaw)],
             normalized="raw", dataset=dataset, method=method, type=type, 
             minN=minN, 
             nrJunctions=pairwise_nrJunctions_raw[upper.tri(pairwise_nrJunctions_raw)]),
          data.table(cor=cormatNorm[upper.tri(cormatNorm)],
             normalized="normalized", dataset=dataset, method=method, type=type, 
             minN=minN,  
             nrJunctions=pairwise_nrJunctions_norm[upper.tri(pairwise_nrJunctions_norm)]))
    }))

ans["SampleCors_minN"] <- list(corDT_minN)
ans


#' # Save numbers to RDS file
#'
saveRDS(ans, outFile)

