#:::::::::::::::::::::::::::::::::::::
#'
#'  FUNCTIONS
#'
#:::::::::::::::::::::::::::::::::::::
#+ define filter functions

#'
#' checks for a "clean" single nucleotide variant
#'
filterForSNP <- function(vcf){
    vcf[VariantAnnotation:::.testForSNV(ref(vcf), alt(vcf))]
}

#'
#' some default QC filter steps
#'
filterQC <- function(vcf, mq=30, qual=90, minAllel=296,
                     filter=c(".", "PASS", "VQSRTrancheSNP99.90to100.00",
                              "VQSRTrancheSNP99.80to99.90"
                     )){

    filterPass   <- rowRanges(vcf)$FILTER %in% filter
    mqPass       <- info(vcf)$MQ >= mq
    qualPass     <- qual(vcf) >= qual

    if("AN" %in% colnames(info(vcf))){
        minAllelPass <- info(vcf)$AN >= minAllel
    } else if("ExAC_AN" %in% colnames(info(vcf))){
        minAllelPass <- info(vcf)$ExAC_AN >= minAllel
        minAllelPass <- sapply(minAllelPass, isTRUE)
    } else {
        minAllelPass <- TRUE
    }

    vcf[na2false(filterPass & mqPass & qualPass & minAllelPass)]
}

#'
#' filter for maximal ALLEL count
#'   * assumes SNPs only
#'
filterMaxAC <- function(vcf, maxAC=4){
    vcf[sapply(info(vcf)$AC, sum, na.rm=TRUE) <= maxAC]
}

#'
#' filters for a splice affecting variant
#'
filterSplicing <- function(vcf, csq=NULL){
    if(is.null(csq)){
        csq <- parseCSQToGRanges(vcf)
    }

    vcf[unique(names(csq[grepl("splice", csq$Consequence)]))]
}

#'
#' CSQ quick filter by regex
#'
filterCSQbyRegex <- function(vcf, filterRegex='splice'){

    # match regex
    elt  <- elementNROWS(info(vcf)[["CSQ"]])
    ulst <- unlist(info(vcf)[["CSQ"]])
    regexMatched <- grepl(filterRegex, ulst)

    # assign match to correct rowID
    rowIDs <- 1:length(vcf)
    matchedRowIDs <- rowIDs[rep(rowIDs, elt)]

    # return only the regex matched rows
    vcf[unique(matchedRowIDs[regexMatched])]
}


#'
#' stats on the snp allel distribution
#'
plotSnpDistribution <- function(gtMatrix){

    gtPerSampleStat <- sapply(1:nrow(gtMatrix), FUN=function(idx) c(
        "Reference"    = sum(gtMatrix[idx,]==0),
        "Heterozygous" = sum(gtMatrix[idx,]==1),
        "Homozygous"   = sum(gtMatrix[idx,]==2),
        "Others"       = sum(gtMatrix[idx,]>2)
    ))

    gtPerLocus <- sapply(1:dim(gtMatrix)[2], FUN=function(idx) c(
        "Reference"    = sum(gtMatrix[,idx]==0),
        "Heterozygous" = sum(gtMatrix[,idx]==1),
        "Homozygous"   = sum(gtMatrix[,idx]==2),
        "Others"       = sum(gtMatrix[,idx]>2)
    ))


    colsOfInt <- 2:3
    gtPerLocusVio      <- lapply(colsOfInt, function(x) gtPerLocus[x,])
    gtPerSampleStatVio <- lapply(colsOfInt, function(x) gtPerSampleStat[x,])
    names(gtPerLocusVio) <- rownames(gtPerLocus)[colsOfInt]
    names(gtPerSampleStatVio) <- rownames(gtPerSampleStat)[colsOfInt]

    par(mfrow=c(2,2))
    boxplot(t(gtPerSampleStat)[,2:4], main="Number of Variants per sample",
            xlab="Varianttype", ylab="#SNPs"
    )
    boxplot(t(gtPerLocus)[,2:4], main="Number of samples per locus",
            xlab="Varianttype", ylab="#samples"
    )
    vioplotx(gtPerSampleStatVio, main="Number of Variants per sample",
             xlab="Varianttype", ylab="#SNPs",
             names=names(gtPerSampleStatVio)
    )
    vioplotx(gtPerLocusVio, main="Number of samples per locus",
             xlab="Varianttype", ylab="#samples",
             names=names(gtPerLocusVio)
    )

    return(list(
        gtPerSampleStat = gtPerSampleStat,
        gtPerLocus = gtPerLocus
    ))
}

readBigVcf <- function(vcfFile, filterRange, size=200,
                    workers=30, tasks=workers*10, ...){

    filterRange <- reduce(granges(filterRange), drop.empty.ranges=TRUE)
    optSize <- min(size, ceiling((length(filterRange)/workers + 1)/2))
    chunks <- chunk(seq_along(filterRange), optSize)
    bpparam <- MulticoreParam(workers, tasks, progressbar = TRUE)

    res <- bplapply(chunks, filter=filterRange, file=vcfFile, BPPARAM=bpparam,
        function(i, filter, file, ...){
            vcfparam <- ScanVcfParam(which=filter[i], ...)
            readVcf(TabixFile(file), param=vcfparam)
        })

    do.call(rbind, res)
}
