#'
#' Function to create test set and bam files for the FraseR package
#'
extractFraseRTestSetAndBamFiles <- function(){

    # generate sample and ranges of interest
    getTestsetGenomicRanges()

    # extract bam files
    system("bash -i ./src/bash/generate_example_data.sh")
}

#'
#'
#' helper function to create the test ranges
#' for the FraseR package
#'
#' Ranges of interest are: TIMMDC1 and MCOLN1
#'   * chr3:119217379-119243937 # TIMMDC1
#'   * chr19:7587512-7598895    # MCOLN1
#'
getTestsetGenomicRanges <- function(){
    library(FraseR)

    rangesOfInterest <- GRanges(
        seqnames = c("chr3", "chr19"),
        ranges = IRanges(
            start = c(119217379, 7587512),
            end = c(119243937, 7598895)
        )
    )

    # get sampledata table
    sampleData <- get_test_set_samples()

    # create FraseRSettings
    testSettings <- FraseRSettings(
        sampleData=sampleData,
        outputFolder=file.path("/s/project/fraser/testset"),
        bamParams=ScanBamParam(which=rangesOfInterest)
    )
    unlink(outputFolder(testSettings), recursive=TRUE)

    # get count data
    fds <- countRNAData(testSettings)
    fds <- calculateSitePSIValue(fds)

    # get only expressed splice sites
    min_expression <- 50
    expressed       <- rowMedians(as.matrix(as.data.table(
        assays(fds@nonSplicedReads)$rawCounts))) >= min_expression
    expressed_other <- rowMedians(as.matrix(as.data.table(
        assays(fds@nonSplicedReads)$rawOtherCounts_sitePSI))) >= min_expression

    # show the counts of expressed junctions
    table(expressed_other | expressed)

    # extract splice sites for expressed ranges
    gr  <- rowRanges(fds@nonSplicedReads)
    grr <- sort(reduce(gr[expressed_other | expressed]))

    # remove non overlapping regions
    grr <- subsetByOverlaps(grr, rangesOfInterest)

    # show reduced grange object
    grr

    # print out the ranges
    cat(paste0(seqnames(grr), ":", start(grr), "-", end(grr)),
        sep="\n", file="src/bash/rangesOfInterest.txt"
    )
}

#'
#' create the sampledata test set from the mitomap for FraseR
#'
#' Samples of interest are: TIMMDC1 and MCOLN1 CLPP + all NHDF samples
#'
get_test_set_samples <- function(){
    load_mitomap()

    # retrive data from mitomap
    sampleData <- MITOMAP_DATATABLE[
            !is.na(RNA_ID)
            & TISSUE == "FIBROBLAST"
            & KNOWN_MUTATION %in% c("TIMMDC1", "MCOLN1", "CLPP", "NHDF")
            & is.na(TRANSDUCED_GENE)
            & GROWTH_MEDIUM == "GLU"
    ]

    # clean table and add bamfile
    sampleData <- sampleData[,.(
        sampleID=RNA_ID,
        bamFile=sapply(RNA_ID,
                get_helmholtz_file, source="RNA", type="BAM.STAR"
        ),
        gene=KNOWN_MUTATION
    )]

    # add group
    sampleData[,group:=as.integer(factor(
            gene, levels = c("TIMMDC1", "MCOLN1", "CLPP")
    ))]

    # sort it and write it out
    sampleData <- sampleData[order(group)]
    write_tsv(sampleData, file="src/bash/samplesOfInterest.txt")

    # change to bamFile object
    sampleData <- sampleData[, bamFile:=unlist(BamFileList(bamFile))]

    # return it
    return(sampleData)
}
