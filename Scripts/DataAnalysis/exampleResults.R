#'---
#' title: Example results with/out outliers
#' author: Ines Scheller
#' wb:
#'  params:
#'   - workers: 1
#'  input:
#'   - inFds: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{AE_impl}/predictedMeans_psiSite.h5"`'
#'  output:
#'   - outPdfs: '`sm expand(config["FIGDATADIR"] + "/resultsOverview/{{dataset}}_example{outlier}_{psiType}__{{AE_impl}}.pdf", psiType=config["psiTypes"], outlier=["Outliers", "NoOutliers"])`'
#'   - wBhtml: '`sm config["FIGDATADIR"] + "/resultsOverview/{dataset}__{AE_impl}.html"`'
#'  type: noindex
#' output:
#'  html_document
#'---

if(FALSE){
    rm(snakemake)
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic", AE_impl="PCA-BB-Decoder")
    parseWBHeader2("./Scripts/DataAnalysis/exampleResults.R", wildcards=wildcards)
    slot(snakemake, "wildcards", check=FALSE) <- wildcards
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input
dataset    <- basename(dirname(snakemake@input$inFds))
workingDir  <- file.path(CONFIG$DATADIR, "datasets")
FraseR_implementation <- snakemake@wildcards$AE_impl
# bpWorkers   <- min(bpworkers(), as.integer(snakemake@params[[1]]$workers))

#+ output
outPdfs <- snakemake@output$outPdfs
outPdfOutliers <- outPdfs[which(grepl(pattern="Outliers", x=outPdfs))]
outPdfNoOutliers <- outPdfs[which(grepl(pattern="NoOutliers", x=outPdfs))]

#'
#' # Get examples plots for the dataset
dataset

#+ load fds
fds <- loadFraseRDataSet(workingDir, dataset)

out <- lapply(psiTypes, function(type) {

    message(date(), ": ", type)

    #+ get pvalues by junction group
    padj   <- padjVals(fds, type, byGroup = TRUE)
    index  <- getSiteIndex(fds, type)
    groups <- unique(index)

    #+ find outlier examples
    outlierIndices <- which(padj < 1e-5)
    outliers <- sample(outlierIndices, min(100, length(outlierIndices)))

    #+ find sites with no outliers
    possibleNonOutlierSites <-which(rowSums(padj < 0.05) == 0)
    non_outliers <- sample(possibleNonOutlierSites, min(100, length(possibleNonOutlierSites)))

    #+ create pdf for outliers
    outlierPdf <- outPdfOutliers[which(grepl(pattern = paste0(strsplit(dataset, "__")[[1]][1], ".*", type), x = outPdfOutliers))]
    message(date(), ": Writing expamples to ", outlierPdf)
    pdf(
        outlierPdf,
        width = 10,
        height = 10,
        onefile = TRUE
    )
    for (ind in outliers) {
        arr.ind <- arrayInd(ind, dim(padj))
        sample <- arr.ind[, 2]
        junctions <- which(index == groups[arr.ind[, 1]])
        plot_list <- lapply(junctions, function(junction) {
            plotJunctionCounts(
                fds,
                type = type,
                site = junction,
                highlightSample = sample,
                title = paste0("Site: ", junction,
                               ", outlier sample: ", sample)
            )
        })
        g <- do.call("grid.arrange", plot_list)
        print(g)
    }
    dev.off()

    #+ create pdf for non-outliers
    nonOutlierPdf <- outPdfNoOutliers[which(grepl(pattern = paste0(strsplit(dataset, "__")[[1]][1], ".*", type), x =
                                                      outPdfNoOutliers))]
    message(date(), ": Writing non-outlier expamples to ", nonOutlierPdf)
    pdf(
        nonOutlierPdf,
        width = 10,
        height = 10,
        onefile = TRUE
    )
    for (site in non_outliers) {
        junctions <- which(index == groups[site])

        plot_list <- lapply(junctions, function(junction) {
            plotJunctionCounts(fds, type = type, site = junction)
        })

        g <- do.call("grid.arrange", plot_list)
        print(g)
    }
    dev.off()

    message(date(), ": Finished ", type)

})

