#'---
#' title: Results of FraseR analysis
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fdsin: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}__{method}/pajdBetaBinomial_psiSite.h5"`'
#'   - html:  '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_stat_calculation.html"`'
#'  output:
#'   - resultTable: '`sm config["DATADIR"] + "/processedData/results/{dataset}/{method}_results.tsv"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}/{method}_results.html"`'
#'  type: noindex
#'---

#+ load config and setup, echo=FALSE
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Minor_Salivary_Gland", method="PCA")
    parseWBHeader2("./Scripts/FraseR/07_extract_results_FraseR.R",
            wildcards=wildcards, rerun=TRUE)
}

#+ input
dataset    <- snakemake@wildcards$dataset
method     <- snakemake@wildcards$method
name       <- paste0(dataset, "__", method)
fdsFile    <- snakemake@input$fdsin
workingDir <- dirname(dirname(dirname(fdsFile)))

#+ load extra functions, echo=FALSE
R.utils::sourceDirectory("../gagneurlab_shared/r/disease")
R.utils::sourceDirectory("../gagneurlab_shared/r/go_enrichment")
MGSA_GO_FULL <- load_mgsaset_for_organism('human')
mimTable <- getFullOmimTable()
opts_chunk$set(fig.width=12, fig.height=8)


#'
#' # Load data
#'
#' Dataset:
#+ echo=TRUE
dataset
method
workingDir

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=name)
bpparam <- MulticoreParam(3, 3)
parallel(fds) <- bpparam

#'
#' ## Extract results
#'
fds <- annotateRanges(fds)
resgr <- results(fds, zScoreCutoff=NA, padjCutoff=0.9, deltaPsiCutoff=0.05)
res   <- as.data.table(resgr)
saveFraseRDataSet(fds)

#'
#' * Add features
#'     * number of samples per gene and variant
res[padjust<=0.1, numSamplesPerGene:=length(unique(sampleID)), by=hgncSymbol]
res[padjust<=0.1, numEventsPerGene:=.N, by="hgncSymbol,sampleID"]
res[padjust<=0.1, numSamplesPerJunc:=length(unique(sampleID)), by="seqnames,start,end"]

#'
#'     * MitoVIP genes
vip_genes <- get_vip_info_table()[,.(hgncSymbol=gene,
        isMitoVIP=ifelse(causal,"causal", "vip"))]
res <- merge(res, vip_genes, all.x=TRUE, "hgncSymbol")

#'
#'     * OMIM phenotypes
#'
res <- merge(res, mimTable, all.x=TRUE, by.x="hgncSymbol", by.y="SYMBOL")

#'
#'     * add colData at the end
res <- merge(res, as.data.table(colData(fds)), by="sampleID")

#'
#' # Results
#'
write_tsv(res, file=snakemake@output$resultTable)
file <- gsub(".html$", paste0("_", dataset, ".tsv"), snakemake@output$wBhtml)
write_tsv(res, file=file)

#'
#' The result table can also be downloaded with the link below.
#'
#+ echo=FALSE, results='asis'
cat(paste0("<a href='./", basename(file), "'>Download result table</a>"))

# get links
res[,genecards:=get_html_link(hgncSymbol, website="genecards", TRUE)]
res[,hgnc:=get_html_link(hgncSymbol, website="hgnc", TRUE)]
res[,omim:=get_html_link(GMIM, website="omim", TRUE)]
res[,entrez:=get_html_link(hgncSymbol, website="entrez", TRUE)]
res[,locus:=get_html_link(paste0(seqnames, ":", start, "-", end), website="locus", TRUE)]

# round numbers
res[,padjust:=signif(padjust, 3)]
res[,deltaPsi:=signif(deltaPsi, 2)]
res[,zScore:=signif(zScore, 2)]
res[,psiValue:=signif(psiValue, 2)]

# set correct order
setcolorder(res, unique(c("sampleID", "genecards", "padjust", "deltaPsi", "type",
        "numSamplesPerGene", "numEventsPerGene", "numSamplesPerJunc",
        "isMitoVIP", "omim", "PMIM", "PINH", "locus", "hgnc", colnames(res))))

#'
#' * Result table
DT::datatable(res[padjust <= 0.1 & abs(deltaPsi) > 0.3],
        options=list(scrollX=TRUE), escape=FALSE)

#'
#' * Sample table
DT::datatable(as.data.table(colData(fds)), options=list(scrollX=TRUE))

#'
#' * Sample correlation
for(type in psiTypes){
    plotCountCorHeatmap(fds, type=type, logit=TRUE, topN=30000, minMedian=5)
    if(method != "BB"){
        plotCountCorHeatmap(fds, type=type, logit=TRUE, topN=30000,
                minMedian=5, normalized=TRUE)
    }
}
