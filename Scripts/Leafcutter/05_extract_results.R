#'---
#' title: Results of rareLeafcutter analysis
#' author: Christian Mertes
#' wb:
#'  input:
#'   - rlds: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/rlds_obj.RDS"`'
#'   - dsFiles: '`sm getInputLeafcutterDS`'
#'  output:
#'   - resultTable: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/results_{dataset}.tsv"`'
#'   - wBhtml: 'Output/html/Leafcutter/{dataset}_results.html'
#'  type: noindex
#'  threads: 10
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic")
    parseWBHeader2("Scripts/Leafcutter/05_extract_results.R",
            wildcards=wildcards, rerun=TRUE)
}


#+ load config and setup, echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
dataset   <- snakemake@wildcards$dataset
rldsFile  <- snakemake@input$rlds
outFile   <- snakemake@output$resultTable
FDR_LIMIT <- snakemake@config$FDR_LIMIT
bpWorkers <- min(bpworkers() * 3, snakemake@threads)


#+ read input files
rls <- readRDS(rldsFile)
BPPARAM <- MulticoreParam(bpWorkers)
register(BPPARAM)

#'
#' get annotation
clusterGeneMap <- clusterGeneMapping(rls)

#'
#' collect results
res <- collectResults(rls, clusterGeneMap=clusterGeneMap,
        plot=FALSE, FDR_CUTOFF=FDR_LIMIT, BPPARAM=BPPARAM)

#'
#' Write result table tsv
write_tsv(res, outFile)

