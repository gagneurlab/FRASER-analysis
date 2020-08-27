#'---
#' title: Run LeafcutterMD Analysis
#' author: Christian Mertes
#' wb:
#'  input:
#'   - counts: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutter_processed/clustered_junctions/leafcutter_perind_numers.counts.gz"`'
#'  output:
#'   - jpvals: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/leafcutter_MD_pVals.txt"`'
#'   - cpvals: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/leafcutter_MD_clusterPvals.txt"`'
#'   - effect: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/leafcutter_MD_effSize.txt"`'
#'   - wBhtml: 'Output/html/Leafcutter/{dataset}/leafcutter_md_testing.html'
#'  threads: 10
#'  type: noindex
#'---

if(FALSE){
    source(file=".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic")
    parseWBHeader2("Scripts/LeafcutterMD/01_run_leafcutterMD.R",
            wildcards=wildcards, rerun=TRUE)
}


#+ load config and setup, echo=FALSE
source(file="./src/r/config.R")
load_all("../rare-disease-leafcutter/")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
rldsFile  <- snakemake@input$rlds
outFile   <- snakemake@output$jpvals
threads   <- snakemake@threads
rldsFile  <- file.path(dirname(dirname(outFile)), "rlds_obj.RDS")

#+ read input files
rls <- readRDS(rldsFile)


#'
#' run differential test per sample
leafcutterMDTesting(rls, threads)


