#'---
#' title: Run Differential Expression Analysis
#' author: Christian Mertes
#' wb:
#'  input:
#'   - counts: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutter_processed/clustered_junctions/leafcutter_perind_numers.counts.gz"`'
#'  output:
#'   - singleRes: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/ds_testing/{condition}_versus_rest/leafcutter_ds_cluster_significance.txt"`'
#'   - wBhtml: 'Output/html/Leafcutter/{dataset}/single_ds/{condition}_results.html'
#'  threads: 10
#'  type: noindex
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic", condition="GTEX-ZVP2")
    parseWBHeader2("Scripts/Leafcutter/04_run_single_ds_leafcutter.R",
            wildcards=wildcards, rerun=TRUE)
}


#+ load config and setup, echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
condition <- snakemake@wildcards$condition
rldsFile  <- snakemake@input$rlds
outFile   <- snakemake@output$singleRes
threads   <- snakemake@threads
rldsFile  <- file.path(dirname(dirname(dirname(outFile))), "rlds_obj.RDS")

#+ read input files
rls <- readRDS(rldsFile)


#'
#' run differential test per sample
dsTesting(condition, config=rls, threads=threads)
