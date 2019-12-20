#'---
#' title: Cluster counts
#' author: Christian Mertes
#' wb:
#'  input:
#'   - counts: '`sm getInputLeafcutterCounts`'
#'  output:
#'   - counts: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutter_processed/clustered_junctions/leafcutter_perind_numers.counts.gz"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/Leafcutter/{dataset}/counting.html"`'
#'  type: noindex
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic")
    parseWBHeader2("Scripts/Leafcutter/03_cluster_counts.R",
            wildcards=wildcards, rerun=TRUE)
}


#+ load config and setup, echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
rldsFile  <- snakemake@input$rlds
outFile   <- snakemake@output$counts
rldsFile  <- file.path(dirname(dirname(dirname(outFile))), "rlds_obj.RDS")

#+ read input files
rls <- readRDS(rldsFile)


#'
#' modify counts
clustering(rls)
leafcutterModifyCountTable(rls)
print_log("finished with clustering")

