#'---
#' title: Counting with Leafcutter
#' author: Christian Mertes
#' wb:
#'  input:
#'   - bamFile: '`sm getInputLeafcutterBams`'
#'  output:
#'   - counts: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutter_processed/junction_files/{sampleID}.junc"`'
#'   - wBhtml: 'Output/html/Leafcutter/{dataset}/single_sample/counting_{sampleID}.html'
#'  threads: 3
#'  type: noindex
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic", sampleID="GTEX-11DXX")
    snakemake <- parseWBHeader2("Scripts/Leafcutter/02_counting_for_leafcutter.R", wildcards=wildcards, debug=TRUE)
    slot(snakemake, "wildcards", check=FALSE) <- wildcards
}


#+ load config and setup, echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
sampleID  <- snakemake@wildcards$sampleID
bamFile   <- snakemake@input$bamFile
outFile   <- snakemake@output$counts
rldsFile  <- file.path(dirname(dirname(dirname(outFile))), "rlds_obj.RDS")

#+ read input files
rls <- readRDS(rldsFile)


#'
#' count bam file and create junc files convert all bamfiles to co
convertBam2Junc(sampleID, config=rls)


