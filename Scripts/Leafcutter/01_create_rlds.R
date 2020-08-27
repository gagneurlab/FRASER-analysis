#'---
#' title: Create Rare leafcutter dataset
#' author: Christian Mertes
#' wb:
#'  input:
#'   - anno: '`sm config["DATADIR"] + "/annotations/{dataset}.tsv"`'
#'  output:
#'   - rlds:   '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/rlds_obj.RDS"`'
#'   - tsv:    '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/rlds_anno.tsv"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/Leafcutter/{dataset}/createRLDS.html"`'
#'  threads: 1
#'  type: noindex
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic")
    parseWBHeader2("Scripts/Leafcutter/01_create_rlds.R",
            wildcards=wildcards, rerun=TRUE)
}


#+ load config and setup, echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
annoFile  <- snakemake@input$anno
outFile   <- snakemake@output$rlds
outTsv    <- snakemake@output$tsv
wDir      <- dirname(outFile)


#+ read input files
anno <- fread(annoFile)

#' Swap names if on GTEx dataset
if(all(c("SMNABTCH", "indivID", "SAMPID", "sampleID") %in% colnames(anno))){
    # rename columns
    anno[,run:=sampleID]
    anno[,sampleID:=indivID]
    anno[,condition:=indivID]

    # check and halt if duplicated ids are present (replicats)
    if(any(duplicated(anno$sampleID))){
        stop("Please check the GTEx tissue for duplicated samples!")
    }
}
rls <- rareLeafcutterSet(anno, workingDir=wDir)


rls
DT::datatable(as.data.table(colData(rls)))


#+ save rlds object
saveRDS(rls, outFile)
write_tsv(anno, outTsv)

