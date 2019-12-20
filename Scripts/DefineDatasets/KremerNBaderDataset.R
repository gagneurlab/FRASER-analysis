#'---
#' title: Define Kremer and Bader et al. datasets
#' author: Christian Mertes
#' wb:
#'  input:
#'    - sampleAnno:  '`sm config["KREMER_SAMPLE_ANNOTATION_FILE"]`'
#'    - fileMapping: "Data/filemapping/Kremer.tsv"
#'  output:
#'   - out:    '`sm config["DATADIR"] + "/annotations/Kremer.tsv"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/annotations/Kremer.html"`'
#'  type: noindex
#' output:
#'  html_document
#'---

if(FALSE){
    snakemake <- readRDS("tmp/snakemake.RDS")
}

#+ load main config, echo=FALSE
source("./src/r/config.R", echo=FALSE)

#+ input
outFile        <-snakemake@output$out
sampleAnnoFile <- snakemake@input$sampleAnno
mappingFile    <- snakemake@input$fileMapping

# dataset name
name <- gsub(".tsv$", "", basename(outFile))

#'
#' # Load and merge data
#'
anno <- fread(sampleAnnoFile)
mapping <- fread(mappingFile)
anno[,sampleID:=RNA_ID]
colData <- merge(anno, mapping)

#'
#' # Create Kremer and Bader et al dataset
#'
#' All samples used in the paper by Kremer-Bader et al.
#' Using the PAPER_MITOMAP object. In total we have 105 fibroblasts.
#'

#'
#' Dataset Name
#+ echo=TRUE
name


#'
#' Setting FIBROBLAST_ID to condition
#+ echo=FALSE
colData[,condition:=FIBROBLAST_ID]

#'
#' Add genome annotation
colData[,genome:="BSgenome.Hsapiens.UCSC.hg19"]

#+ echo=FALSE
finalTable <- colData

#'
#' ## Final sample table `r name`
#'

#+ savetable
setcolorder(finalTable, unique(c(
    "sampleID", "condition", "bamFile", colnames(finalTable))))

DT::datatable(finalTable, options=list(scrollX=TRUE))

dim(finalTable)
write_tsv(finalTable, file=outFile)

