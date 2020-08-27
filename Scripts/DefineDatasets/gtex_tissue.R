#'---
#' title: Create GTEx tissue object
#' author: Christian Mertes
#' wb:
#'  input:
#'    - sampleAnno:  '`sm config["GTEX_SAMPLE_ANNOTATION_FILE"]`'
#'    - fileMapping: "Data/filemapping/GTEx.tsv"
#'    - annoPheno:   '`sm config["GTEX_PHENO_ANNOTATION_FILE"]`'
#'  output:
#'    - out:    '`sm config["DATADIR"] + "/annotations/{tissue,[^/]+}.tsv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/annotations/{tissue,[^/]+}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

#+ load main config, echo=FALSE
source("./src/r/config.R", echo=FALSE)

# for interactive run
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(tissue="Lung")
    parseWBHeader2("./Scripts/DefineDatasets/gtex_tissue.R",
            wildcards=wildcards, rerun=TRUE)
}


#+ input
outFile       <- snakemake@output$out
annoFile      <- snakemake@input$sampleAnno
annoPhenoFile <- snakemake@input$annoPheno
mappingFile   <- snakemake@input$fileMapping

#+ dataset name
name <- gsub(".tsv$", "", basename(outFile))

#'
#' # Load and merge Annotations
#'
name
anno      <- fread(annoFile)
annoPheno <- fread(annoPhenoFile)
mapping   <- fread(mappingFile)
sraDT     <- getSRAProjectTable()

#' Clean input variables
anno <- anno[!grepl('^K-562', SAMPID)] # remove k562 cell line
anno <- anno[!grepl('_rep[0-9]+$', SAMPID)] # remove replicates
anno[,SMTSD:=gsub("_$", "", gsub("[\\s-()]+", "_", SMTSD, perl=TRUE))]
anno[,indivID:=gsub("^([^-]+-[^-]+)-.*", "\\1", SAMPID, perl=TRUE)]

setnames(annoPheno, "SUBJID", "indivID")

sraDT[,body_site:=gsub("_$", "", gsub("[\\s-()]+", "_", body_site, perl=TRUE))]
sraDT2merge <- sraDT[body_site==name, .(SAMPID=submitted_sample_id, run)]
sraDT2merge <- merge(mapping, sraDT2merge, by.x="sampleID", by.y="run")

#' Subset and merge
colData <- merge(anno[SMTSD == name], annoPheno)
colData <- merge(colData, sraDT2merge, by.x="SAMPID", by.y="SAMPID")

#'
#' ## GTEx tissue: `r name`
#'
name

#'
#' ### Clean data
#'
#' * Remove unwanted samples from dataset
#'
colData <- colData[SMAFRZE == "USE ME"]
colData <- colData[SMRIN >= 5.7]

#'
#' Setting `submitted_subject_id` to condition
#+ setting-gtex-ids, echo=FALSE
colData[,condition:=sampleID]

#'
#' Add genome annotation
colData[,genome:="BSgenome.Hsapiens.1000genomes.hs37d5"]

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

