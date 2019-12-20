
#'---
#' title: Run FraseR Benchmark on GTEx data
#' author: Christian Mertes
#' output:
#'   html_document:
#'     pandoc_args: ["+RTS", "-K64m","-RTS"]
#'     toc: yes
#'     css: ../../js-css-html/fraser-markdown-report-style.css
#'---
#' render function, echo=FALSE, cache=FALSE
ren <- function(){
    # render it
    library(knitr)
    library(rmarkdown)
    source("src/r/config.R")

    project_webserver_dir <- "/s/public_webshare/project/fraser"
    rscript <- "src/r/benchmark/runGTExBenchmark.R"
    outdir <- file.path(project_webserver_dir, dirname(rscript))
    render(rscript, output_format = 'html_document', output_dir = outdir)
}

#+ load packages, echo=FALSE, cache=FALSE
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(ensemblVEP)
    library(vioplotx)
    library(data.table)
    library(Matrix)
})
opts_chunk$set(echo=FALSE, fig.width=12, fig.height=8, cache=TRUE)

#'
#' # Loading data
#'
#' Load GTEx FraseR object
#+ loading fraser data
fraserDir  <- "/s/project/fraser/analysis/datasets"
fraserName <- "gtex-fibroblast-transformed"
fraserDir  <- "/s/project/fraser/analysis/006-gtex"
fraserName <- "gtexfull"
fraserDir
fraserName
fds <- loadFraseRDataSet(fraserDir, fraserName)
fds
fdsres <- results(fds, zscoreCut=2, dPsiCut=0.01, fdrCut=1)
fdsres <- fdsres[abs(fdsres$zscore) >= 2]
fdsres <- fdsres[abs(fdsres$deltaPsi) >= 0.01]
fdsres <- fdsres[order(fdsres$p.adj)]

#'
#' Load GTEx sample table
#'
#+ load gtex sample table
gtexdt <- getSRAProjectTable(study_accession='SRP012682', mc.cores=15)

#'
#' Load GTEx variant object
#'
#+ load variants
GTEX_FILTERED_VARIANTS_RDS
gtexVcfObj <- readRDS(GTEX_FILTERED_VARIANTS_RDS)
gtexMatGT <- gtexVcfObj$gtexMatGT
dim(gtexMatGT)

#'
#' Generate sample map for VCF object
gesusuid <- sapply(strsplit(rownames(gtexMatGT), "-"),
    function(x) paste(x[1:2], collapse="-")
)
sampleMap <- gtexdt[submitted_subject_id %in% gesusuid & run %in% samples(fds),
    .(run, subject=submitted_subject_id)
]
setkey(sampleMap,subject)

#' only take samples which are in our dataset
gtexMatGT <- gtexMatGT[sampleMap[,subject],]

singltons <- apply(gtexMatGT, 2, function(x){
    freq <- table(x)
    if(any(!c("0", "1") %in% names(freq))){
        return(FALSE)
    }
    freq["0"] == length(x)-1 & freq["1"] == 1
})
sum(singltons)

gtSampleIdx <- apply(gtexMatGT[,singltons], 2, which.max)
subjectID <- rownames(gtexMatGT)[gtSampleIdx]
setkey(sampleMap,subject)
sampleMap <- sampleMap[subjectID]
sampleMap[,vcfName:=rownames(gtexMatGT)[gtSampleIdx]]
sampleMap

#'
#' Number of private variants per sample
barplot(sort(table(gtSampleIdx)), names.arg = 1:dim(gtexMatGT)[1],
    xlab="Sample rank", ylab="Frequency"
)

#'
#' Match variants to junctions/splice sites
#' +-2 bp around variant should the splice site be
#'
#+ prepare data for matching
fdsres
fdsres$hasVariant <- FALSE
fdsres$hasVariantAndSample <- FALSE
fdsres$subjectID <- as.character(NA)
fdsres$runID <- as.character(NA)

vcfgranges <- granges(gtexVcfGT[colnames(gtexMatGT[,singltons])])
subjGr <- resize(shift(vcfgranges, shift = -2), 4)
mcols(subjGr) <- sampleMap
subjGr

#+ match data
for(i in c(TRUE, FALSE)){
    ov <- findOverlaps(flank(granges(fdsres),1, start=i), subjGr)
    fdsres[from(ov)]$hasVariant <- TRUE
    fdsres[from(ov)]$subjectID <- sampleMap[to(ov), subject]
    fdsres[from(ov)]$runID <- sampleMap[to(ov), run]
    na2false(fdsres[from(ov)]$sampleID == sampleMap[to(ov), run])
}

gtSample <- sapply(gtexMatGT, which.max)

#' * Number of varians on significant junctions
sum(fdsres[fdsres$p.adj < 0.1]$hasVariant)
#' * Number of variants on significant junctions with correct sample
sum(fdsres[fdsres$p.adj < 0.1]$hasVariantAndSample)

#'
#' ## ROC
#' There is a bug in it. Need to check it. Looks like random
#+ generate roc
tmpres <- data.table(
    pvalue       = fdsres$pvalue,
    hasVar       = fdsres$hasVariant,
    hasSampleVar = fdsres$hasVariantAndSample,
    type         = fdsres$type
)
tmpres <- tmpres[order(pvalue)]
tmpres <- tmpres[pvalue < 0.01]
#tmpres <- tmpres[type != "psiSite"]
plot.roc(tmpres[order(pvalue), hasVar], tmpres[order(pvalue), pvalue])

par(mfrow=c(1,3))
sapply(c("psi3", "psi5", "psiSite"), function(psiType){
    tmp2plot <- tmpres[type==psiType]
    tmp2plot <- tmp2plot[order(pvalue)]
    tmp2plot[nrow(tmp2plot),hasVar:=TRUE]
    plot.roc(tmp2plot[order(pvalue), hasVar], tmp2plot[order(pvalue), pvalue], main=psiType)
})


#
# testedVars <- as.data.table(
#     mcols(fds, type="psiSite")[to(ov),c("spliceSiteID", "psiSite_tested")]
# )
# testedVars[,psi3_tested:=sdf]
# lapply(c("psi3", "psi5"), function(psiType, getID){
#     psiType="psi3"
#     getID=function(type) ifelse(type=="psi3", "startID", "endID")
#
#     spliceSiteID <- mcols(fds, type=psiType)[,getID(psiType)]
#     match <- testedVars[,spliceSiteID] %in% spliceSiteID
#
#     ans <- logical(nrow(testedVars))
#     ans[match] <-
#     $startID
# })
