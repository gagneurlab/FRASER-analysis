#'---
#' title: Create SNP benchmark set for FraseR
#' author: Christian Mertes
#' output:
#'   html_document:
#'     pandoc_args: ["+RTS", "-K64m","-RTS"]
#'     toc: yes
#'     css: ../../../../gagneurlab_shared/r/knitr_plugins/add_content_table.css
#'---
#' render function, echo=FALSE
ren <- function(){
    # render it
    library(knitr)
    library(rmarkdown)
    project_webserver_dir <- "/s/public_webshare/project/fraser"
    rscript <- "src/r/benchmark/getSnpForBenchmark.R"
    outdir <- file.path(project_webserver_dir, dirname(rscript))
    render(rscript, output_format = 'html_document', output_dir = outdir)
}

#+ load packages, echo=FALSE
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(ensemblVEP)
    library(vioplotx)
    library(data.table)
    library(Matrix)
    library(FraseR)
})
source("src/r/config.R")
opts_chunk$set(echo=FALSE, fig.width=12, fig.height=8, cache=TRUE)

#'
#' file names
#'
gtexFileRoot <- "/s/project/sra-download/files/dbGaP-11206"
wgsVcfRoot   <- file.path(gtexFileRoot, "phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1")
# old files (148 samples)
wgsVcfAnnoFile <- file.path(wgsVcfRoot, "GTEx_Analysis_20150112_WholeGenomeSeq_VarSitesAnnot.vcf.gz")
wgsVcfGenoFile <- file.path(wgsVcfRoot, "GTEx_Analysis_20150112_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz")
# new files (635 samples)
wgsVcfAnnoFile <- "/s/project/fraser/benchmark/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Indiv_GATK_HaplotypeCaller.siteOnly.HIGH.vep.vcf.gz"
wgsVcfGenoFile <- "/s/project/sra-download/files/dbGaP-11206/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz"


#'
#' filtering parameters for the vcf parsing
#'
scanVcfParam <- ScanVcfParam()
scanVcfParamChr18CSQ <- ScanVcfParam(fixed=NA, info="CSQ", geno=NA, which=GRanges("18:1-536870912"))
scanVcfParamChr18ALL <- ScanVcfParam(which=GRanges("18:1-536870912"))


#'
#' set the parameters to be used in this analysis
#'
vcfFile      <- wgsVcfAnnoFile
scanVcfParam <- scanVcfParam
# scanVcfParam <- scanVcfParamChr18ALL

#+ show defined file names, echo=TRUE
vcfFile
scanVcfParam


#'
#' load it
#'
if(!"vcf" %in% ls()){
    system.time({ vcf <- readVcf(vcfFile, param=scanVcfParam) })
}


#:::::::::::::::::::::::::::::::::::::
#'
#'  Filter the VCF file
#'
#:::::::::::::::::::::::::::::::::::::
#+ start filtering the variants, echo=TRUE
vcfDims <- list()
# default filtering (QC, max allele count)
vcfDims[["all"]] <- length(vcf)
vcfFilt <- filterQC(vcf, minAllel=1)
vcfDims[["QC"]] <- length(vcfFilt)

# keep this number high so we get still private events in a subpopulation
vcfFilt <- filterMaxAC(vcfFilt, 50)
vcfDims[["maxAC"]] <- length(vcfFilt)
vcfFilt <- vcfFilt[sapply(info(vcfFilt)$AC, length) == 1]
vcfDims[["diploidOnly"]] <- length(vcfFilt)

# specific splicing selection
vcfFilt <- filterCSQbyRegex(vcfFilt, 'splice')
vcfDims[["spliceregion"]] <- length(vcfFilt)
vcfFilt <- filterCSQbyRegex(vcfFilt, 'splice_[ad]')
vcfDims[["splicesite"]] <- length(vcfFilt)

vcfFiltsnp <- filterForSNP(vcfFilt)
vcfDims[["snpOnly"]] <- length(vcfFiltsnp)

gc()


vcfDims

#:::::::::::::::::::::::::::::::::::::
#'
#'  Extract genotypes for the given sites
#'
#:::::::::::::::::::::::::::::::::::::

# vcfOnlyGeno <- readBigVcf(wgsVcfGenoFile, vcfFiltsnp, worker=40)
vcfgt <- readBigVcf(wgsVcfGenoFile, vcfFilt, worker=40)

# remove entries which overlap the region, but are of no interest
vcfgt <- vcfgt[names(vcfgt) %in% names(vcfFilt)]

# get GT matrix
gtMatrix <- as(genotypeToSnpMatrix(vcfgt)$genotype, "numeric")

# remove non-diploid entries
gtMatrix <- gtMatrix[,!apply(gtMatrix, 2, function(x) any(is.na(x)))]
dim(gtMatrix)

gtStats <- plotSnpDistribution(gtMatrix)
attach(gtStats)

#' variants to be used in the benchmark
table(data.table(t(gtPerLocus)[,2:3])[Homozygous <= 5 & Heterozygous <= 5])

#' Save filtered variants to file
#+ save variants
GTEX_FILTERED_VARIANTS_RDS
saveRDS(list(gtexVcfGT=vcfgt, gtexVcfAnnoSS=vcfFilt, gtexMatGT=gtMatrix),
        GTEX_FILTERED_VARIANTS_RDS)

