#'---
#' title: Create SNP benchmark set for SCARED
#' author: Christian Mertes, Max Herzog, Daniel Bader
#' output:
#'   html_document:
#'     toc: yes
#'     toc_float: yes
#'---

#' get kerberos ticket
if(dir.exists("../gagneurlab_shared")){
    suppressPackageStartupMessages({
        source("../gagneurlab_shared/r/lab_orga/kerberos-api.R")})
    kerberos.init()
}

#+ load packages, echo=FALSE
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(ensemblVEP)
    library(data.table)
    library(Matrix)
    library(FraseR)
})

sourceFolder <- function(dir){
    DEV_NULL <- sapply(
        list.files(dir, pattern="\\.R$", full.names=TRUE),
        source
    )
}
# source extra functions
sourceFolder("src/r/functions")


#'
#'
#' # LOAD VARIANTS
#'
#'
#'
#' file names
#'
gtexFileRoot <- "/s/project/sra-download/files/dbGaP-11206"
wgsVcfRoot   <- file.path(gtexFileRoot, "phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1")
wesVcfRoot   <- file.path(gtexFileRoot, "")
wgsVcfAnnoFile <- file.path(wgsVcfRoot, "GTEx_Analysis_20150112_WholeGenomeSeq_VarSitesAnnot.vcf.gz")
wgsVcfGenoFile <- file.path(wgsVcfRoot, "GTEx_Analysis_20150112_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz")

# output
FILE_GTEX_FILTERED_VARIANTS_SCARED_RDS <- file.path('/s/project/scared/Data/GTEx/RDS/gtex_high_impact_variants.RDS')
# FILE_GTEX_FILTERED_VARIANTS_SCARED_RDS <- file.path('/s/project/scared/Data/GTEx/RDS/gtex_high_impact_variants_chr18.RDS')

#'
#' filtering parameters for the vcf parsing
#'
vcfFile      <- wgsVcfAnnoFile
scanVcfParam <- ScanVcfParam()
scanVcfParamChr01ALL <- ScanVcfParam(which=GRanges("1:1-536870912"))
scanVcfParamChr18ALL <- ScanVcfParam(which=GRanges("18:1-536870912"))

#'
# for debuging only
# scanVcfParam <- scanVcfParamChr18ALL # Comment this line to get the full dataset

#+ show defined file names, echo=TRUE
vcfFile
scanVcfParam


#'
#' load VCF, ca 4min (chr18)
#'
message("Loading annotated VCF")
system.time({
    raw_vcf_annotated <- readVcf(vcfFile, param=scanVcfParam)
})
message("...done")



#'
#'
#'  ## Filter the VCF file
#'
#'

#+ start filtering the variants, echo=TRUE
vcfDims <- list()
vcfDims[["all"]] <- length(raw_vcf_annotated)

# QC filter
vcfqc <- filterQC(raw_vcf_annotated, minAllel=1)
vcfDims[["QC"]] <- length(vcfqc)

# filter for HIGH effect mutations, see: https://www.ensembl.org/info/genome/variation/predicted_data.html
vcfqchigh <- filterCSQbyRegex(vcfqc, '\\|HIGH\\|')
vcfDims[["HIGH"]] <- length(vcfqchigh)

gc()
vcfDims




#'
#'
#'  ## Extract genotypes for the given sites from GT file
#'
#'
# def regions of filtered variants to load smaller file
scan_vcf_param_filtered <- ScanVcfParam(which=granges(vcfqchigh))

# duration ca 30sec (chr18)
message("Loading genotyped VCF")
system.time({
    raw_vcf_gt  <- readVcf(wgsVcfGenoFile, param=scan_vcf_param_filtered)
})
message("...done")
# restrict to filtered variants by index
vcfgt <- raw_vcf_gt[names(vcfqchigh),]

dim(vcfqchigh)
dim(vcfgt)



#'
#' ## Save filtered variants to file
#'
#+ save variants

message("Save VCF objects")
saveRDS(
    list(vcfGT=vcfgt, vcfAnno=vcfqchigh),
    FILE_GTEX_FILTERED_VARIANTS_SCARED_RDS
)


