#'---
#' title: Create GTEx SNP set for enrichment analysis
#' author: Christian Mertes
#' wb:
#'   threads: 10
#'   input:
#'     - subvcfannotbi: '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}Variants.vcf.gz.tbi"`'
#'     - subvcfgttbi:   '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}_gt.vcf.gz.tbi"`'
#'   output:
#'     - variantTable:  '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}VariantsTable.tsv"`'
#'     - wBhtml:        '`sm config["htmlOutputPath"] + "/GTEx_variant_enrichment/{snptype}_extraction.html"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---


# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(snptype="rareSplicing")
    parseWBHeader2("Scripts/GTExEnrichment/ExtractSNPsForEnrichment.R", wildcards=wildcards)
    slot(snakemake, name="wildcards", check=FALSE) <- wildcards
}

#'
#' Load config
source("./src/r/config.R")
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(ensemblVEP)
})
register(MulticoreParam(15))

#'
#' Read input
inputVCFAnno <- gsub('.tbi$', '', snakemake@input$subvcfannotbi)
inputVCFGT   <- gsub('.tbi$', '', snakemake@input$subvcfgttbi)
outputTSV    <- snakemake@output$variantTable
curType      <- snakemake@wildcards$snptype
MAF_LIMIT    <- 0.05

#'
#' # Load annotation vcf
#+ read vcf
vcf <- readVcf(inputVCFAnno)
dim(vcf)

#' Remove multi allele loci and not passed calls
vcf <- vcf[filt(vcf) == 'PASS']
alleleCount <- sapply(info(vcf)$AF, length)
vcf <- vcf[alleleCount == 1]
dim(vcf)

#' apply MAF limit
af <- sapply(info(vcf)$AF, max)
vcf <- vcf[af < MAF_LIMIT]
dim(vcf)
hist(log10(af))

#' parse CSQ field
#+ csq parsing
csq <- parseCSQToGRanges(vcf)
csqdt <- as.data.table(csq)
csqdt[,variantID:=names(csq)]
length(unique(csqdt$variantID))

#' apply MAF limit
csqdt[,MAF:=maf2number(csq$GMAF, naMAF=1/max(info(vcf)$AN))]
csqdt <- csqdt[MAF < MAF_LIMIT]
length(unique(csqdt$variantID))
hist(log10(csqdt$MAF))

#' remove unwanted transcripts and impact variants
#' https://www.ensembl.org/info/genome/variation/predicted_data.html
if(curType == "rareSplicing"){
    csqdt <- csqdt[grepl("splice_", Consequence)]
    csqdt <- csqdt[BIOTYPE %in% c('protein_coding', 'lincRNA',
            'retained_intron', 'nonsense_mediated_decay',
            'processed_transcript')]
} else if(curType == "rareModerateNHigh"){
    csqdt <- csqdt[IMPACT %in% c('MODERATE', 'HIGH')]
    csqdt <- csqdt[BIOTYPE %in% c('protein_coding', 'lincRNA')]
} else {
    stop("Please add your variant type to the script!")
}
length(unique(csqdt$variantID))

#'
#' # Extract genotypes for all variants of interest
#'
#+ get genotype
gr <- granges(vcf[unique(csqdt$variantID)])
vcfParam <- ScanVcfParam(info=c('AF'), geno=c('GT', 'AD'))
vcfgt <- readVcfParallel(inputVCFGT, gr, 50, nPerChunk=250, vcfParam=vcfParam)
vcfgt <- vcfgt[unique(csqdt$variantID)]
dim(vcfgt)


#:::::::::::::::::::::::::::::::::::::
#'
#' # Make tidy data
#'
#:::::::::::::::::::::::::::::::::::::
#'
#'
#' ## Table: variant, sample, genotype
#+ add genotype table
geno_dt <- as.data.table(geno(vcfgt)$GT, keep.rownames=T)
setnames(geno_dt, 'rn', 'variantID')
geno_dt <- melt(geno_dt, id.vars='variantID',
        variable.name='subjectID', value.name='GT')

#' Available genotypes
geno_dt[,.N,by=GT]

#' Get homozygous & heterozygous variants
geno_dt <- geno_dt[!GT %in% c('.', '0/0')]
hist(geno_dt[,.N,by='subjectID'][,N], main='Number of variants per sample')
hist(geno_dt[,.N,by='variantID'][,N], main='Number of sampels per variant')

#'
#' ## Table:  variant pos, ref, alt
#'
#' ALT: contains all annotated alternatives
#+ add range table
rowranges_dt <- data.table(
    seqnames=as.character(seqnames(vcfgt)),
    start=start(vcfgt),
    end=end(vcfgt),
    strand=as.character(strand(vcfgt)),
    variantID=names(vcfgt),
    QUAL=qual(vcfgt),
    FILTER=filt(vcfgt),
    REF=as.character(ref(vcfgt)))

rowranges_dt[,ALT:=unlist(mclapply(alt(vcfgt),
        function(x){ as.character(x[1]) }, mc.cores=10))]

# merge
setkey(geno_dt, NULL)
geno_dt <- merge(geno_dt, unique(rowranges_dt))


#'
#' ## Table: variant, sample, allelic depth
#'
#' * This table helps to select specific ALT alleles for further use
#' * Allelic depths for the ref and alt alleles in the order listed
#+ add allele depth
allele_depth_dt <- as.data.table(geno(vcfgt)$AD, keep.rownames=TRUE)
setnames(allele_depth_dt, 'rn', 'variantID')
allele_depth_dt <- melt(
    allele_depth_dt,
    id.vars = 'variantID',
    variable.name = 'subjectID',
    value.name = 'AllelicDepth'
)

# merge and reduce and simplify data
geno_dt <- merge(geno_dt, allele_depth_dt, by=c('variantID', 'subjectID'))
geno_dt[, AlleleDepthREF:=sapply(AllelicDepth, function(x) as.numeric(x[1]))]
geno_dt[, AlleleDepthALT:=sapply(AllelicDepth, function(x) as.numeric(x[2]))]
geno_dt[,AllelicDepth:=NULL]


#'
#' ## Table: variant info
#+ add info table
info_dt <- as.data.table(info(vcf))
info_dt[,variantID:=names(vcf)]
info_dt[,CSQ:=NULL]
info_dt[,AF:=unlist(AF)]
info_dt[,AC:=unlist(AC)]
info_dt[,MLEAC:=unlist(info_dt$MLEAC)]
info_dt[,MLEAF:=unlist(info_dt$MLEAF)]

# merge
geno_dt <- merge(geno_dt, info_dt, by='variantID')



#'
#' ## CSQ: ensembl VEP info
#'
#' CSQ is column in info(raw_vcf_annotated) with multiple values separated by "|"
#'

#' 1 variant can have multiple annotations from ENSEMBL VEP tool
csqdt[,numTranscripts:=.N, by=variantID]


#'
#' MERGE and keep all entries --> *multiple rows per variant and sample*
#'
geno_dt <- merge(geno_dt, csqdt, allow.cartesian=TRUE,
        by=c('variantID', 'seqnames', 'start', 'end', 'strand'))


#' Overview annotated info:
#'
colnames(geno_dt)
geno_dt

#'
#' # Save files
#+ save files
write.table(geno_dt, file=outputTSV, sep='\t', quote=FALSE, row.names=FALSE)

