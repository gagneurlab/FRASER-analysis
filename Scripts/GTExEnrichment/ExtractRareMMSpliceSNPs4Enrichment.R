#'---
#' title: Extract rare MMSplice variants for enrichment
#' author: Christian Mertes
#' wb:
#'   input:
#'     - variantTable:  "/s/project/gtex-processed/mmsplice-scripts/data/processed/mmsplice_rare_variants.csv"
#'   output:
#'     - variantTable: '`sm config["DATADIR"] + "/GTEx_variant_enrichment/rareMMSplice_filtered_VariantsTable.tsv.gz"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")

# Interactive mode (debuging)
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2(
        "Scripts/GTExEnrichment/ExtractRareMMSpliceSNPs4Enrichment.R",
        options=options, rerun=TRUE)
}

snptype <- "rareMMSplice"
vcfFile <- snakemake@input$variantTable
outFile <- snakemake@output$variantTable
tissues <- snakemake@config$EnrichmentTissues

MAF_LIMIT <- 0.05

vcfFile
outFile
snptype


###########################################
#'
#' # VCF parsing
#'
###########################################

#' Read in vcf file
#+ read vcf
variants <- fread(vcfFile)

# cleanup table
variants <- variants[,.(variantID=ID, subjectID=sample, MAF=maf,
        delta_logit_psi, Gene=gene_id, SYMBOL=gene_name,
        simple_conseq="MMSplice", IMPACT="MMSplicing")]

hist(variants[,.(MAF=max(MAF)), by=c('subjectID', 'variantID')][,log10(MAF)],
     breaks=30, main='MAF distribution over all variants of interest',
     xlab='log10(MAF)')
abline(v=log10(MAF_LIMIT), col='red')


#' MAF filter
table(variants[,MAF < MAF_LIMIT])
variants <- variants[MAF < MAF_LIMIT]


#' Number of variants per sample
hist(variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,N],
        main='Number of variants per sample', xlab='Number of variants')
variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,
        .(mean=round(mean(N), 2), sd=round(sd(N), 2))]

#' 
#' * Total number of variants/subject combinations considered in analysis
nrow(variants[,.N,by="variantID,subjectID"])

#'
#' * Delta Psi distribution in logit space
#' 
hist(variants[,delta_logit_psi], breaks=100)

#' Keep only the gene level and most severe variant annotation
#' (remove transcript level)
variants <- variants[order(subjectID, Gene, -abs(delta_logit_psi))]
dupVars <- duplicated(variants[,.(subjectID, Gene)])
table(dupVars)
variantsByGene <- variants[!dupVars]
sort(table(variantsByGene[,.N,by=c('variantID', 'simple_conseq')][,
        simple_conseq], useNA='always'))


#'
#' # Final var overview
#'

#' * Final number of variants per sample
hist(variantsByGene[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,N])
variantsByGene[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,
        .(mean=round(mean(N), 2), sd=round(sd(N), 2))]

#' Number of variants sample combinations
nrow(variantsByGene[,.N,by=c('variantID', 'subjectID')])

#' * Number of affected samples:
length(unique(variantsByGene$subjectID))

#' * Number of affected genes:
length(unique(variantsByGene$SYMBOL))
length(unique(variantsByGene$Gene))

#'
#' # Save variant table
#' 
fwrite(variantsByGene[,.(variantID, subjectID, MAF, Gene, 
        SYMBOL, simple_conseq, IMPACT, Consequence="MMSplice")], outFile)


#'
#' # Stats for paper
#'
#' * Number of tissues
length(tissues)

tissueDT <- rbindlist(lapply(tissues, function(x){
    dt <- fread(paste0("./Data/paperPipeline/annotations/", x, ".tsv"))
    dt[,.(tissue=x, subjectID=indivID)]
}))
tissueDT


#' ## Number of variants by individual
vdt <- variants[,.N,by="variantID,subjectID"][,.N,by=subjectID]
vdt[,.(mean=round(mean(N), 2), sd=round(sd(N), 2))]

#' ## Number of variants by tissue
merge(vdt, tissueDT)[,.(N=sum(N)), by=tissue][,
        .(mean=round(mean(N), 2), sd=round(sd(N), 2))]

#' ## Number of gene-sample combi by individual
vdt <- variantsByGene[,.N,by="variantID,subjectID"][,.N,by=subjectID]
vdt[,.(mean=round(mean(N), 2), sd=round(sd(N), 2))]

#' ## Number of gene-sample combi by tissue
merge(vdt, tissueDT)[,.(N=sum(N)), by=tissue][,
        .(mean=round(mean(N), 2), sd=round(sd(N), 2))]


