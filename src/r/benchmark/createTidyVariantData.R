#'---
#' title: Create SNP benchmark set for SCARED
#' author: Christian Mertes, Max Herzog, Daniel Bader
#' output:
#'   html_document:
#'     pandoc_args: ["+RTS", "-K64m","-RTS"]
#'     toc: yes
#'     toc_float: yes
#'---

#+ load packages, echo=FALSE
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(ensemblVEP)
    library(data.table)
})


#+ FILES
FILE_GTEX_FILTERED_VARIANTS_SCARED_RDS <- file.path('/s/project/scared/Data/GTEx/RDS/gtex_high_impact_variants_chr18.RDS')
FILE_GTEX_FILTERED_VARIANTS_TIDY_TABLE <- file.path('/s/project/scared/Data/GTEx/gtex_high_impact_variants_chr18.tsv')
# all
FILE_GTEX_FILTERED_VARIANTS_SCARED_RDS <- file.path('/s/project/scared/Data/GTEx/RDS/gtex_high_impact_variants.RDS')
FILE_GTEX_FILTERED_VARIANTS_TIDY_TABLE <- file.path('/s/project/scared/Data/GTEx/gtex_high_impact_variants.tsv')

#'
#' ## Save filtered variants to file
#'
#+ save variants

filtered_vcf <- readRDS(
    FILE_GTEX_FILTERED_VARIANTS_SCARED_RDS
)
vcfgt <- filtered_vcf$vcfGT
vcfqchigh <- filtered_vcf$vcfAnno


#:::::::::::::::::::::::::::::::::::::
#'
#' # Make tidy data
#'
#:::::::::::::::::::::::::::::::::::::
#'
#'
#' ## Table: variant, sample, genotype
#'
geno_dt <- as.data.table(geno(vcfgt)$GT, keep.rownames = T)
setnames(geno_dt, 'rn', 'variant_id')
geno_dt <- melt(geno_dt, id.vars = 'variant_id', variable.name = 'submitted_sample_id', value.name = 'GT')

#' Available genotypes
geno_dt[,.N,by=GT]

#' **get homozygous variants**
homo_high_var_dt <- geno_dt[GT=='1/1']



#'
#' ## Table: variant, sample, allelic depth
#'
#' * This table helps to select specific ALT alleles for further use
#' * Allelic depths for the ref and alt alleles in the order listed
#'
allele_depth_dt <- as.data.table(geno(vcfgt)$AD, keep.rownames = T)
setnames(allele_depth_dt, 'rn', 'variant_id')
allele_depth_dt <- melt(
    allele_depth_dt,
    id.vars = 'variant_id',
    variable.name = 'submitted_sample_id',
    value.name = 'AllelicDepth'
)

# merge and reduce data
homo_high_var_dt <- merge(homo_high_var_dt, allele_depth_dt)

#' simplify column
homo_high_var_dt[,
    length_AD:=lapply(AllelicDepth, function(x) as.numeric(length(x))),
    by=variant_id
]
homo_high_var_dt[,
    allelic_depth:=lapply(AllelicDepth, function(x) paste(as.numeric(x), collapse=',')),
    by=variant_id
]
homo_high_var_dt[,AllelicDepth:=NULL]

#' with homozygous variants we would expect all to be 2
table(homo_high_var_dt$length_AD)


#'
#' ## Table:  variant pos, ref, alt
#'
#' ALT: contains all annotated alternatives
#'
rowranges_dt <- as.data.table(rowRanges(vcfgt))
setnames(rowranges_dt, 'paramRangeID', 'variant_id')

#'
#' Column "ALT" is DNAStringSet.
#'
rowranges_dt[,alt_length:= lapply(ALT, function(x){
        length(as.character(x))
    }),
    by=variant_id
]
rowranges_dt[,alt_char:= lapply(ALT, function(x){
        paste(as.character(x), collapse = ',')
    }),
    by=variant_id
]
rowranges_dt[,ALT:=NULL]

# merge
setkey(homo_high_var_dt, NULL)
homo_high_var_dt <- merge(homo_high_var_dt, rowranges_dt)



#'
#' ## Table: variant info
#'
info_dt <- as.data.table(info(vcfqchigh))
info_dt[,variant_id:=rownames(vcfqchigh)]
info_dt[,CSQ:=NULL]

# merge
homo_high_var_dt <- merge(homo_high_var_dt, info_dt)




#'
#' ## CSQ: ensembl VEP info
#'

#' CSQ is column in info(raw_vcf_annotated) with multiple values separated by "|"
csqRanges <- parseCSQToGRanges(vcfqchigh)

#' Filter again for "HIGH" impact variants.
#' Necessary, because table before was selected by locus,
#' now we filter annotations per locus additionally
csqRanges <- csqRanges[csqRanges$IMPACT == "HIGH"]

#' 1 variant can have multiple annotations from ENSEMBL VEP tool
ensembl_vep_dt <- as.data.table(csqRanges)
ensembl_vep_dt[,variant_id:=names(csqRanges)]
ensembl_vep_dt[, num_vep_entries_per_variant:=.N, by=variant_id]

#' Overview annotated info:
#'
colnames(ensembl_vep_dt)

#' Info about multi annotations
#+ vep anno per variant
# barplot(
#     table(unique(ensembl_vep_dt[,.(variant_id, num_vep_entries_per_variant)])[,2]),
#     xlab='Num annotations',
#     ylab='Num variants'
# )
ensembl_vep_dt[num_vep_entries_per_variant==11, ]


#'
#' MERGE and keep all entries --> **multiple rows per variant and sample**
#'
result_dt <- merge(homo_high_var_dt, ensembl_vep_dt, allow.cartesian=T)

# SAVE
write.table(result_dt,
    file=FILE_GTEX_FILTERED_VARIANTS_TIDY_TABLE,
    sep='\t',
    quote=FALSE,
    row.names = FALSE
)


#'
#' # Overview results
#'

#' Number variants
length(unique(result_dt$variant_id))

#' Number affected samples
length(unique(result_dt$submitted_sample_id))

#' Number affected genes
length(unique(result_dt$SYMBOL))

