#'---
#' title: Create variant table from GTEx for SCARED benchmark
#' author: Christian Mertes, Daniel Bader
#' output:
#'   html_document:
#'     toc: yes
#'     toc_float: yes
#'---

sourceFolder <- function(dir){
    DEV_NULL <- sapply(
        list.files(dir, pattern="\\.R$", full.names=TRUE),
        source
    )
}
# source extra functions
sourceFolder("src/r/functions")
suppressPackageStartupMessages({
    library(data.table)
})


#+ FILES
FILE_GTEX_FILTERED_VARIANTS_TIDY_TABLE <- file.path(
    # '/s/project/scared/Data/GTEx/gtex_high_impact_variants_chr18.tsv'
    '/s/project/scared/Data/GTEx/gtex_high_impact_variants.tsv'
)
FILE_GTEX_SAMPLE_MATCH_DNA2RNA_TABLE <- file.path(
    '/s/project/scared/Data/GTEx/gtex_high_impact_dna2rna_table.tsv'
)

#'
#'
#' # Get "subjects" with variants
#'
#'
#' * get gtex sample informations
gtexdt <- getSRAProjectTable(study_accession='SRP012682', mc.cores=15)
result_dt <- fread(FILE_GTEX_FILTERED_VARIANTS_TIDY_TABLE)

#' * merge subject id to variant table
result_dt <- merge(result_dt, gtexdt[,.(submitted_subject_id, submitted_sample_id)])

#' * get subjects
subjects <- sort(unique(result_dt$submitted_subject_id))


#'
#' ## Get RNA_ID
#'
#' * subset gtex samples for skin
#'
body_sites_skin=c(
    "Cells - Transformed fibroblasts",
    "Skin - Not Sun Exposed (Suprapubic)",
    "Skin - Sun Exposed (Lower leg)"
)
samples_skin <- gtexdt[
    submitted_subject_id %in% subjects & body_site %in% body_sites_skin,
    .(submitted_subject_id, submitted_sample_id, molecular_data_type, body_site, rna_run=run)
]

#' * get overview by site and data type
samples_skin[,.N, by=molecular_data_type]
samples_skin[,.N, by=body_site]

#'
#' * merge for DNA to RNA sample matching
#'
gtex_dna2rna_table <- merge(
    samples_skin,
    unique(result_dt[,.(submitted_subject_id, DNA_ID=submitted_sample_id)])
    )[order(submitted_subject_id)]
setnames(gtex_dna2rna_table, 'submitted_sample_id', 'RNA_ID')

#' * get also DNA run IDs for BAM filename
#'
gtex_dna2rna_table <- merge(gtex_dna2rna_table,
    gtexdt[,.(submitted_sample_id, dna_run=run)],
    by.x = "DNA_ID",
    by.y = "submitted_sample_id"
)


# SAVE
write.table(
    gtex_dna2rna_table,
    file=FILE_GTEX_SAMPLE_MATCH_DNA2RNA_TABLE,
    quote=F,
    sep='\t',
    row.names = F
)
