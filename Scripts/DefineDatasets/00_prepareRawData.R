#'---
#' title: Get/Create input files FraseR
#' author: Christian Mertes
#' wb:
#'  output:
#'    - KREMER_SAMPLE_ANNOTATION_FILE: '`sm config["KREMER_SAMPLE_ANNOTATION_FILE"]`'
#'    - GTEX_SAMPLE_ANNOTATION_FILE:   '`sm config["GTEX_SAMPLE_ANNOTATION_FILE"]`'
#'    - GTEX_PHENO_ANNOTATION_FILE:    '`sm config["GTEX_PHENO_ANNOTATION_FILE"]`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---
#'

#+ load config
source('src/r/config.R')

#+ input
gtexAnnoFile <- snakemake@output$GTEX_SAMPLE_ANNOTATION_FILE
GTEX_VERSION <- basename(dirname(gtexAnnoFile))

#'
#' Download raw data files for the full analysis
#'

GTEX_URLS <- list(
    V7 = list(
        anno  = "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt",
        pheno = "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
    ),
    V6P = list(
        anno  = "https://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
        pheno = "https://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
    ),
    V6 = list(
        anno = "https://storage.googleapis.com/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
        pheno = "https://storage.googleapis.com/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
    ),
    inhouse = list(
        anno = "/data/nasif12/home_if12/mertes/projects/OUTRIDER-analysis-final/Data/GTEx_SRA_anno_file.tsv",
        pheno = "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
    )
)

URLS <- c(
    GTEX_SAMPLE_ANNOTATION_FILE   = GTEX_URLS[[GTEX_VERSION]][['anno']],
    GTEX_PHENO_ANNOTATION_FILE    = GTEX_URLS[[GTEX_VERSION]][['pheno']],
    KREMER_SAMPLE_ANNOTATION_FILE = 'https://media.nature.com/original/nature-assets/ncomms/2017/170612/ncomms15824/extref/ncomms15824-s8.txt'
)

FILE_NAMES <- snakemake@output[names(URLS)]

for(i in names(URLS)){
    file2save <- FILE_NAMES[[i]]
    message('Downloading file:', i)
    message('From: ', URLS[i])
    message('To: ', file2save)
    download.file(URLS[i], destfile=file2save, method='wget', extra='-c')
}
