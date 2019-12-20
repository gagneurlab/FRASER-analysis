#'---
#' title: Full GTEx Rare variant enrichtment
#' author: Christian Mertes
#' wb:
#'   threads: 2
#'   input:
#'     - rds:         '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}.RDS",  dataset=config["EnrichmentTissues"], snptype=["rareSplicing", "rareMMSplice"], deltaPsi=["0.0", "0.1"])`'
#'     - html: '`sm expand(config["htmlOutputPath"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}.html", dataset=config["EnrichmentTissues"], snptype=["rareSplicing", "rareMMSplice"], deltaPsi=["0.0", "0.1"])`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    parseWBHeader2("Scripts/GTExEnrichment/GTExFullSNPEnrichment.R", rerun=TRUE)
}


#' 
#' # TODO make links to htmls
#' 
1 + 2

