#'---
#' title: GTEx Rare variant enrichtment Outlier extraction (per tissue)
#' author: Christian Mertes
#' wb:
#'   threads: 9
#'   input:
#'     - leafcutter:    '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/results_{dataset}.tsv"`'
#'     - leafcutterMD:  '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/results_{dataset}.tsv"`'
#'     - leafcutterMD_html: '`sm config["htmlOutputPath"] + "/Leafcutter/{dataset}_results_MD.html"`'
#'     - spot:          '`sm config["DATADIR"] + "/processedData/spot/{dataset}/spot__fullResults.tsv"`'
#'     - annotation:    '`sm config["DATADIR"] + "/annotations/{dataset}.tsv"`'
#'     - datasets_all:  '`sm expand(config["DATADIR"] + "/datasets/savedObjects/{{dataset}}__{method}/padjBetaBinomial_psiSite.h5", method=config["GTExMethods"])`'
#'     - results_all:   '`sm expand(config["DATADIR"] + "/processedData/results/{{dataset}}/{method}_results.tsv", method=config["GTExMethods"])`'
#'   output:
#'     - table:         '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__outlierStatus__{deltaPsi}.tsv.gz"`'
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_variant_enrichment/{dataset}__outlierStatus__{deltaPsi}.html"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")
devtools::load_all("../rare-disease-leafcutter/")

# Interactive mode (debuging)
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Uterus", deltaPsi="0_0")
    options   <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/GTExEnrichment/03_ExtractOutlierStatusPerTissue.R",
            wildcards=wildcards, options=options, rerun=TRUE)
}

curAEVersion <- snakemake@config$AE_IMPLEMENTATION
tissue     <- snakemake@wildcards$dataset
deltaPsi   <- as.numeric(gsub("_", ".", snakemake@wildcards$deltaPsi))
rlc_file   <- snakemake@input$leafcutter
lcmd_file  <- snakemake@input$leafcutterMD
spot_file  <- snakemake@input$spot
fds_files  <- snakemake@input$datasets_all
res_files  <- snakemake@input$results_all
anno_file  <- snakemake@input$annotation
threads    <- snakemake@threads

outtable   <- snakemake@output$table

MIN_DELTA_PSI <- deltaPsi
MIN_COVERAGE  <- 10
BPPARAM       <- MulticoreParam(min(4, threads))
register(BPPARAM)

#' 
#' # User input data
#' 
anno_file
rlc_file
lcmd_file
spot_file
res_files
outtable
BPPARAM

###########################################
#'
#' # Read in outlier calls
#'
###########################################

###########################################
#'
#' ## Read Leafcutter calls
#' 
#+ leafcutter readin
rlds <- readRDS(file.path(dirname(rlc_file), "rlds_obj.RDS"))
clusterGeneMap <- clusterGeneMapping(rlds)
tmp_dt <- collectResults(rlds, clusterGeneMap=clusterGeneMap,
        plot=FALSE, FDR_CUTOFF=2, BPPARAM=BPPARAM, save=FALSE)
tmp_dt <- tmp_dt[!is.na(p)]

#' collect info:
#' subjectID, geneID, p_leafcutter, fdr_leafcutter
#'   --> pvalues are holm corrected within genes
enrich_rlc <- rbindlist(bplapply(tmp_dt[,unique(sampleIDs)], 
        data=tmp_dt[!is.na(genes)], BPPARAM=BPPARAM, 
        FUN=function(x, data){
            data <- data[sampleIDs == x]
            data <- data[,.(subjectID=sampleIDs, geneID=genes, p=gene_p)]
            data <- data[,.(tissue=tissue, Leafcutter_p=min(p.adjust(p, method="holm"))), by="subjectID,geneID"]
            data
}))
enrich_rlc[,Leafcutter_fdr:=p.adjust(Leafcutter_p, method="BY"), by="subjectID"]
tibble(enrich_rlc)

##########################################
#'
#' ## Read LeafcutterMD calls
#' 
#+ leafcutterMD readin
tmp_dt <- fread(lcmd_file)
tmp_dt[abs(maxEffect) < MIN_DELTA_PSI, c("pvalue_gene", "padj"):=list(1, 1)]
enrich_lcmd <- tmp_dt[,.(subjectID=sample, geneID, tissue=tissue,
        LeafcutterMD_p=pvalue_gene, LeafcutterMD_fdr=padj, Leafcutter_effect=maxEffect)]
tibble(enrich_lcmd)

###########################################
#'
#' ## Read SPOT calls
#' 
#+ spot readin
tmp_dt <- fread(spot_file)
enrich_spot <- tmp_dt[,.(subjectID=SAMPLE_ID, geneID=GENE_ID, tissue=tissue,
        SPOT_p=gene_p, SPOT_fdr=gene_fdr)]
# reduce it to uniq rows per subject/gene
enrich_spot <- enrich_spot[!is.na(SPOT_fdr)]
tibble(enrich_spot)


###########################################
#'
#' ## Read FRASER calls
#'
#+ fds readin
loadFds_obj <- FALSE
ncpus <- 1
enrich_fds_obj <- mclapply(fds_files, load_fraser_enrichment_data,
        mc.cores=threads, internCPUs=ncpus, minCoverage=MIN_COVERAGE,
        minDeltaPsi=MIN_DELTA_PSI, debug=loadFds_obj)
enrich_fds <- lapply(enrich_fds_obj, "[[", "enrich_obj")
fds_ls     <- lapply(enrich_fds_obj, "[[", "fds")
res_ls     <- lapply(enrich_fds_obj, "[[", "res")
names(enrich_fds) <- basename(dirname(fds_files))
names(fds_ls) <- basename(dirname(fds_files))
names(res_ls) <- basename(dirname(fds_files))


#'
#' ### use subject ID aka individual identifier instead of run id
#' 
#+ anno readin
anno <- fread(anno_file)
enrich_fds <- lapply(enrich_fds, anno=anno, function(x, anno){
    setnames(x, "subjectID", "run")
    x <- merge(x, anno[,.(run=sampleID, indivID)], by="run")
    setnames(x, "indivID", "subjectID")
    x[,run:=NULL]
    x
})
enrich_fds <- lapply(enrich_fds, function(x) { 
    x[,tissue:=tissue]
    setcolorder(x, unique(c("subjectID", "geneID", "tissue", colnames(x))))
    x
})


#' 
#' Merge all data sets
#' 
list2merge <- c(list(rlc=enrich_rlc, lcmd=enrich_lcmd, spot=enrich_spot), enrich_fds)
merged_dt <- Reduce(x=list2merge, f=function(x, y){ 
        merge(x, y, by=c('subjectID', 'geneID', "tissue"), all=TRUE) })


#' 
#' write out table
#' 
fwrite(merged_dt, outtable)

