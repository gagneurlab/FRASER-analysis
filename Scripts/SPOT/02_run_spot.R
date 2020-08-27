#'---
#' title: Run SPOT
#' author: Ines Scheller
#' wb:
#'  input:
#'   - intron_clusters: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutter_processed/clustered_junctions/leafcutter_perind_numers.counts.gz"`'
#'  output:
#'   - mahalanobis_dist: '`sm config["DATADIR"] + "/processedData/spot/{dataset}/spot__md.txt.gz"`'
#'   - pvalues:  '`sm config["DATADIR"] + "/processedData/spot/{dataset}/spot__emperical_pvalue.txt.gz"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/SPOT/{dataset}/run_spot.html"`'
#'  threads: 30
#'  type: noindex
#'---

#+ load config and setup, echo=FALSE
source("./src/r/config.R")

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic")
    parseWBHeader2("Scripts/SPOT/02_run_spot.R",
            wildcards=wildcards, rerun=TRUE)
}


#+ input
intron_clusters  <- snakemake@input$intron_clusters
mahalanobis_dist <- snakemake@output$mahalanobis_dist
pvalues          <- snakemake@output$pvalues
output_prefix    <- paste0(dirname(mahalanobis_dist), "/spot_")
spot_file        <- "../SPOT/spot.py"
threads          <- snakemake@threads

BPPARAM <- MulticoreParam(threads, 4000)
register(BPPARAM)

#' 
#' function to run spot in parallel
#' 
run_spot_parallel <- function(clusters, mat){
    # prepare input
    tmpmat <- mat[gsub(".*:", "", CLUSTER_ID) %in% clusters,]
    tinfile <- tempfile()
    fwrite(tmpmat, file=tinfile, sep="\t")
    
    # tmp output prefix
    toutfile <- tempfile()
    
    # run spot
    cmd <- sprintf(paste0("python %s --juncfile %s --outprefix %s"),
            spot_file, tinfile, toutfile)
    message("Running SPOT command: ", cmd)
    system(cmd)
    
    # read results
    pvals <- fread(paste0(toutfile, "_emperical_pvalue.txt"))
    mds   <- fread(paste0(toutfile, "_md.txt"))
    
    return(list(p=pvals, md=mds))
}


#' 
#' Read in count matrix
#' 
mat <- fread(intron_clusters)
setnames(mat, "V1", "CLUSTER_ID")
dim(mat)
if(any(grepl("^PSEUDO_COUNT", colnames(mat)))){
    col <- grep("^PSEUDO_COUNT", colnames(mat), value=TRUE)
    mat[,c(col):=list(NULL)]
}
colnames(mat)

#' 
#' Run SPOT in parallel
clus <- mat[,gsub(".*:", "", CLUSTER_ID)]
unique(clus)

chunks2process <- chunk(unique(clus), chunk.size=50)
res <- bplapply(chunks2process, run_spot_parallel, mat=mat, BPPARAM=BPPARAM)


#' 
#' save results
#' 
res_md   <- rbindlist(lapply(res, "[[", "md"))
res_pval <- rbindlist(lapply(res, "[[", "p"))

fwrite(res_md,   sep="\t", file=mahalanobis_dist)
fwrite(res_pval, sep="\t", file=pvalues)
