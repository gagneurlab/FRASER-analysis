#'---
#' title: Results of leafcutterMD analysis
#' author: Christian Mertes
#' wb:
#'  input:
#'   - rlds:   '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/rlds_obj.RDS"`'
#'   - jpvals: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/leafcutter_MD_pVals.txt"`'
#'   - cpvals: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/leafcutter_MD_clusterPvals.txt"`'
#'   - effect: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/leafcutter_MD_effSize.txt"`'
#'  output:
#'   - resultTable: '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/results_{dataset}.tsv"`'
#'   - wBhtml: 'Output/html/Leafcutter/{dataset}_results_MD.html'
#'  type: noindex
#'  threads: 4
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Uterus")
    parseWBHeader2("Scripts/LeafcutterMD/02_extract_results.R",
            wildcards=wildcards, rerun=TRUE)
}

#+ load config and setup, echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
dataset    <- snakemake@wildcards$dataset
rldsFile   <- snakemake@input$rlds
jpvalsFile <- snakemake@input$jpvals
cpvalsFile <- snakemake@input$cpvals
effectFile <- snakemake@input$effect
outFile    <- snakemake@output$resultTable
FDR_LIMIT  <- snakemake@config$FDR_LIMIT
threads    <- snakemake@threads
BPPARAM    <- MulticoreParam(threads)
register(BPPARAM)


#+ read input files
rls <- readRDS(rldsFile)
cpvals <- fread(cpvalsFile)
colnames(cpvals)[1] <- "clusterID"
jpvals <- fread(jpvalsFile)
colnames(jpvals)[1] <- "clusterID"
effect <- fread(effectFile)
colnames(effect)[1] <- "clusterID"

#'
#' get annotation
clusterGeneMap <- clusterGeneMapping(rls, countFile=getUnmodCountFile(rls))
dtClusterGeneMap <- as.data.table(clusterGeneMap)[,
        .(chr=gsub(":.*", "", clu), clu=gsub(".*:", "", clu), genes)][,
        .(geneID=unlist(strsplit(genes, ","))),by="chr,clu"]


#' 
#' build results
#' 
resp <- as.data.table(pivot_longer(jpvals, -clusterID, names_to="sample", values_to="pvalue"))
rese <- as.data.table(pivot_longer(effect, -clusterID, names_to="sample", values_to="effect"))
resc <- as.data.table(pivot_longer(cpvals, -clusterID, names_to="sample", values_to="cpvalue"))


#' 
#' merge with genes
res <- merge(resc, dtClusterGeneMap, by.x="clusterID", by.y="clu", allow.cartesian=TRUE)


#' 
#' correct p values
#' 
#' * FWER correction by gene
#' * BY corecction genome-wide
#' 
res_padj <- rbindlist(bplapply(res[,unique(sample)], data=res[!is.na(geneID)], 
        BPPARAM=BPPARAM, FUN=function(x, data){
            data <- data[sample == x]
            data <- data[,.(pvalue_gene=min(p.adjust(cpvalue, method="holm"))), by="sample,geneID"]
            data
}))
res_padj[,padj:=p.adjust(pvalue_gene, method="BY"), by=sample]


#' 
#' add effect 
#' 
#' For any junctions with  p-value < 1e-3
#' take: max(abs(effect))
#' 
resp[,effect:=rese[,effect]]
resps <- resp[pvalue < 1e-3]
resps[,clu:=gsub(".+:", "", clusterID)]
res_junctions <- merge(resps, dtClusterGeneMap, by="clu", allow.cartesian=TRUE)
effect2merge <- res_junctions[,.(maxEffect=effect[max(abs(effect)) == abs(effect)][1]), by="sample,geneID"]

res_all <- merge(res_padj, effect2merge, by=c("sample", "geneID"), all.x=TRUE)
res_all[is.na(maxEffect), maxEffect:=0]

DT::datatable(res_all[padj < 0.1])

#'
#' Write result table tsv
write_tsv(res_all, outFile)

