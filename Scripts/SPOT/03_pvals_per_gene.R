#'---
#' title: Get gene level pvalues from SPOT results
#' author: Ines Scheller
#' wb:
#'  input:
#'   - rls:              '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/rlds_obj.RDS"`'
#'   - mahalanobis_dist: '`sm config["DATADIR"] + "/processedData/spot/{dataset}/spot__md.txt.gz"`'
#'   - pvalues:          '`sm config["DATADIR"] + "/processedData/spot/{dataset}/spot__emperical_pvalue.txt.gz"`'
#'  output:
#'   - gene_pvalues:  '`sm config["DATADIR"] + "/processedData/spot/{dataset}/spot__fullResults.tsv"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/SPOT/{dataset}/spot_gene_pvals.html"`'
#'  type: noindex
#'---

#+ load config and setup, echo=FALSE
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")

if(FALSE){
    load_wbuild()
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="Skin_Not_Sun_Exposed_Suprapubic")
    parseWBHeader2("Scripts/SPOT/03_pvals_per_gene.R",
            wildcards=wildcards, rerun=TRUE)
}

#+ input
rls_file         <- snakemake@input$rls
mahalanobis_dist <- snakemake@input$mahalanobis_dist
spot_pvals       <- snakemake@input$pvalues
dataset          <- snakemake@wildcards$dataset

#+ output
outFile <- snakemake@output$gene_pvalues

#+ read in distance/effect
mdvals <- fread(mahalanobis_dist)
mdvals

#+ read in pvalues
pvals <- fread(spot_pvals)
dim(pvals)
colnames(pvals)

# merge p values with effect
res <- melt(pvals, id.vars="CLUSTER_ID", value.name="p", 
        variable.name = "SAMPLE_ID")
res_md <- melt(mdvals, id.vars="CLUSTER_ID", value.name="md", 
            variable.name = "SAMPLE_ID")
res[,md:=res_md$md]

#+ get cluster to gene mapping
rls <- readRDS(rls_file)
clusterGeneMap <- clusterGeneMapping(rls, countFile=getUnmodCountFile(rls))

# one line per gene
dtClusterGeneMap <- as.data.table(clusterGeneMap)[,
        .(chr=gsub(":.*", "", clu), CLUSTER_ID=gsub(".*:", "", clu), genes)][,
        .(GENE_ID=unlist(strsplit(genes, ","))),by="chr,CLUSTER_ID"]


#+ add genes for the clusters 
res_merged <- merge(res, dtClusterGeneMap, all.x=TRUE, by="CLUSTER_ID", allow.cartesian=TRUE)

#+ adjust p values per gene
res_by_gene <- res_merged[,.(
        gene_md=max(md), 
        gene_p=min(p.adjust(p, "holm"), 1, na.rm=TRUE)),
            by="SAMPLE_ID,GENE_ID"]
res_by_gene[,gene_fdr:=p.adjust(gene_p, "BY"), by="SAMPLE_ID"]
res_by_gene

#+ write results file
fwrite(res_by_gene, file=outFile)
