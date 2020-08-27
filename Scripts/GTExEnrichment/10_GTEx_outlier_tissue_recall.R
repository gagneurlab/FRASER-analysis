#'---
#' title: GTEx Rare variant enrichtment Outlier extraction (per tissue)
#' author: Christian Mertes
#' wb:
#'   threads: 9
#'   input:
#'     - variants:      '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}_filtered_VariantsTable.tsv.gz"`'
#'     - annotation:    '`sm config["GTEX_SAMPLE_ANNOTATION_FILE"]`'
#'     - outlierCalls:  '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__outlierStatus__{{deltaPsi}}.tsv.gz", dataset=config["EnrichmentTissues"])`'
#'   output:
#'     - plots:         '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}__{deltaPsi}_reproducability.RDS"`'
#'     - allres:        '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}__{deltaPsi}_reproducability.tsv.gz"`'
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_variant_enrichment/{snptype}__{deltaPsi}_reproducability.html"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")
library(gplots)

# Interactive mode (debuging)
if(FALSE){
    load_wbuild()
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(snptype="rareSplicing", deltaPsi="0.0")
    options   <- c("--configfile", "wbuild_small.yaml")
    parseWBHeader2("Scripts/GTExEnrichment/10_GTEx_outlier_tissue_recall.R",
                   wildcards=wildcards, options=options, rerun=TRUE)
}

minTissuesTested  <- 10
goodMethods       <- c("BB_p", "LeafcutterMD_p", "SPOT_p", "PCA_p")
variant_file      <- snakemake@input$variants
annotation_file   <- snakemake@input$annotation
outlierCall_files <- snakemake@input$outlierCalls

outfile    <- snakemake@output$plots
outResFile <- snakemake@output$allres

#' 
#' The user input
#' 
goodMethods
annotation_file
variant_file
outlierCall_files
outResFile


#' 
#' prepare variant table
#' 
variants <- fread(variant_file)
vars <- variants[,
        .(subjectID, geneID=SYMBOL, 
                IMPACT=factor(IMPACT, levels=c("HIGH", "MODERATE", "LOW")))][
        order(subjectID, geneID, IMPACT)]
vars <- vars[!duplicated(vars, by=c("subjectID", "geneID"))]
tibble(vars)
variants[,table(IMPACT)]

#'
#' preprare subsetting (sample / method)
#' 
anno <- fread(annotation_file)
annoTable <- table(anno[grepl("GTEX-", SAMPID),.(SAMPID, SUBJECT=gsub("^", "GTEX-", 
        gsub("-.*", "", gsub("GTEX-", "", SAMPID))), SMTSD)][,.(SUBJECT, SMTSD)])
annoTable[annoTable > 0] <- 1
goodTissues <- names(which(colSums(annoTable) >= 20))
annoTable <- annoTable[,goodTissues]

plot(ecdf(rowSums(annoTable)), main="ECDF of number of tissues per subject")
grid()
goodSamples <- names(which(rowSums(annoTable) >= 20))
length(goodSamples)
heatmap.2(annoTable[goodSamples,], trace="none", srtCol=-90, offsetCol=0, adjCol=0)

goodMethods

#' 
#' prepare tissue prediction table
#' 
BPPARAM <- MulticoreParam(10, progressbar=TRUE)
transdtls <- bplapply(outlierCall_files, BPPARAM=BPPARAM, FUN=function(f){
    x <- fread(f)
    cur_tissue <- x[,unique(tissue)]
    x <- x[subjectID %in% goodSamples]
    x <- x[!is.na(geneID) & geneID != ""]
    dtpl <- x[,c("subjectID", "geneID", goodMethods), with=FALSE]
    
    dttrans <- pivot_longer(dtpl,
            cols=all_of(goodMethods),
            names_to="Method",
            values_to=cur_tissue) %>%
        # mutate(!!cur_tissue := replace(.[[cur_tissue]], is.na(.[[cur_tissue]]), 2)) %>%
        as.data.table()
    dttrans[!is.na(cur_tissue)]
})

ugenes   <- unique(unlist(lapply(transdtls, function(x) unique(x[["geneID"]]))))
usubject <- unique(unlist(lapply(transdtls, function(x) unique(x[["subjectID"]]))))
umethod  <- unique(unlist(lapply(transdtls, function(x) unique(x[["Method"]]))))

length(ugenes)
length(usubject)
length(umethod)
length(ugenes) * length(usubject) * length(umethod)
dt2merge <- as.data.table(expand.grid(subjectID=factor(usubject), geneID=factor(ugenes), Method=factor(umethod)))
setkey(dt2merge, subjectID, geneID, Method)

mergeResTables <- function(x, dt2m, withMethod=TRUE){
    x[,subjectID:=factor(subjectID, levels=usubject)]
    x[,geneID:=factor(geneID, levels=ugenes)]
    if(isTRUE(withMethod)){
        x[,Method:=factor(Method, levels=umethod)]
        setkey(x, subjectID, geneID, Method)
    } else {
        setkey(x, subjectID, geneID)
    }
    ans <- x[dt2m]
    ans
}
mergels <- bplapply(transdtls, FUN=mergeResTables, dt2m=dt2merge, BPPARAM=BPPARAM)

res <- dt2merge
for(i in seq_along(mergels)){
    tissue_name <- colnames(mergels[[i]])[4]
    message(date(), "Run iteration: ", i, " and tissue: ", tissue_name)
    res[,c(tissue_name):=list(mergels[[i]][[tissue_name]])]
}

#' merge variants
res <- mergeResTables(vars, res, withMethod=FALSE)

#' 
#' clean up
#' 
res <- res[!is.na(Method)]
res[,totalTested:=NA_integer_]
setcolorder(res, unique(c("subjectID", "geneID", "Method", "totalTested", "IMPACT"), colnames(res)))
res

resMat <- as.matrix(res[,-(1:5)])
total <- rowSums(!is.na(resMat))
res[,totalTested:=total]

gt1 <- ggplot(res[totalTested > 0,.(Method, totalTested)], aes(totalTested, fill=Method)) + 
    geom_histogram(position="dodge") + 
    scale_fill_brewer(palette="Dark2") +
    ylab("Number of tested events") + 
    xlab(bquote("Number of tissues event is tested")) + 
    theme_bw() + 
    grids() + 
    scale_y_log10()
gt1


resMat4p <- resMat[total >= minTissuesTested,]
dt2p <- res[totalTested >= minTissuesTested, .(subjectID, geneID, Method, IMPACT, totalTested)]
dt2p

dt2p$hits3 <- rowSums(resMat4p < 1e-3, na.rm = TRUE)
dt2p$hits5 <- rowSums(resMat4p < 1e-5, na.rm = TRUE)
dt2p$hits7 <- rowSums(resMat4p < 1e-7, na.rm = TRUE)
dt2p$hits9 <- rowSums(resMat4p < 1e-9, na.rm = TRUE)
dt2p[,spliceVariant:=ifelse(!is.na(IMPACT), "rare splice variant", "no rare variant")]
setcolorder(dt2p, unique(c("subjectID", "geneID", "Method", "totalTested", "IMPACT"), colnames(dt2p)))
setkey(dt2p, subjectID, geneID, Method)
dt2p


g1 <- ggplot(dt2p[,.(hits5=ifelse(hits5 >= 1, hits3, 0), Method, spliceVariant)][hits5 != 0]
            , aes(x=hits5, fill=Method)) + 
    geom_bar(position="dodge") + 
    theme_bw() + 
    facet_wrap(~spliceVariant) + 
    scale_fill_brewer(palette="Dark2") + 
    ylab("Number of events") + 
    xlab("Number of tissues outlier is present")
g1 + 
    theme_cowplot() + 
    grids()

g2 <- ggplot(dt2p[,.(hits7=ifelse(hits7 >= 1, hits3, 0), Method, spliceVariant)][hits7 != 0],
            aes(x=hits7, fill=Method)) + 
    geom_bar(position="dodge") + 
    theme_bw() + 
    facet_wrap(~spliceVariant) + 
    scale_fill_brewer(palette="Dark2") + 
    ylab("Number of events") + 
    xlab("Number of tissues outlier is present")
g2

getPercentage <- function(dt, col, rep_col){
    numHits <- dt[,.(hits=get(col), rep=get(rep_col))][,ifelse(hits > 0, rep, 0)]
    if(is.factor(numHits)){
        numHits <- levels(numHits)
    } else {
        numHits <- as.character(seq_len(max(numHits)))
    }
    nameFac <- factor(numHits, levels=numHits)
    dt2p2 <- dt[,.(hits=factor(get(col), levels=numHits), Method)][,
            .(counts=as.vector(table(hits)), names=nameFac), by=Method]
    dt2p2[, freq:=counts/sum(counts), by=Method]
    dt2p2
}

plotPercentage <- function(dt, value){
    ggplot(dt, aes(y=freq*100, x=names, fill=Method)) + 
        geom_bar(stat="identity", position="dodge") + 
        labs(x=paste0("Number of tissues outlier is present (p < 1e-", value, ")"),
            y="Percentage within Method") + 
        scale_fill_brewer(palette="Dark2") + 
        theme_cowplot() + 
        grids()
}

dt2p9p <- getPercentage(dt2p, "hits9", "hits3")
dt2p9p[names == "1"]
gg9p <- plotPercentage(dt2p9p, "9")
gg9p

dt2p7p <- getPercentage(dt2p, "hits7", "hits3")
dt2p7p[names == "1"]
gg7p <- plotPercentage(dt2p7p, "7")
gg7p

dt2p5p <- getPercentage(dt2p, "hits5", "hits3")
gg5p <- plotPercentage(dt2p5p, "5")
gg5p

ggarrange(gg9p, gg7p, gg5p, ncol=1, common.legend=TRUE)


#' 
#' # saving the results
#' 
allobject <- list(
    ggplots = list(
        g1,
        g2,
        gg9p,
        gg7p,
        gg5p
    ),
    datatables = list(
        dt2p=dt2p,
        dt2p9p=dt2p9p,
        dt2p7p=dt2p7p,
        dt2p5p=dt2p5p
    )
)
saveRDS(allobject, outfile)
fwrite(res, sep="\t", file=outResFile)
