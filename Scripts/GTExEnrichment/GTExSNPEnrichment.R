#'---
#' title: GTEx Rare variant enrichtment
#' author: Christian Mertes
#' wb:
#'   threads: 9
#'   input:
#'     - variantTable:  '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}_filtered_VariantsTable.tsv.gz"`'
#'     - leafcutter:    '`sm config["DATADIR"] + "/processedData/leafcutter/{dataset}/results_{dataset}.tsv"`'
#'     - annotations:   '`sm config["DATADIR"] + "/annotations/{dataset}.tsv"`'
#'     - datasets_all:  '`sm expand(config["DATADIR"] + "/datasets/savedObjects/{{dataset}}__{method}/pajdBetaBinomial_psiSite.h5", method=config["GTExMethods"])`'
#'     - results_all:   '`sm expand(config["DATADIR"] + "/processedData/results/{{dataset}}/{method}_results.tsv", method=config["GTExMethods"])`'
#'   output:
#'     - rds:           '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}.RDS"`'
#'     - ggplotFile:    '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}_ggplot.RDS"`'
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}.html"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")
load_all("../rare-disease-leafcutter/")

# Interactive mode (debuging)
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(snptype="rareSplicing", dataset="Uterus", deltaPsi="0_0")
    options   <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/GTExEnrichment/GTExSNPEnrichment.R",
            wildcards=wildcards, options=options, rerun=TRUE)
}

# stop(paste0(date(), ": Stop this and dont run it ..."))
PRODUCTION <- FALSE
useRANDOM <- FALSE
curAEVersion <- snakemake@config$AE_IMPLEMENTATION

tissue     <- snakemake@wildcards$dataset
snptype    <- snakemake@wildcards$snptype
deltaPsi   <- as.numeric(snakemake@wildcards$deltaPsi)
vcfFile    <- snakemake@input$variantTable
rlc_files  <- snakemake@input$leafcutter
fds_files  <- snakemake@input$datasets_all
res_files  <- snakemake@input$results_all
anno_files <- snakemake@input$annotations
threads    <- snakemake@threads

outRDS     <- snakemake@output$rds
ggplotFile <- snakemake@output$ggplotFile

MIN_DELTA_PSI <- deltaPsi
MIN_COVERAGE  <- 10
MAF_LIMIT     <- 0.05
BPPARAM       <- MulticoreParam(min(4, threads))
register(BPPARAM)

vcfFile
rlc_files
res_files

###########################################
#'
#' # Load VCF file
#'
###########################################

#' Read in vcf file
#+ read vcf
variants <- fread(vcfFile)

###########################################
#'
#' # Read in outlier calls
#'
###########################################

#'
#' ## Read Leafcutter calls
#+ leafcutter readin
names(rlc_files) <- basename(dirname(rlc_files))
enrich_rlc <- lapply(names(rlc_files), function(x){
    rlds <- readRDS(file.path(dirname(rlc_files[[x]]), "rlds_obj.RDS"))
    clusterGeneMap <- clusterGeneMapping(rlds)
    ans1 <- collectResults(rlds, clusterGeneMap=clusterGeneMap,
            plot=FALSE, FDR_CUTOFF=2, BPPARAM=BPPARAM, save=FALSE)

    ans <- ans1[,.(subjectID=sampleIDs, geneID=genes, p=gene_p)][,
        .(p=min(p, na.rm=TRUE), tissue=x), by="subjectID,geneID"][
        order(p)]
    names(ans) <-  c('subjectID', 'geneID', "Leafcutter_p", "tissue")
    ans[,-"tissue"]
})
names(enrich_rlc) <- names(rlc_files)


#'
#' ## Read all FraseR objects
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
#' ## Read annotations
#+ anno readin
anno <- fread(anno_files)
enrich_fds <- lapply(enrich_fds, anno=anno, function(x, anno){
    setnames(x, "subjectID", "run")
    x <- merge(x, anno[,.(run=sampleID, indivID)], by="run")
    setnames(x, "indivID", "subjectID")
    x[,run:=NULL]
    x
})


#'
#' ## Get tested genes (both FraseR + Leafcutter)
# tested genes
testedGenes <- intersect(
        enrich_fds[[1]][,geneID],
        enrich_rlc[[1]][,geneID])

#'
#' # Merge data sets into one table
#+ merge data
tissueData <- CJ(subjectID=anno$indivID, geneID=testedGenes)
var2merge <- unique(variants[, .(subjectID, geneID=SYMBOL, simple_conseq, MAF, IMPACT)])
var2merge <- var2merge[!duplicated(paste(subjectID, geneID))]

features2merge <- c(list(tested=tissueData, vars=var2merge), enrich_fds, enrich_rlc)
featuresSubset <- Reduce(x=features2merge, function(x, y){
        merge(x, y, by=c('subjectID', 'geneID'), all.x=TRUE) })
methodZ <- grep('_z$', names(featuresSubset), value=TRUE)
methodP <- grep('_n?p$', names(featuresSubset), value=TRUE)


#'
#' ## Number of variants we test for enrichment
#'
ggplots <- list()
ggplots[['NumVariants']] <- ggplot(featuresSubset[,.(Consequence=simple_conseq)], aes(Consequence)) + geom_bar() +
    scale_y_log10() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggplots[['NumVariants']]
table(featuresSubset[, !is.na(simple_conseq)])
featuresSubset[, sum(!is.na(simple_conseq))/.N]*100

#'
#' # Overview of enrichment
#+ Zscore overview per variant type
#+ zscore distribution, fig.width=14
dt <- melt(featuresSubset, measure.vars=methodP, variable.name='Method', value.name='Zscore')
ggplots[['ZscorePerVarType']] <- ggplot(dt[abs(Zscore) > 2], aes(simple_conseq, Zscore, col=Method)) +
    stat_summary(fun.data=function(x, ...){c(y=dt[,min(Zscore, na.rm=TRUE)*1.2], label=length(x))},
                 position=position_dodge(width=.85), angle=90, geom="text") +
    geom_boxplot(position=position_dodge(width=.85)) +
    geom_violin( position=position_dodge(width=.85)) +
    grids() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
ggplots[['ZscorePerVarType']]

#'
#' #Add random samples
if(isTRUE(useRANDOM)){
    featuresSubset[,'RAND_p':=runif(.N, 0, 1)]
    featuresSubset[,'RAND_z':=rnorm(.N, 0, 5)]
    methodP <- unique(c(methodP, "RAND_p"))
    methodZ <- unique(c(methodZ, "RAND_z"))
}
featuresSubset[, back_simple_conseq:=simple_conseq]


#'
#' # Enrichments (all)
#'
enrich_table <- list(
    list(name="Pval 1e-3", methods=methodP, isZscore=FALSE, cutOff=1e-3),
    list(name="Pval 1e-4", methods=methodP, isZscore=FALSE, cutOff=1e-4),
    list(name="Pval 1e-5", methods=methodP, isZscore=FALSE, cutOff=1e-5),
    list(name="Zscore 2", methods=methodZ, isZscore=TRUE, cutOff=2),
    list(name="Zscore 3", methods=methodZ, isZscore=TRUE, cutOff=3),
    list(name="Zscore 5", methods=methodZ, isZscore=TRUE, cutOff=5)
)

enrich_final_ls <- list()
for(et in enrich_table){
    enrichdt <- rbindlist(lapply(et$methods, function(x){
        dt <- calculateEnrichment(featuresSubset, cols=x, cutOff=et$cutOff,
                isZscore=et$isZscore)
        dt$dt[,.(cutoff, nRareEvent, total, fraction, nNA, Method=x,
                 enrichment=dt$enrichment, min.ci=dt$min.ci, max.ci=dt$max.ci)]
    }))
    enrich_final_ls[[paste0(tissue, ": ", et$name)]] <- enrichdt

    ggplots[[paste0('enrichAll_', et$name)]] <-
        ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
            geom_point() +
            geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
            coord_flip() +
            labs(title=paste0('Enrichment (All, ', tissue, ', ', et$name, ')'))
}

names2plot <- paste0('enrichAll_', sapply(enrich_table, "[[", "name"))
etl <- length(enrich_table)/2
order2plot <- rep(1:etl, each=2) + rep(c(0, etl), etl)
grid.arrange(ncol=2, grobs=ggplots[names2plot][order2plot])


#'
#' # Enrichment (rare: MAF < 0.01)
mafsub <- copy(featuresSubset)
mafsub[MAF > 0.01, simple_conseq:=NA]
for(et in enrich_table){
    enrichdt <- rbindlist(lapply(et$methods, function(x){
        dt <- calculateEnrichment(mafsub, cols=x, cutOff=et$cutOff,
                isZscore=et$isZscore)
        dt$dt[,.(cutoff, nRareEvent, total, fraction, nNA, Method=x,
                enrichment=dt$enrichment, min.ci=dt$min.ci, max.ci=dt$max.ci)]
    }))
    enrich_final_ls[[paste0(tissue, ": ", et$name, " + MAF < 0.01")]] <- enrichdt

    ggplots[[paste0('enrichMAF_0.01_', et$name)]] <-
        ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
            geom_point() +
            geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
            coord_flip() +
            labs(title=paste0('Enrichment (MAF < 0.01, ', tissue, ', ', et$name, ')'))
}

names2plot <- paste0('enrichMAF_0.01_', sapply(enrich_table, "[[", "name"))
etl <- length(enrich_table)/2
order2plot <- rep(1:etl, each=2) + rep(c(0, etl), etl)
grid.arrange(ncol=2, grobs=ggplots[names2plot][order2plot])


#'
#' * P-value scatter plot for AE versus Leafcutter and PCA
#+ scatter p-value
ggplots[['scatterAE_PCA']] <- ggplot(featuresSubset,
            aes(x=-log10(PCA_p), y=-log10(get(paste0(curAEVersion, '_p'))))) +
    geom_hex(bins=50) +
    geom_abline(intercept = 0, slope = 1, col='red') +
    scale_fill_gradient(trans = "log") +
    ylab(paste0("-log10(", curAEVersion, ")")) +
    ggtitle('Scatter FraseR versus PCA genewise P-values') +
    grids(color="white") +
    coord_fixed()

ggplots[['scatterAE_Leafcutter']] <- ggplot(featuresSubset,
            aes(x=-log10(Leafcutter_p), y=-log10(get(paste0(curAEVersion, '_p'))))) +
    geom_hex(bins=50) +
    geom_abline(intercept = 0, slope = 1, col='red') +
    scale_fill_gradient(trans = "log") +
    ylab(paste0("-log10(", curAEVersion, ")")) +
    ggtitle('Scatter FraseR versus Leafcutter genewise P-values') +
    grids(color="white") +
    coord_fixed()

ggplots[['scatterAE_PCA']]
ggplots[['scatterAE_Leafcutter']]


#'
#' # Recall Rank plots all (HIGH/MODERATE)
#'
#+ create recall data 1
#' add random
methods2plot <- union(methodP, methodZ)
rrdt <- rbindlist(bplapply(methods2plot, dt=featuresSubset, BPPARAM=BPPARAM,
        function(x, dt){
            totalHits <- sum(!is.na(dt[,simple_conseq]))
            dt <- calculateRecallRank(dt, x, grepl('_n?p$', x))
            dt <- data.table(Method=x, dt[,.(
                    nTrueHits=get(paste0(x, "_recall")),
                    recall=get(paste0(x, "_recall"))/totalHits,
                    rank=get(paste0(x, '_rank')))])
            dt[,Type:=ifelse(grepl('_p$', Method), 'P-value', 
                    ifelse(grepl('_np$', Method), 'P-value norm', 'Z-score'))]
            dt[,Method:=gsub('_n?[pz]$', '', Method)]
        }))
rrdt

#'
#' Merge with cutoffs from enrichment
#'
dt_tmp <- rbindlist(lapply(names(enrich_final_ls[1:5]), function(x){
    Cutoff <- strsplit(x, " ")[[1]][3]
    enrich_final_ls[[x]][cutoff == TRUE, .(
            rank=total,
            Method=gsub("_n?[pz]", "", Method),
            Type=ifelse(grepl('_p$', Method), 'P-value', 
                    ifelse(grepl('_np$', Method), 'P-value norm', 'Z-score')),
            Cutoff=Cutoff)]
}))
dt4cutoffs <- merge(dt_tmp, rrdt)[order(Method, Type)]


#'
#' The recall plots
#'
#+ recall plots 1
maxRank <- 5000
ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
            rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(tissue, '\n', snptype, " rank < ", maxRank)) +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed")
ggplots[[paste0('recall_n=', maxRank)]]

maxRank <- 50000
ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(tissue, '\n', snptype, " rank < ", maxRank)) +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed")
ggplots[[paste0('recall_n=', maxRank)]]

maxRank <- rrdt[,.(x=floor(log10(max(rank))))][,(x+1)*10^x]
tmplines <- rrdt[,.(x=0:max(rank), y=0:max(rank)*1/max(rank), Type="a")][x>1000 | x == 0]
ggplots[["recall_all"]] <- plotRecallRankForEnrichment(rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(tissue, '\n', snptype, " full data")) +
    grids(color="white") +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=4) +
    scale_y_log10() +
    scale_x_log10() +
    geom_line(data=tmplines, aes(x=x, y=y), col="firebrick", linetype="dashed")
ggplots[["recall_all"]]


#+ Save results
saveRDS(file=outRDS, object=list(
        enrich_fds=enrich_fds, enrich_rlc=enrich_rlc, variants=variants,
        featuresSubset=featuresSubset, enrich_final_ls=enrich_final_ls,
        recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))
saveRDS(ggplots, ggplotFile)

