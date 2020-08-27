#'---
#' title: GTEx Rare variant enrichtment
#' author: Christian Mertes
#' wb:
#'   threads: 9
#'   input:
#'     - variantTable:    '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}_filtered_VariantsTable.tsv.gz"`'
#'     - multiTOC:        '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{snptype}__{deltaPsi}_reproducability.RDS"`'
#'     - results_all:     '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__outlierStatus__{deltaPsi}.tsv.gz"`'
#'     - resultsByMethod: '`sm [ config["DATADIR"] + "/processedData/leafcutter/{dataset}/results_{dataset}.tsv", config["DATADIR"] + "/processedData/leafcutter/{dataset}/leafcutterMD_testing/results_{dataset}.tsv", config["DATADIR"] + "/processedData/spot/{dataset}/spot__fullResults.tsv", expand(config["DATADIR"] + "/processedData/results/{{dataset}}/{method}_results.tsv", method=config["GTExMethods"]) ]`'
#'   output:
#'     - rds:             '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}.RDS"`'
#'     - ggplotFile:      '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}_ggplot.RDS"`'
#'     - wBhtml:   '`sm config["htmlOutputPath"] + "/GTEx_variant_enrichment/{dataset}__{snptype}__{deltaPsi}.html"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")

# Interactive mode (debuging)
if(FALSE){
    load_wbuild()
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(snptype="rareSplicing", dataset="Uterus", deltaPsi="0.1")
    options   <- c("--configfile", "wbuild.yaml")
    
    wildcards <- list(snptype="rareSplicing", dataset="Adrenal_Gland", deltaPsi="0.0")
    options   <- c("--configfile", "wbuild_small.yaml")
    
    parseWBHeader2("Scripts/GTExEnrichment/04_GTExSNPEnrichment.R",
            wildcards=wildcards, options=options, rerun=TRUE)
}

curAEVersion <- snakemake@config$AE_IMPLEMENTATION

tissue     <- snakemake@wildcards$dataset
snptype    <- snakemake@wildcards$snptype
deltaPsi   <- as.numeric(snakemake@wildcards$deltaPsi)
vcfFile    <- snakemake@input$variantTable
multiTOCFile <- snakemake@input$multiTOC
res_all_files  <- snakemake@input$results_all
res_by_method_files <- snakemake@input$resultsByMethod
anno_files <- snakemake@input$annotations
threads    <- snakemake@threads

outRDS     <- snakemake@output$rds
ggplotFile <- snakemake@output$ggplotFile

MIN_DELTA_PSI <- deltaPsi
MAF_LIMIT     <- 0.05
MULTI_TISSUE_CUTOFF <- c(2, 10)
BPPARAM       <- MulticoreParam(min(4, threads))
register(BPPARAM)

vcfFile
multiTOCFile
res_all_files
res_by_method_files
outRDS

###########################################
#'
#' # Load VCF file
#'
###########################################

#' Read in vcf file
#+ read vcf
variants <- fread(vcfFile)
tibble(variants)

###########################################
#'
#' # Read in multi tissue outlier calls
#'
###########################################

dt_mtoc <- readRDS(multiTOCFile)$datatable$dt2p
dt_mtoc_wide <- pivot_wider(
    dt_mtoc[,.(subjectID, geneID, Method, numMultiCall=hits3)],
    names_from="Method",
    names_glue="{Method}_{.value}",
    values_from=numMultiCall,
    values_fill=FALSE) %>%
    as.data.table()
tibble(dt_mtoc_wide)

###########################################
#'
#' # Read in outlier calls
#'
###########################################

res_all <- fread(res_all_files)
tibble(res_all)
methodZ <- grep('_z$', names(res_all), value=TRUE)
methodP <- grep('_n?p$', names(res_all), value=TRUE)

total_test <- sapply(methodP, function(x) sum(!is.na(res_all[,get(x)])))
sort(total_test)

###########################################
#'
#' # Set unified factors across tables
#'
###########################################

usamples <- unique(c(res_all$subjectID, levels(dt_mtoc_wide$subjectID)))
ugenes   <- unique(c(res_all$geneID, levels(dt_mtoc_wide$geneID)))
res_all[ ,subjectID :=factor(subjectID, levels=usamples)]
res_all[ ,geneID    :=factor(geneID,    levels=ugenes)]
dt_mtoc_wide[ ,subjectID :=factor(as.character(subjectID, levels=usamples))]
dt_mtoc_wide[ ,geneID    :=factor(as.character(geneID,    levels=ugenes))]
variants[,subjectID :=factor(as.character(subjectID, levels=usamples))]
variants[,geneID    :=factor(as.character(SYMBOL,    levels=ugenes))]
variants[,IMPACT    :=factor(IMPACT, levels=c("HIGH", "MODERATE", "LOW"))]

###########################################
#'
#' # Merge variants with outlier calls
#'
###########################################

#'
#' ## Get tested genes for all methods (FRASER, SPOT, LEAFCUTTER)
#' 
overlapGenesTested <- !rowAnyNAs(as.matrix(res_all[,c(methodZ, methodP), with=FALSE]))
res_all <- res_all[overlapGenesTested]
length(unique(res_all[,subjectID]))
length(unique(res_all[,geneID]))

#'
#' # Merge data sets into one table
#+ merge data
var2merge <- unique(variants[, .(subjectID, geneID, simple_conseq, MAF, IMPACT)])
setkey(var2merge, subjectID, geneID, IMPACT)
var2merge <- var2merge[!duplicated(var2merge, by=c("subjectID", "geneID"))]

featuresSubset <- merge(merge(
    res_all, var2merge, all.x=TRUE, by=c("subjectID", "geneID")),
        dt_mtoc_wide, by=c("subjectID", "geneID"), all.x=TRUE)
featuresSubset[, back_simple_conseq:=simple_conseq]

featuresSubset[, table(is.na(IMPACT))]
featuresSubset[, table(is.na(tissue))]
featuresSubset[, table(is.na(tissue) | is.na(IMPACT))]


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
dt <- melt(featuresSubset, measure.vars=methodZ, variable.name='Method', value.name='Zscore')
ggplots[['ZscorePerVarType']] <- ggplot(dt[abs(Zscore) > 2], aes(simple_conseq, Zscore, col=Method)) +
    stat_summary(fun.data=function(x, ...){c(y=dt[,min(Zscore, na.rm=TRUE)*1.2], label=length(x))},
                 position=position_dodge(width=.85), angle=90, geom="text") +
    geom_boxplot(position=position_dodge(width=.85)) +
    geom_violin( position=position_dodge(width=.85)) +
    grids() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
ggplots[['ZscorePerVarType']]



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
        dt1 <- calculateEnrichment(featuresSubset, cols=x, cutOff=et$cutOff,
                isZscore=et$isZscore)
        ans <- dt1$dt[, .(cutoff, nRareEvent, total, fraction, nNA, Method=x,
                enrichment=dt1$enrichment, min.ci=dt1$min.ci, max.ci=dt1$max.ci,
                multiCall=1)]
        
        if(paste0(x, "_numMultiCall") %in% colnames(featuresSubset)){
            for(i in MULTI_TISSUE_CUTOFF){
                dtmc <- calculateEnrichment(featuresSubset, 
                        cols=x, cutOff=et$cutOff,
                        isZscore=et$isZscore, multiCall=i)
                ans <- rbind(ans, dtmc$dt[,.(
                        cutoff, nRareEvent, total, fraction, nNA, Method=x, 
                        enrichment=dtmc$enrichment, min.ci=dtmc$min.ci, 
                        max.ci=dtmc$max.ci, multiCall=i)])
            }
        }
        ans
    }))
    enrich_final_ls[[paste0(tissue, ": ", et$name)]] <- enrichdt

    ggplots[[paste0('enrichAll_', et$name)]] <-
        ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
            geom_point() +
            facet_grid(rows="multiCall", scales="free_y") + 
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

ggplots[['scatterAE_LeafcutterMD']] <- ggplot(featuresSubset,
            aes(x=-log10(LeafcutterMD_p), y=-log10(get(paste0(curAEVersion, '_p'))))) +
    geom_hex(bins=50) +
    geom_abline(intercept = 0, slope = 1, col='red') +
    scale_fill_gradient(trans = "log") +
    ylab(paste0("-log10(", curAEVersion, ")")) +
    ggtitle('Scatter FraseR versus LeafcutterMD genewise P-values') +
    grids(color="white") +
    coord_fixed()

ggplots[['scatterAE_SPOT']] <- ggplot(featuresSubset,
            aes(x=-log10(SPOT_p), y=-log10(get(paste0(curAEVersion, '_p'))))) +
    geom_hex(bins=50) +
    geom_abline(intercept = 0, slope = 1, col='red') +
    scale_fill_gradient(trans = "log") +
    ylab(paste0("-log10(", curAEVersion, ")")) +
    ggtitle('Scatter FraseR versus SPOT genewise P-values') +
    grids(color="white") +
    coord_fixed()

ggplots[['scatterAE_PCA']]
ggplots[['scatterAE_LeafcutterMD']]
ggplots[['scatterAE_SPOT']]


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
saveRDS(file=outRDS, object=list(variants=variants, 
        featuresSubset=featuresSubset,
        enrich_final_ls=enrich_final_ls,
        recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))
saveRDS(ggplots, ggplotFile)

