#'---
#' title: Number of outlier events per method (Figure S11)
#' author: Christian Mertes
#' wb:
#'  threads: 5
#'  input:
#'   - stats:       '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_stats.RDS", dataset=config["EnrichmentTissues"])`'
#'   - lefacutter:  '`sm expand(config["DATADIR"] + "/processedData/leafcutter/{dataset}/final_results.RDS", dataset=config["EnrichmentTissues"])`'
#'   - anno:        '`sm expand(config["DATADIR"] + "/annotations/{dataset}.tsv", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS11_outlier_events_by_method.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS11_outlier_events_by_method.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS11_total_outlier_numbers.R",
            options=options, rerun=TRUE)
    threads <- 24
}

#+ input
fdsStatFiles       <- snakemake@input$stats
leafcutterResFiles <- snakemake@input$lefacutter
annoFiles          <- snakemake@input$anno
datasets           <- gsub(".tsv", "", basename(annoFiles))

names(annoFiles)    <- datasets
names(fdsStatFiles) <- datasets

threads <- snakemake@threads
workDir <- file.path(CONFIG$DATADIR, "datasets")
FraseR_implementation <- CONFIG$AE_IMPLEMENTATION
FDR_LIMIT <- 0.1

#+ output
outPdf <- snakemake@output$outPdf
outPng <- snakemake@output$outPng

#'
#' Input
length(datasets)
annoFiles[1:5]
fdsStatFiles[1:5]
leafcutterResFiles[1:5]
outPng


#' 
#' ## Load Leafcutter
# f <- leafcutterResFiles[1]
res_leafcutter <- rbindlist(mclapply(leafcutterResFiles, mc.cores=threads, 
    function(f){
            resdt <- readRDS(f)
            tissue <- basename(dirname(f))
            anno <- unique(fread(annoFiles[tissue])[,.(indivID, sampleID)])
            res <- unique(resdt[!is.na(genes),.(conditionID, genes, gene_fdr)])
            res[gene_fdr < FDR_LIMIT][order(gene_fdr)]
            res1 <- res[,.N, by=c("conditionID", "genes")][,.(
                    NumberOfEvents=.N),by="conditionID"]
            ans <- merge(all.x=TRUE,
                    anno[,.(sampleID, indivID, tissue=tissue)], 
                    res1[,.(indivID=conditionID, NumberOfEvents)])
            ans[is.na(NumberOfEvents), NumberOfEvents:=0]
            ans[,-"indivID"]
    }))


#' 
#' ## Load FRASER plotting data
#' 
fdsls <- lapply(fdsStatFiles, readRDS)
dt <- rbindlist(lapply(fdsls, "[[", "allResNumbers"))

#' 
#' filter and naming
#'
dt <- dt[quantile %in% c("50%", "90%") & minCov != 5]
dt[,labelCut:=paste("atop(",
        ifelse(is.na(padjCut), "", paste("FDR < ", padjCut, " , ")),
        ifelse(is.na(zsCut) | zsCut == 0,   "", paste("'|z score|' > ", zsCut, " , ")),
        ifelse(is.na(dPsiCut), "", paste("Delta*psi > ", dPsiCut, " , ")),
        #ifelse(is.na(minCov),  "", paste("coverage >= ", minCov)),
        " )")]
dt

dt[,type:=factor(type, levels=c("psi5", "psi3", "psiSite"))]
levels(dt$type) <- c("psi[5]", "psi[3]", "theta")
# dt$type <- rep(rep(c("psi[5]", "psi[3]", "theta"), each=2), 42)

g1 <- ggplot(dt, aes(y=outlier, x=quantile)) + 
    geom_violin() + 
    facet_grid(type ~ labelCut, labeller=label_parsed) + 
    scale_y_log10() + 
    theme_bw() +
    ylab("Outlier calls") + 
    xlab("Quantile per tissue")
g1



#'
#' Assemble figures
#'
#+ assemble figure for psi5
g <- g1
g

#+ save heatmap figure
factor <- 0.37
outPng
ggsave(outPng, g, width = 17*factor, height = 13*factor)
ggsave(outPdf, g, width = 17*factor, height = 14*factor)

