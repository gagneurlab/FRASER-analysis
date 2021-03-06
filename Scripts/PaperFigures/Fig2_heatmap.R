#'---
#' title: Paper figure 02 (Heatmaps)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - fds:    '`sm expand(config["DATADIR"] + "/datasets/savedObjects/{dataset}__" + config["AE_IMPLEMENTATION"] + "/delta_psiSite.h5", dataset=config["heatmap_tissues"])`'
#'   - stats:  '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_stats.RDS", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/Figure2_heatmap_{psiType}.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/Figure2_heatmap_{psiType}.pdf"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/PaperFigures/heatmaps/heatmapFigure_{psiType}.html"`'
#'  type: noindex
#'  threads: 5
#'---

#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=10, fig.height=10)

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(psiType="psi5")
    opt <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/Fig2_heatmap.R",
            wildcards=wildcards, options=opt, rerun=TRUE)
}

#+ input
statFiles   <- snakemake@input$stats
fdsFiles    <- snakemake@input$fds[1:3]
datasets    <- basename(dirname(fdsFiles))
workingDir  <- dirname(dirname(dirname(fdsFiles[1])))
ptype       <- snakemake@wildcards$psiType
BPPARAM     <- MulticoreParam(5)

#+ output
outPng      <- snakemake@output$outPng
outPdf      <- snakemake@output$outPdf

#'
#' # Load datasets
#+ echo=TRUE
datasets
outPng

#+ echo=FALSE
fds_ls <- lapply(datasets, loadFraserDataSet, dir=workingDir)
names(fds_ls) <- datasets

stats_ls <- lapply(statFiles, readRDS)
names(stats_ls) <- basename(dirname(statFiles))


#'
#'# create annotation coloring scheme
#'
annotation_cols_2plot <- list(
    DTHHRDY   = c("DTHHRDY", "Set2"),
    AGE       = c("AGE",     "white royalblue"),
    GENDER    = c("GENDER",  "Paired"),
    SMRIN     = c("RIN",     "white seagreen"),
    SMNABTCHT = c("BATCH",   "Accent"),
    SMCENTER  = c("CENTER",  "Set2"))

annotation_full_col <- rbindlist(lapply(fds_ls, function(fds){
    as.data.frame(colData(fds)[,names(annotation_cols_2plot)]) }))
colorSets <- sapply(annotation_cols_2plot, "[[", 2)

#'
#' Correct naming
#'
annotation_full_col <- SMNABTCHT4plot(annotation_full_col)
annotation_full_col <- SMRIN4plot(annotation_full_col)
annotation_full_col <- GENDER4plot(annotation_full_col)
colnames(annotation_full_col) <- sapply(annotation_cols_2plot, "[[", 1)
names(colorSets) <- colnames(annotation_full_col)


#+ manually create annotation colors
needed_colors <- apply(annotation_full_col, 2, function(x){ length(levels(as.factor(x)))})
annotation_colors <- lapply(names(needed_colors), function(name){
    nrValues <- needed_colors[name]
    if(!grepl(" ", colorSets[name])){
        res <- brewer.pal(max(3, nrValues), colorSets[name])
        if(nrValues < 3){
            res <- res[1:nrValues]
        }
        if(any(is.na(annotation_full_col[,get(name)]))){
            res <- c(res, "white")
        }
    } else {
        res <- colorRampPalette(strsplit(colorSets[name], " ")[[1]])(nrValues)
    }
    names(res) <- levels(factor(paste0("", annotation_full_col[,get(name)])))
    res
})
names(annotation_colors) <- names(colorSets)


#'
#' # Create heatmaps
#'
#+ creating heatmaps
topN <- 30000
topJ <- 10000
heatmap <- bplapply(names(fds_ls), BPPARAM=BPPARAM, function(tissue) {
    
    message(date(), ": ", tissue)
    fds <- fds_ls[[tissue]]
    # fds <- fds[seqnames(fds) == "21", 1:30]
    
    colData <- colData(fds)[,names(annotation_cols_2plot)]
    colData <- SMNABTCHT4plot(colData)
    colData <- SMRIN4plot(colData)
    colData <- GENDER4plot(colData)
    colnames(colData) <- sapply(annotation_cols_2plot, "[[", 1)
    
    for(i in colnames(colData)){
        colData[,i] <- factor(colData[,i], levels=levels(
            factor(annotation_full_col[,get(i)])))
    }
    colData(fds) <- cbind(colData(fds)[,c("sampleID", "bamFile")], colData)
    
    x <- plotCountCorHeatmap(
        fds = fds,
        type = ptype,
        logit = TRUE,
        topN = topN,
        topJ = topJ,
        plotType = "junctionSample",
        sampleClustering = NA,
        normalized = FALSE,
        annotation_col = colnames(colData),
        annotation_row = NA,
        annotation_colors = annotation_colors,
        sampleCluster = NA,
        plotMeanPsi = FALSE,
        plotCov = FALSE,
        annotation_legend = TRUE
    )
    
    x
})
heatmap[[1]]


#'
#' # Correlation reduction
#'
corDT <- rbindlist(lapply(stats_ls, "[[", "SampleCors"))[type == ptype]
levels <- corDT[,sort(unique(dataset))[c(4:6,1:3,7:length(unique(dataset)))]]
corDT[,dataset:=factor(dataset, levels=levels)]
levels(corDT$dataset) <- dName4plot(levels(corDT$dataset))
corDT[,method:=factor(method)]
corDT[,normalized:=factor(normalized, levels=c("raw", "normalized"))]


# corDT <- corDT[sample(.N, 100000),]

corBoxplots <- ggplot(corDT, aes(y=abs(cor), x=dataset, fill=normalized)) +
    geom_boxplot() +
    theme_bw() +
    theme(
            axis.text.x = element_text(angle=35, hjust=1),
            axis.text   = element_text(size=11),
            axis.title  = element_text(face="bold", size=14),
            legend.text = element_text(size=14) ) +
    ggtitle("") +
    ylab("|Sample correlation|") +
    xlab("GTEx tissue") +
    scale_fill_brewer(palette = "Dark2", direction = -1,
            name="Controlled", labels=c("No", "Yes"))
corBoxplots


#'
#' Assemble figure
#'
#+ assemble full figure
get_legend_grob <- function(gt, spacing=20){
    gt <- gtable_remove_grobs(gt, c(
        "main", "col_tree", "row_tree", "matrix", "col_annotation"))
    gt <- gtable_add_rows(gt, unit(1, "bigpts"), pos=0)
    gt$heights[1] <- unit(spacing, "bigpts")
    gt
}

get_heatmap_grob <- function(gt, spacing=10){
    gt <- gtable_remove_grobs(gt, c(
        "main", "legend", "col_annotation_names", "annotation_legend"))
    gt$widths[3] <- unit(1, "grobwidth", data=gt$grobs[[3]])
    gt <- gtable_add_cols(gt, unit(1, "bigpts"))
    gt$widths[4] <- unit(spacing, "bigpts")
    gt
}


displayNames <- dName4plot(datasets, TRUE)
displayNames[grepl("Skin not Sun Exposed", displayNames)] <- "Suprapubic Skin"
g <- ggarrange(nrow=3, heights=c(8, 0.3, 5),
        ggarrange(nrow=1, widths=c(1,10,10,10,4.5), labels=c("", letters[1:3], ""),
                font.label=list(size=20, face="bold"),
                grid.text("Introns                          ", rot=90, gp=gpar(fontsize=14)),
                get_heatmap_grob(heatmap[[1]][4]$gtable),
                get_heatmap_grob(heatmap[[2]][4]$gtable),
                get_heatmap_grob(heatmap[[3]][4]$gtable, spacing=0),
                get_legend_grob(heatmap[[3]][4]$gtable, spacing=60)),
        ggarrange(ncol=5, widths=c(1.3, 9, 9, 9, 5),
                grid.text(""),
                grid.text(paste("  ",      
                        displayNames[1], "samples"), gp=gpar(fontsize=14)),
                grid.text(paste("        ",
                        displayNames[2], "samples"), gp=gpar(fontsize=14)),
                grid.text(paste("                     ",
                        displayNames[3], "samples"), gp=gpar(fontsize=14)),
                grid.text("")),
        ggarrange(ncol=1, labels=letters[4], label.y=1.05,
                font.label=list(size=20, face="bold"),
                corBoxplots + theme(legend.position="right")))
g


#+ save heatmap figure
size  <- 13
ratio <- 0.97
outPng
ggsave(outPng, g, width = size * ratio, height = size * 1/ratio)
ggsave(outPdf, g, width = size * ratio, height = size * 1/ratio)
