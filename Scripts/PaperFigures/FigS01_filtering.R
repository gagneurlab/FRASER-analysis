#'---
#' title: Paper figure S01 junction filtering
#' author: Christian Mertes
#' wb:
#'  input:
#'   - rawFDS:    '`sm expand(config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/fds-object.RDS", dataset=["Kremer"] + config["heatmap_tissues"])`'
#'   - html:      '`sm expand(config["htmlOutputPath"] + "/FraseR/{dataset}_filterExpression.html", dataset=["Kremer"] + config["heatmap_tissues"])`'
#'   - statsGTEx: '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_stats.RDS", dataset=config["EnrichmentTissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS1_filtering.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS1_filtering.pdf"`'
#'  threads: 5
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list()
    opt <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS01_filtering.R",
            wildcards=wildcards, options=opt, rerun=TRUE)
}

#+ input
fdsFiles    <- snakemake@input$rawFDS
statsFiles  <- snakemake@input$statsGTEx
datasets    <- basename(dirname(fdsFiles))
workingDir  <- dirname(dirname(dirname(fdsFiles[1])))
BPPARAM     <- MulticoreParam(5)

#+ output
outPng      <- snakemake@output$outPng
outPdf      <- snakemake@output$outPdf

#'
#' # Load datasets
#+ echo=TRUE
datasets

#+ echo=FALSE
fds_ls <- lapply(datasets, loadFraseRDataSet, dir=workingDir)
names(fds_ls) <- datasets

stats <- lapply(statsFiles, readRDS)
names(stats) <- dName4plot(basename(dirname(statsFiles)))

#'
#' # Plot filter expression per dataset
plot_ls <- lapply(fds_ls, plotFilterExpression)
plot_ls <- lapply(names(fds_ls), function(x){ 
        plot_ls[[x]] + ggtitle(dName4plot(x)) })
names(plot_ls) <- datasets

idx <- grepl("raw-Skin_Not_Sun_Exposed_Suprapubic", names(plot_ls))
if(any(idx)){
    plot_ls[[which(idx)]] <- plot_ls[[which(idx)]] + ggtitle("Suprapubic Skin")
}

#'
#' Plot Number of expressed genes
#'
statsdt <- as.data.table(pivot_longer(
    as.data.table(sapply(stats, "[[", "ImportantNumbers"), keep.rownames=TRUE),
    -rn, names_to="tissue", values_to="value"))
nameMapping <- c(
        "NJunctions"="Number of expressed\nintrons",
        "NSites"="Number of expressed\nsplice sites",
        "NRawJunctions"="Total number of observed\nintrons",
        "Nsamples"="Number of samples")
dt2plot <- statsdt[rn %in% names(nameMapping)]
dt2plot[,rn:=factor(nameMapping[rn], levels=nameMapping)]
dt2plot <- dt2plot[order(tissue, decreasing=TRUE)]
dt2plot[,tissue:=factor(tissue, levels=sort(unique(tissue), decreasing=TRUE))]
g1 <- ggplot(dt2plot, aes(x=tissue, y=value, fill=rn)) +
    geom_bar(stat="identity") +
    coord_flip() +
    facet_wrap(~rn, scales="free_x", nrow=1) +
    theme_bw() +
    theme(legend.position="none", axis.title.x = element_blank()) +
    scale_fill_brewer(palette="Dark2") +
    xlab("")
g1

#'
#' Assemble figure
#'
#+ assemble full figure
g <-ggarrange(ncol=1, nrow=2, legend="none", heights = c(1,2.5),
    ggarrange(labels=LETTERS[1:4], nrow=1, align="hv", legend="none",
        plot_ls[[1]] + annotation_custom(
            get_legend(plot_ls[[1]]),
                xmin=log10(2000),  xmax=log10(100000),
                ymin=log10(200000), ymax=log10(2000000)),
        plot_ls[[2]],
        plot_ls[[3]],
        plot_ls[[4]]),
    ggarrange(labels=LETTERS[5], g1, legend="none"))
g

#+ save heatmap figure
factor <- 0.6
ggsave(outPng, g, width = 17*factor, height = 16*factor)
ggsave(outPdf, g, width = 17*factor, height = 16*factor)

