#'---
#' title: ASHG 2019 Presentation figures
#' author: Christian Mertes
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
source("./src/r/config.R")


#'
#' # Correlation
#'
tissue <- "Cells_Transformed_fibroblasts"
method <- "PCA-BB-Decoder"
wdir <- "/s/project/fraser/snakemake_pipeline/Data/paperPipeline/datasets"

fds <- loadFraseRDataSet(wdir, paste0(tissue, "__", method))

annotation_row <- as.data.frame(colData(fds)[,c('SMCENTER'), drop=FALSE])
annotation_col <- as.data.frame(colData(fds)[,c('AGE', 'SMRIN'), drop=FALSE])
annotation_col$AGE <- as.integer(gsub("-.*", "", colData(fds)$AGE))

colnames(annotation_row) <- c("Center")
colnames(annotation_col) <- c("Age", "RIN")

#' Uncorrected data
plotCountCorHeatmap(fds, logit=TRUE, normalized=FALSE, sampleClustering=NA,
        annotation_row=annotation_row, annotation_col=annotation_col,
        )

#' Corrected data
plotCountCorHeatmap(fds, logit=TRUE, normalized=TRUE, sampleClustering=NA,
        annotation_row=annotation_row, annotation_col=annotation_col)



resDir <- "/data/ouga04b/ag_gagneur/project_local/fraser/snakemake_pipeline/Data/paperPipeline/processedData/results"
statsls <- lapply(
    list.files(resDir, pattern="PCA-BB.*.RDS", full.names=TRUE, recursive=TRUE),
    readRDS)

corDT <- rbindlist(lapply(statsls, "[[", "SampleCors"))[type == 'psi5']
corDT[,dataset:=factor(dataset)]
levels(corDT$dataset) <- toTitleCase(gsub("_", " ", levels(corDT$dataset)))
corDT[,normalized:=factor(normalized, levels=c("raw", "normalized"))]

# remove long name in the beginning
corDT <- corDT[!grepl("Gland|Omentum|Subcutaneous|leum|Junction|Lower", dataset)]

corBoxplots <- ggplot(corDT, aes(y=abs(cor), x=dataset, fill=normalized)) +
    geom_boxplot() +
    theme_bw() +
    theme(
        axis.text.x=element_text(angle = 20, hjust=1),
        axis.text = element_text( size = 14),
        axis.title = element_text(face="bold", size=14),
        legend.text = element_text(size=14) ) +
    ggtitle("") +
    ylab("|Sample correlation|") +
    xlab("GTEx tissue") +
    scale_fill_brewer(palette = "Dark2", direction = -1,
                      name="Controlled", labels=c("No", "Yes"))

#' export with 1100 x 350
corBoxplots


stats

type <- "psi5"
tissue <- "Cells_Transformed_fibroblasts"
method <- "PCA-BB-Decoder"
wdir <- "/s/project/fraser/snakemake_pipeline/Data/paperPipeline/datasets"

fds <- loadFraseRDataSet(wdir, paste0(tissue, "__", method))




