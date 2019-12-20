# TODO: Add comment
#
# Author: mertes
###############################################################################

# source main config
source("./src/r/config.R")
source("../mitomultiomics/src/r/functions/autoIGV.R")
load_all(PKG_ROOT)

#' get variants
vars <- readRDS(GTEX_FILTERED_VARIANTS_RDS)
gtexMatGT <- vars$gtexMatGT

#' get dataset
fds <- loadFraseRDataSet("/s/project/fraser/analysis/datasets", "gtex-fibroblast-transformed")
sampleData <- as.data.table(colData(fds))

#' get only sample of interest
gtexMatGT <- gtexMatGT[rownames(gtexMatGT) %in% sampleData[,condition],]

#' get only private events
gtexMatGT <- gtexMatGT[,apply(gtexMatGT, 2, function(x) sum(x > 0) == 1)]

#' get only samples of interest
samples2Viz <- sampleData[unique(c(which(grepl("SRR10790", sampleID)), 1:100))]
gt2viz <- gtexMatGT[rownames(gtexMatGT) %in% samples2Viz[,condition],]
gt2viz <- gt2viz[,apply(gt2viz, 2, function(x) sum(x > 0) == 1)]
dim(gt2viz)

#' get var string
vardt <- rbindlist(lapply(strsplit(colnames(gt2viz), "_"),
        function(x) as.data.table(t(data.table(x)))))

vardt[,condition:=apply(gt2viz, 2, function(x) names(which(x > 0)))]
vardt[,sampleID:=sapply(condition, function(x) sampleData[x==condition, sampleID])]
vardt[,V2:=as.integer(V2)]
coords <- vardt[, .(chr=paste0("chr", V1), start=V2-2, end=V2+2,
        maximum = -1, coverage.only = F,
        batch_file_prefix = paste("gtex-fib-trans", V1, V2, sep="_"),
        suffix.to.file = paste("gtex-fib-trans", V1, V2, sep="_"),
        sampleID=sampleID)]

# data to show
bamPaths <- samples2Viz[, .(paths=bamFile, names=sampleID, show_junctions=TRUE)]

# run autoIGV
res <- sapply(1:nrow(coords), coords=coords, bamPaths=bamPaths,
                function(idx, coords, bamPaths){
        bams <- unique(rbind(bamPaths[names == coords[idx, sampleID]], bamPaths[1:15]))[1:10]
        AutoIGV(coords[idx], indents=c(2500,2500), bams_paths=bams,
        folder_name="/s/project/fraser/analysis/datasets/benchmark/variant-igv-pics")
})
