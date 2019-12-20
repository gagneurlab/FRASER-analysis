#'---
#' title: Paper numbers
#' author: Christian Mertes
#' wb:
#'  input:
#'   - annoGTEx:    '`sm expand(config["DATADIR"] + "/annotations/{dataset}.tsv", dataset=config["EnrichmentTissues"])`'
#'   - statsGTEx:   '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_stats.RDS", dataset=config["EnrichmentTissues"])`'
#'   - annoKremer:  '`sm config["DATADIR"] + "/annotations/Kremer.tsv"`'
#'   - statsKremer: '`sm config["DATADIR"] + "/processedData/results/Kremer/" + config["AE_IMPLEMENTATION"] + "_stats.RDS"`'
#'   - resKremer:   '`sm config["DATADIR"] + "/processedData/results/Kremer/" + config["AE_IMPLEMENTATION"] + "_results.tsv"`'
#'   - fdsKremer:   '`sm config["DATADIR"] + "/datasets/savedObjects/Kremer__" + config["AE_IMPLEMENTATION"] + "/pvaluesBetaBinomial_psiSite.h5"`'
#'  threads: 5
#'---

#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=10, fig.height=10)

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list()
    opt <- c("--configfile", "wbuild_small.yaml")
    parseWBHeader2("Scripts/PaperFigures/PaperNumbers.R",
            wildcards=wildcards, options=opt, rerun=TRUE)
}

#+ input
statFiles        <- snakemake@input$statsGTEx
annoFiles        <- snakemake@input$annoGTEx
statFileKremer   <- snakemake@input$statsKremer
annoFileKremer   <- snakemake@input$annoKremer
resFileKremer    <- snakemake@input$resKremer
datasets         <- gsub(".tsv", "", basename(annoFiles))
names(annoFiles) <- datasets
names(statFiles) <- datasets

#'
#' # Load datasets
#+ echo=TRUE
datasets
annos  <- lapply(annoFiles, fread)
stats  <- lapply(statFiles, readRDS)
statsK <- readRDS(statFileKremer)
annoK  <- fread(annoFileKremer)
resK   <- fread(resFileKremer)
fdsK   <- loadFraseRDataSet(file=snakemake@input$fdsKremer)

#' # Numbers
#'
#' ## Number of analized samples and tissues
#'
statsdt <- as.data.table(pivot_longer(
    as.data.table(sapply(stats, "[[", "ImportantNumbers"), keep.rownames=TRUE),
    -rn, names_to="tissue", values_to="value"))
ggplot(statsdt, aes(x=tissue, y=value, fill=tissue)) +
    geom_bar(stat="identity") +
    facet_wrap(~rn, scales="free") +
    theme(axis.text.x=element_blank(), legend.position="bottom")

#' * Number of tissues
length(unique(statsdt[,tissue]))

#' * Number of total sample across tissues
sum(statsdt[rn=="Nsamples", value])

#' * Number of individuals analyzed
length(unique(rbindlist(annos)[,indivID]))

#' * Number of non-sun-exposed skin samples
nrow(annos[["Skin_Not_Sun_Exposed_Suprapubic"]])

#' 
#' * Number of donors and acceptors and percentage of noval ones
#' 
getDNAPercentage <- function(dt, type){
    ans <- sapply(unique(dt$tissue), function(x){
            known <- dt[rn==paste0("known", type) & tissue == x, value]
            total <- dt[rn==paste0("N", type)     & tissue == x, value] 
            c(percentNoval=100 * (1 - known/total), noval=total - known, 
                    total=total)})
    list(all=sort(ans["percentNoval",]), mean=mean(ans["percentNoval",]), 
            sd=sd(ans["percentNoval",]), novelSites=ans["noval",], 
            totalMean=mean(ans["total",]), totalSD=sd(ans["total",]),
            novalSitesMean=mean(ans["noval",]), novalSitesSD=sd(ans["noval",]))
}
getDNAPercentage(statsdt, "Starts")
getDNAPercentage(statsdt, "Ends")

#' 
#' * Optimal q's + sd
statsdt[rn %in% paste0("Q_", psiTypes)][,.(
        mean=round(mean(value), 2), sd=round(sd(value),2)), by=rn]

#'
#' ## Sample correlations
#'
rbindlist(lapply(stats, "[[", "SampleCors"))[,.(
        mean=mean(abs(cor)),
        sd=sd(abs(cor))), by="dataset,method,type,normalized"][,.(
                mean=round(mean(mean, na.rm=TRUE), 3),
                sd=round(sd(mean, na.rm=TRUE), 3)), by="method,type,normalized"]

#'
#' ## Kremer data set and results
#'
statsK$ImportantNumbers[c("NJunctions","NSites", "NRawJunctions", "NRawSites")]
resK <- results(fdsK, padjCutoff=0.05, zScoreCutoff=NA, deltaPsiCutoff=0.3)

#'
#' ### TAZ diagnostics
#' 
sampleID <- "MUC1398"
geneID   <- "TAZ"
plotVolcano(fdsK, sampleID, "psi5", aggregate=TRUE, basePlot=FALSE)

#'
#' ## GTEx data set summary
#'
rowMeans(sapply(stats, "[[", "ImportantNumbers"))[c("NJunctions","NSites", "NRawJunctions", "NRawSites")]
apply(sapply(stats, "[[", "ImportantNumbers"), 1, sd)[c("NJunctions","NSites", "NRawJunctions", "NRawSites")]




