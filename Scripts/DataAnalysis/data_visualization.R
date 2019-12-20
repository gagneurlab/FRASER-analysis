#'---
#' title: Dataset-wise splicing correlation
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fds_raw: '`sm config["DATADIR"] + "/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'   - fds_fil: '`sm config["DATADIR"] + "/datasets/savedObjects/{dataset}/pajdBinomial_psiSite.h5"`'
#'  output:
#'   - wBhtml: "Output/html/DataAnalysis/data_viz/{dataset}_data_viz.html"
#'  type: noindex
#'---

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/DataAnalysis/splicing_correlation.R", dataset="example")

    fds_gtex_skin_ns <- loadFraseRDataSet("/s/project/gtex-processed/splicing_map", "gtex-skin-notsunexposedsuprapubic", TRUE)
    fds_kremer       <- loadFraseRDataSet("/s/project/fraser/analysis/datasets", "kremer-bader-et-al", TRUE)
    fds_prokisch_all <- loadFraseRDataSet("/s/project/fraser/analysis/datasets", "prokisch_batch5", TRUE)

}

#+ source main config
source("./src/r/config.R")

#+ input
fdsFile    <- snakemake@input$fds_fil
dataset    <- snakemake@wildcards$dataset
workingDir <- dirname(dirname(dirname(fdsFile)))

#' # Load dataset
dataset
workingDir
fds_raw <- loadFraseRDataSet(workingDir, paste0("raw-", dataset))
fds     <- loadFraseRDataSet(workingDir, dataset)

#'
#' ## Number of samples
#'
ncol(fds)

#'
#' ## Number of junctions
#'
junction_numbers <- c(
    raw_junc = nrow(fds_raw),
    fil_junc = nrow(fds),
    raw_ss   = nrow(nonSplicedReads(fds_raw)),
    fil_ss   = nrow(nonSplicedReads(fds)))
data <- data.table(type=names(junction_numbers), count=junction_numbers)
ggplot(data, aes(type, count, fill=type)) +
    geom_bar(stat="identity") +
    scale_y_log10()

#'
#' ## PSI distribution
#'
print_freq <- 0.01
data <- rbindlist(lapply(psiTypes, function(i){
    n <- prod(dim(K(fds, type=i)))
    probs <- c(print_freq, 1-print_freq)
    selection <- sample(c(TRUE, FALSE), n, replace=TRUE, prob=probs)
    data.table(type=i, logitPsi=qlogis(assay(fds, i)[selection])) } ))
range <- range(data[is.finite(logitPsi),logitPsi]) + c(-1,1)
data[is.infinite(logitPsi),logitPsi:=ifelse(logitPsi < 0, range[1], range[2])]
ggplot(data, aes(type, logitPsi, fill=type)) + geom_violin()


#'
#' ## Variance across junctions
#'
data <- rbindlist(lapply(psiTypes, function(i){
    data.table(type=i, logitPsiSD=colSds(x(fds, type=i, all=TRUE))) } ))
ggplot(data, aes(type, logitPsiSD, fill=type)) + geom_violin()

#'
#' ## Coverage across junctions
#'
data <- rbindlist(lapply(psiTypes, function(i){
    data.table(type=i, meanK=rowMeans(K(fds, type=i))) } ))
ggplot(data, aes(type, meanK + 1, fill=type)) + geom_violin() +
    scale_y_log10()

data <- rbindlist(lapply(psiTypes, function(i){
    data.table(type=i, meanTotal=rowMeans(N(fds, type=i))) } ))
ggplot(data, aes(type, meanTotal + 1, fill=type)) + geom_violin() +
    scale_y_log10()

#'
#' Negative binomial fit of junction coverage
#' 
fitNegBinom <- function(par, counts){
  return( -mean(dnbinom(counts, mu=par[1], size=par[2], log=TRUE)) )
}
negBinomFit <- rbindlist(lapply(psiTypes, function(i){
  fit <- optim(c(500, 0.5), fitNegBinom, counts=round(data[type == i, meanTotal]))
  return( data.table(type=i, mu=fit$par[1], size=fit$par[2]) ) } )) 
DT::datatable(negBinomFit)

#'
#' ## Correlation
#'
known_factors_row <- c("condition")
known_factors_col <- c("FIBROBLAST_ID", "SEX", "BATCH", "TISSUE",
        "GROWTH_MEDIUM", "RNA_HOX_GROUP", "RNA_BATCH_GROUP",
        "SMCENTER", "SMRIN", "AGE", "DTHHRDY", "GENDER", "SMATSSCR")

known_factors_col <- known_factors_col[known_factors_col %in% colnames(colData(fds))]
known_factors_row <- known_factors_row[known_factors_row %in% colnames(colData(fds))]

if(length(known_factors_col) == 0){
    known_factors_col <- NA
}
if(length(known_factors_row) == 0){
    known_factors_row <- NA
}

plist <- lapply(psiTypes, plotCountCorHeatmap, fds=fds, logit=TRUE, topN=100000,
        annotation_col=known_factors_col, annotation_row=known_factors_row)

#'
#' ## Junctions grouped by shared donors/acceptors 
#'
dt <- data.table(
    chr = as.factor(seqnames(fds)),
    start = start(fds),
    end = end(fds),
    strand = as.factor(strand(fds)) ) 
groups <- rbind(data.table(groupsize = dt[,length(end), by=c("chr", "start", "strand")]$V1, groupedBy="donor"), 
                data.table(groupsize = dt[,length(start), by=c("chr", "end", "strand")]$V1, groupedBy="acceptor") )

summary(groups[groupedBy=="donor",groupsize])
summary(groups[groupedBy=="acceptor",groupsize])

ggplot(groups, aes(groupsize, fill=groupedBy)) + geom_histogram(alpha=0.7, position="identity", binwidth=1) #+ scale_y_log10()

#'
#' ## Distribution of H, D, rho after autoencoder fit
#'
rho <- rbindlist(lapply(psiTypes, function(i){
  data.table(rho=rho(fds, i), type=i) } ))
ggplot(rho, aes(x=rho, fill=type)) + geom_density(alpha = 0.7) + scale_x_log10()

E <- rbindlist(lapply(psiTypes, function(i){
  data.table(E=as.vector(E(fds, i)), type=i) } ))
ggplot(E, aes(x=E, fill=type)) + geom_density(alpha = 0.7) 

D <- rbindlist(lapply(psiTypes, function(i){
  data.table(D=as.vector(D(fds, i)), type=i) } ))
ggplot(D, aes(x=D, fill=type)) + geom_density(alpha = 0.7)

H <- rbindlist(lapply(psiTypes, function(i){
  data.table(H=as.vector(H(fds, i)), type=i) } ))
ggplot(H, aes(x=H, fill=type)) + geom_density(alpha = 0.7)



