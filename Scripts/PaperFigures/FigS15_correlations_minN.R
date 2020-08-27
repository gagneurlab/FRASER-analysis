#'---
#' title: Paper figure S15 (correlations with minN)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - stats:  '`sm expand(config["DATADIR"] + "/processedData/results/{dataset}/" + config["AE_IMPLEMENTATION"] + "_stats_minN.RDS", dataset=config["heatmap_tissues"])`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS15_correlations_minN.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS15_correlations_minN.pdf"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/PaperFigures/heatmaps/correlationFigure_minN.html"`'
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
    parseWBHeader2("Scripts/PaperFigures/FigS15_correlations_minN.R",
                   wildcards=wildcards, options=opt, rerun=TRUE)
}

#+ input
statFiles   <- snakemake@input$stats

#+ output
outPng      <- snakemake@output$outPng
outPdf      <- snakemake@output$outPdf

#'
#' # Load datasets
#+ echo=TRUE
outPng

#+ echo=FALSE
stats_ls <- lapply(statFiles, readRDS)
names(stats_ls) <- basename(dirname(statFiles))



#'
#' # Correlation reduction
#'
corDT <- rbindlist(lapply(stats_ls, "[[", "SampleCors_minN"))
levels <- corDT[,sort(unique(dataset))]
corDT[,dataset:=factor(dataset, levels=levels)]
levels(corDT$dataset) <- dName4plot(levels(corDT$dataset))
corDT[,method:=factor(method)]
corDT[,normalized:=factor(normalized, levels=c("raw", "normalized"))]
corDT[,minN:=factor(minN)]

par(mfrow = c(3, 3))
for(ptype in psiTypes){
    for(dset in corDT[,unique(dataset)]){
        LSD::heatscatter(corDT[type == ptype & dataset == dset & minN == 1, cor], 
                         corDT[type == ptype & dataset == dset & minN == 100, cor], 
                         xlim=c(-1, 1), ylim=c(-1, 1), main=paste0(dset, ": ", ptype),
                         xlab="sample correlations for minN=1", 
                         ylab="sample correlations for minN=100")
        abline(0, 1, lty="dashed")
    }
}
par(mfrow = c(1, 1))

corDT[,type:=factor(type, levels=c("psi5", "psi3", "psiSite"))]
levels(corDT$type) <- c("psi[5]", "psi[3]", "theta")
means <- corDT[,mean(nrJunctions),by="dataset,type,minN,normalized"]
means

corBoxplots <- ggplot(corDT, aes(y=abs(cor), x=minN, fill=normalized)) +
    facet_grid(type ~ dataset, labeller=labeller(type = label_parsed)) +
    geom_boxplot() +
    geom_text(data=means, aes(label=round(V1), y=0.95), size=3.4) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle=35, hjust=1),
        axis.text   = element_text(size=11),
        axis.title  = element_text(face="bold", size=14),
        legend.text = element_text(size=14) ) +
    ggtitle("") +
    ylab("|Sample correlation|") +
    xlab("Minimal total coverage") +
    scale_fill_brewer(palette = "Dark2", direction = -1,
                      name="Controlled", labels=c("No", "Yes"))
corBoxplots


#
# compare rand index between the clusterings
#
randIndex <- function(cluster1, cluster2){
    n <- length(cluster1)
    inSameClusterInBoth <- function(i,j, clu1, clu2){
        return( (clu1[i] == clu1[j]) == (clu2[i] == clu2[j]) )
    }
    sameClusterFunc <- Vectorize(inSameClusterInBoth, 
                                 vectorize.args=list("i","j"))
    clusterRes <- outer(1:n, 1:n, sameClusterFunc, 
                        clu1=cluster1, clu2=cluster2)
    
    ab <- sum(clusterRes[upper.tri(clusterRes)])
    rI <- (ab)/(choose(n, 2))
    
    return(rI)
}
getRandIndex <- function(dt, ds, ptype, minN1=0, minN2=100, nClust=5){
    dt1 <- dt[dataset == ds & type == ptype & minN == minN1,]
    dt2 <- dt[dataset == ds & type == ptype & minN == minN2,]
    
    n <- (1 + sqrt(1 + 8*nrow(dt1)) ) / 2
    cor1 <- matrix(1, nrow=n, ncol=n)
    cor1[upper.tri(cor1)] <- dt1[,cor]
    cor1[lower.tri(cor1)] <- t(cor1)[lower.tri(t(cor1))]
    cor2 <- matrix(1, nrow=n, ncol=n)
    cor2[upper.tri(cor2)] <- dt2[,cor]
    cor2[lower.tri(cor2)] <- t(cor2)[lower.tri(t(cor2))]
    
    c1 <- cutree(hclust(dist(cor1)), k=nClust)
    c2 <- cutree(hclust(dist(cor2)), k=nClust)
    
    rI <- randIndex(c1, c2)
    return(rI)
}

corDT_raw <- corDT[normalized == "raw",]
datasets <- corDT[,unique(dataset)]

#' compare clustering between min coverage 0 and 10
rI_0_10 <- sapply(datasets, function(ds, cor_dt){
    sapply(c("psi[3]", "psi[5]", "theta"), getRandIndex, dt=cor_dt, ds=ds, 
           minN1=0, minN2=10)
}, cor_dt=corDT_raw )
rI_0_10

#' compare clustering between min coverage 0 and 100
rI_0_100 <- sapply(datasets, function(ds, cor_dt){
    sapply(c("psi[3]", "psi[5]", "theta"), getRandIndex, dt=cor_dt, ds=ds, 
            minN1=0, minN2=100)
}, cor_dt=corDT_raw )
rI_0_100

#' compare clustering between min coverage 10 and 100
rI_10_100 <- sapply(datasets, function(ds, cor_dt){
    sapply(c("psi[3]", "psi[5]", "theta"), getRandIndex, dt=cor_dt, ds=ds, 
           minN1=10, minN2=100)
}, cor_dt=corDT_raw )
rI_10_100


#+ save figure
size  <- 9
ratio <- 1.05
outPng
ggsave(outPng, corBoxplots, width = size * ratio, height = size * 1/ratio)
ggsave(outPdf, corBoxplots, width = size * ratio, height = size * 1/ratio)
