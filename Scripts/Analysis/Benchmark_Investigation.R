#'---
#' title: Benchmark Investigation
#' author: Christian Mertes
#' wb:
#'  input: 
#'   - a: "wbuild.yaml"
#' output:
#'  html_document
#'---

#+ load config
source('./src/r/config.R')
register(MulticoreParam(10))

#' 
#' # Data import
#' 
#' * Skin not sun exposed psi5
#' 
#' ## FRASER object
fds <- loadFraseRDataSet(file="../FraseR-analysis/Data/paperPipeline/datasets/inject_uniformDistr/psi5/PCA-BB-Decoder/savedObjects/Skin_Not_Sun_Exposed_Suprapubic/fds-object.RDS")

#' ## Injection precision recall object
prefix <- "../FraseR-analysis/Data/paperPipeline/processedData/precRec/Skin_Not_Sun_Exposed_Suprapubic/inject_uniformDistr/psi5/"
dt <- fread(paste0(prefix, "PCA-BB-Decoder_plotData_byJunctionGroup_8.tsv.gz"))

#' adapt data.table to factors
dt[,Nsamples:=factor(Nsamples)]
dt[,geneMeanBin:=factor(geneMeanBin)]
dt[,inj_value:=factor(inj_value)]
dt[,dPsiBin:=factor(dPsiBin)]
dt[,correction:=factor(correction)]
dt[,dataset:=factor(dataset)]
dt[,pmethod:=factor(pmethod)]
dt[,junctionMeanBin:=factor(junctionMeanBin)]
dt[,inj:=factor(inj)]

sort(sapply(dt, class))
dt[1:50]

#' 
#' ## Numbers
#' Dimensions
dim(K(fds, "psi5"))


#' Number of alternative junctions expressed in dataset
table(data.table(index=getSiteIndex(fds, "psi5"))[,.N,by=index][,N])


#' Global Injections
otable <-table(assay(fds, "trueOutliers_psi5"))
otable
signif(otable/sum(otable), 3)

#' Expression mean
junc_rmea <- rowMeans2(N(fds, "psi5"))
junc_rmed <- rowMedians(N(fds, "psi5"))

par(mfrow=c(2,1))

hist(log10(junc_rmea), breaks=50, main="Hist of mean log10", xlim=c(-1,6))
abline(v=log10(quantile(junc_rmea, c(1:3/3))), col="firebrick")

hist(log10(junc_rmed), breaks=50, main="Hist of median log10", xlim=c(-1,6))
abline(v=log10(quantile(junc_rmed, c(1:3/3))), col="firebrick")

par(mfrow=c(1,1))

#' # Precision recall plots
#' 
#' 
#' ## z score only
dt <- dt[order(-abs(zScore))]
dt[,zrank:=.I]
dt[,zrecall:=cumsum(trueOutlier != 0)/sum(trueOutlier != 0)]
dt[,zprecision:=cumsum(trueOutlier != 0)/.I]
dt[,zcutoff:=cut(abs(zScore), c(0,2,3,Inf), right=FALSE)]
z14p <- unique(dt[,.(rec=round(zrecall, 5), prec=round(zprecision, 5), cutoff=zcutoff, score=round(abs(zScore), 2), type="z score only")])


#' ## z score filtered (dPsi, n, #zeros)
dt <- dt[,zScore2:=ifelse(n <= 5 | abs(dPsi) < 0.1 | nZeros > 60, 0, zScore)]
dt <- dt[order(-abs(zScore2))]
dt[,z2rank:=.I]
dt[,z2recall:=cumsum(trueOutlier != 0)/sum(trueOutlier != 0)]
dt[,z2precision:=cumsum(trueOutlier != 0)/.I]
dt[,z2cutoff:=cut(abs(zScore2), c(0,2,3,Inf), right=FALSE)]
z24p <- unique(dt[,.(rec=round(z2recall, 5), prec=round(z2precision, 5), cutoff=z2cutoff, score=round(abs(zScore2), 2), type="z score + n<=5 + dpis<0.1 + nzero>60")])


#' ## p values only
dt <- dt[order(pvalues)]
dt[,p1rank:=.I]
dt[,p1recall:=cumsum(trueOutlier != 0)/sum(trueOutlier != 0)]
dt[,p1precision:=cumsum(trueOutlier != 0)/.I]
dt[,p1cutoff:=cut(padjust, c(0,0.05,0.1,2), right=FALSE)]
p14p <- unique(dt[,.(rec=round(p1recall, 5), prec=round(p1precision, 5), cutoff=p1cutoff, score=round(-log10(pvalues), 2), type="p value only")])


#' ## p values filtered (dPsi, n, #zeros)
dt <- dt[,pvalues2:=ifelse(n <= 5 | abs(dPsi) < 0.1 | nZeros > 60, 1, pvalues)]
dt[,p2rank:=.I]
dt[,p2recall:=cumsum(trueOutlier != 0)/sum(trueOutlier != 0)]
dt[,p2precision:=cumsum(trueOutlier != 0)/.I]
dt[,p2cutoff:=cut(ifelse(pvalues2 == 1, 1, padjust), c(0,0.05,0.1,2), right=FALSE)]
p24p <- unique(dt[,.(rec=round(p2recall, 5), prec=round(p2precision, 5), cutoff=p2cutoff, score=round(-log10(pvalues2), 2), type="p value + n<=5 + dpis<0.1 + nzero>60")])


#' ## p value filtered on alternative junctions (#alternatives, dPsi, n, #zeros)
dts <- dt[nAlternatives > 1]
dts[,pvalues3:=ifelse(n <= 5 | abs(dPsi) < 0.1 | nZeros > 60, 1, pvalues)]
dts[,p3rank:=.I]
dts[,p3recall:=cumsum(trueOutlier != 0)/sum(trueOutlier != 0)]
dts[,p3precision:=cumsum(trueOutlier != 0)/.I]
dts[,p3cutoff:=cut(ifelse(pvalues3 == 1, 1, padjust), c(0,0.05,0.1,2), right=FALSE)]
p34p <- unique(dts[,.(rec=round(p3recall, 5), prec=round(p3precision, 5), cutoff=p3cutoff, score=round(-log10(pvalues3), 2), type="p value + nAlt>1 + n<=5 + dpis<0.1 + nzero>60")])



#' 
#' # Full precision-recall plot
#' 
dt4plotting <- rbindlist(list(z14p, z24p, p14p, p24p, p34p))
dt4plotting
ggplot(dt4plotting, aes(rec, prec, color=type, linetype=cutoff)) + geom_line() +
    scale_linetype_manual(values=c(3,2,1,1,2,3))

#'
#' # Distributions of values between precision cutoff
#' 
ggplot(dts, aes(log10(n+1), fill=p3precision > 0.5)) + 
    geom_histogram(bins=100) + 
    scale_y_log10()


ggplot(dts, aes(log10(geneMean), fill=p3precision > 0.5)) + 
    geom_histogram(bins=100) + 
    scale_y_log10()


ggplot(dts, aes(nZeros, fill=p3precision > 0.5)) + 
    geom_histogram(bins=100) + 
    scale_y_log10()


ggplot(dts, aes(abs(dPsi), fill=p3precision > 0.)) + 
    geom_histogram(bins=100) +
    scale_y_log10()


ggplot(dt, aes(abs(dPsi), fill=p2precision > 0.5)) + 
    geom_histogram(bins=100) +
    scale_y_log10()

if(FALSE){
    
    dttmp4p <- unique(dttmp)
    dttmp4p
    
    ggplot(dttmp4p, aes(rec, prec, color=value)) + geom_line() +
        scale_linetype_manual(values=c(1,2,3))
}
