library(data.table)
library(ggplot2)
library(parallel)
library(cowplot)
library(ggpubr)

#' 
#' Bootstrap code from here:
#' https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf
#' 

predictLabels <- function(n=1e4, name='example'){
  score <- abs(rnorm(n))
  p <- 1-rank(score)/(length(score))
  label <- sapply(p, function(px) {
    px <- pmin(0.85, pmax(0.15, px))
    sample(c(T, F), 1, prob=c(1-px, px)) })
  data.table(score, label, name=as.character(name))
}


calculatePrecisionRecall <- function(score, label, decr=TRUE, total=sum(label)){
  dt <- data.table(score=score, label=label)
  dt <- dt[order(score, decreasing=decr)]
  dt[,rank      := 1:.N]
  dt[,TP        := cumsum(label)]
  dt[,c('TP', 'rank'):=list(max(TP), max(rank)), by=score]
  dt[,precision := TP/rank]
  dt[,recall    := TP/total]
  dt
}

bootstrapPR <- function(score, label, total=sum(label), decr=TRUE, n=500, mc.cores=10){
  stopifnot(length(score) == length(label))
  if(length(score) == 0){
    warning('score is empty')
    return(data.table(
      score=NA_real_, label=FALSE, rank=NA_integer_, TP=NA_integer_,
      precision=NA_real_, recall=NA_real_, lower=NA_real_, upper=NA_real_))
  }
  
  d <- mclapply(0:n, mc.cores=mc.cores, function(i){
    if(i == 0){
      mySampleId <- seq_along(score)
    } else {
      mySampleId <- sample(length(score), replace=TRUE)
    }
    d <- calculatePrecisionRecall(score=score[mySampleId], 
        label=label[mySampleId], decr=decr, total=sum(label[mySampleId]))
    d[,bootstrep:=i]
    d
  })
  td <- d[[1]]
  d1 <- rbindlist(d[1:n + 1])
  
  d1[,tPrecision:=rep(td$precision,n)]
  d1[,pos:=rep(seq_row(td), n)]
  
  # estimating the quantiles (just taking the nth element)
  d1 <- d1[order(rank, precision - tPrecision)]
  lowerDiff <- d1[, .(lDiff=(precision - tPrecision)[max(1, floor(length(precision)*0.025))]), by=rank]
  upperDiff <- d1[, .(uDiff=(precision - tPrecision)[ceiling(length(precision)*0.975)]), by=rank]
  
  ans <- merge(merge(td, lowerDiff, all.x=TRUE), upperDiff, all.x=TRUE)
  
  ans[, lower:=precision + lDiff]
  ans[, upper:=precision + uDiff]
  
  ans
}

if(FALSE){
  
  data <- rbindlist(lapply(1:3, function(x) predictLabels(n=1e3, name=x)))
  data
  ggplot(data=data, aes(x=label, y=score, col=name)) + 
    geom_boxplot(position=position_dodge(0.85)) + 
    geom_violin(position=position_dodge(0.85))
  
  plotdt <- data[, bootstrapPR(score=score, label=label, TRUE), by=name]
  p1 <- ggplot(plotdt, aes(recall, precision, col=name)) + geom_line() +
    grids() + ylim(c(0,1)) + xlim(c(0,1))
  p1
  p2 <- ggplot(plotdt, aes(recall, precision, col=name)) + geom_line() +
    geom_ribbon(aes(x=recall, 
                    ymin = pmax(0, pmin(1, lower)), 
                    ymax = pmax(0, pmin(1, upper)), fill=name), 
                alpha=0.2, col=NA) +
    grids() + ylim(c(0,1)) + xlim(c(0,1))
  p2
  plot_grid(p1, p2)
  
  
  library(microbenchmark)
  
  scores <- data.table(s=rnorm(1e6))
  f1 <- function(){
    scores <- scores[order(s)]
    scores[,rank:=1:.N]
    scores[,rank:=max(rank), by=s]
  }
  
  f2 <- function(){
    scores <- scores[order(s)]
    scores[,rank:=1:.N]
    #scores[,rank:=max(rank), by=s]
  }
  
  f3 <- function() { rank(scores$s) }
  
  f4 <- function() { rank(scores$s, ties.method = 'max') }
  
  
  res <- microbenchmark(
    data.table.max=f1(), 
    data.table.random=f2(), 
    rank.default=f3(), 
    rank.max=f4(),
    times=20)
  res
  autoplot(res)
}
