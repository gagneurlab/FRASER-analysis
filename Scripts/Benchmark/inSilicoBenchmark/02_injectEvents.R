#'---
#' title: Create Benchmarksets - Inject insilico events
#' author: Christian Mertes
#' wb:
#'  input:
#'   - datasetdone: 'Output/html/Benchmark/inSilicoBenchmark/inSilico_n{nsamples}_{dataset}.html'
#'  output:
#'   - countsJ: '`sm config["DATADIR"] + "/savedObjects/inSilico_efreq{efreq}_n{nsamples}_dp{dp}_{dataset}/rawCountsJ.h5"`'
#'   - countsS: '`sm config["DATADIR"] + "/savedObjects/inSilico_efreq{efreq}_n{nsamples}_dp{dp}_{dataset}/rawCountsSS.h5"`'
#'   - injections: 'Data/benchmark/insilicoInjection/efreq{efreq}_n{nsamples}_dp{dp}_{dataset}_event_table.tsv'
#'   - wBhtml: 'Output/html/Benchmark/inSilicoBenchmark/efreq{efreq}_n{nsamples}_dp{dp}_{dataset}_injections.html'
#'  type: noindex
#'---

#+ echo=FALSE
source("./src/r/config.R")

#'
#' # Dataset
#+ echo=TRUE
set.seed(42)
snakemake@wildcards$dataset
inpname <- gsub(".html$", "", basename(snakemake@input$datasetdone))
outname <- basename(dirname(snakemake@output$countsJ))
wddir   <- dirname(dirname(dirname(snakemake@output$countsJ)))
fds     <- loadFraseRDataSet(wddir, inpname)
ns      <- as.integer(snakemake@wildcards$nsamples)
ef      <- as.double(snakemake@wildcards$efreq)
dp      <- as.double(snakemake@wildcards$dp)
ne      <- ceiling(length(fds)*ef*ns/10)*10

#'
#' # Create injection set
#'
injections <- data.table(
    psiType=sample(c("psi5", "psi3", "psiSite"), ne, replace=TRUE),
    junctionID=sample(1:length(fds), ne),
    sampleIdx=sample(1:dim(fds)[2], ne, replace=TRUE),
    overExpression=sample(c(TRUE, FALSE), ne, replace=TRUE),
    psiSiteType=NA_character_
)
injections[psiType=="psiSite", psiSiteType:=sample(c("5", "3"), .N, replace=TRUE)]

inj <- data.table()
for(i in c("psi5", "psi3", "psiSite")){
    newgr <- as.data.table(granges(fds)[injections[psiType==i,junctionID]])
    newgr1 <- newgr[, .(seqnames, start, end, width, strand, startID, endID,
                        hgnc_symbol)]
    if(i != "psiSite"){
        newgr1[,c("alpha", "beta"):=list(
            alpha=newgr[,get(paste0(i, "_alpha"))],
            beta =newgr[,get(paste0(i, "_beta" ))])]
    }
    newgr1 <- cbind(newgr1, injections[psiType==i])
    inj <- rbind(inj, newgr1, fill=TRUE)
}
inj[,sampleID:=samples(fds)[sampleIdx]]

for(i in 1:nrow(inj)){
    if(inj[i,psiType!="psiSite"]){
        next
    }
    spliceID <- inj[i, ifelse(psiSiteType == "5", startID, endID)]
    idx <- mcols(nonSplicedReads(fds))$spliceSiteID == spliceID
    if(sum(idx) != 1){
        warning(paste(i, sum(idx), "not good idx"))
        next
    }
    dt <- as.data.table(mcols(nonSplicedReads(fds)[idx]))
    gr <- granges(nonSplicedReads(fds)[idx])
    inj[i,c("start", "end", "width", "alpha", "beta", "junctionID"):=list(
        start(nonSplicedReads(fds))[idx],
        end(nonSplicedReads(fds)[idx]),
        width(nonSplicedReads(fds)[idx]),
        alpha=dt$psiSite_alpha,
        beta=dt$psiSite_beta,
        junctionID=which(idx))]
}

#'
#' Remove NAs from the table
#'
inj <- inj[!is.na(alpha)]
inj <- inj[sample(1:nrow(inj), min(nrow(inj), ne/10))]
inj[,type:=psiType]
grinj <- makeGRangesFromDataFrame(inj, keep.extra.columns = TRUE)

#'
#' Write benchmark set
#+ write tsv
write_tsv(inj, snakemake@output$injections)
inj


#'
#' save only injected junctions
#+ save new fds
newfds <- fds[na2false(mcols(fds, type="j")$startID %in% inj$startID) &
        na2false(mcols(fds, type="j")$endID %in% inj$endID)]
name(newfds) <- outname
newfds <- saveFraseRDataSet(newfds)


#'
#' Inject into count table
#+ start injection
cts <- list(
    "rawJ"      = as(counts(newfds, type="psi5",    side="ofI"), "matrix"),
    "rawS"      = as(counts(newfds, type="psiSite", side="ofI"), "matrix"),
    "psi5_o"    = as(counts(newfds, type="psi5",    side="oth"), "matrix"),
    "psi3_o"    = as(counts(newfds, type="psi3",    side="oth"), "matrix"),
    "psiSite_o" = as(counts(newfds, type="psiSite", side="oth"), "matrix"))
injMeth <- "deltaPsi"
for(i in 1:nrow(inj)){
    type <- inj[i, psiType]
    sidx <- inj[i, sampleIdx]
    if(type == "psiSite"){
        jidx <- from(findOverlaps(nonSplicedReads(newfds), grinj[i], type = "equal"))
    } else {
        jidx <- from(findOverlaps(newfds, grinj[i], type = "equal"))
    }
    medianPsi <- rowMedians(as.matrix(assays(newfds)[[type]][jidx,]))

    size <- max(5,
            cts[[ifelse(type=="psiSite", "rawS", "rawJ")]][jidx,sidx] +
            cts[[paste0(type, "_o")]][jidx,sidx])
    k    <- bbmean(size, inj[i, alpha], inj[i, beta])
    kvar <- bbvariance(size, inj[i, alpha], inj[i, beta])


    if(injMeth == "deltaPsi"){
        newk <- size*(medianPsi-dp*ifelse(medianPsi<0.5,-1,1))
    } else if(injMeth == "sd"){
        newk <- k - kvar*sd
    }
    newk <- max(newk, 0)
    newo <- max(size - newk, 0)

    cts[[ifelse(type=="psiSite", "rawS", "rawJ")]][jidx,sidx] <- ceiling(newk)
    cts[[paste0(type, "_o")]][jidx,sidx] <- ceiling(newo)
}

counts(newfds, type="psi5",    side="ofI") <- cts[["rawJ"]]
counts(newfds, type="psiSite", side="ofI") <- cts[["rawS"]]
counts(newfds, type="psi5",    side="oth") <- cts[["psi5_o"]]
counts(newfds, type="psi3",    side="oth") <- cts[["psi3_o"]]
counts(newfds, type="psiSite", side="oth") <- cts[["psiSite_o"]]

parallel(newfds) <- MulticoreParam(30, 150, progressbar=TRUE)
register(parallel(newfds))

#'
#' remove old values
#+ start recalculating values
for(i in grep("^((delta|zscore|pvalue)_)?psi",
            assayNames(newfds), value=TRUE, perl=TRUE)){
    assays(newfds)[[i]] <- NULL
}

#' calculate PSI values
newfds <- calculatePSIValues(newfds, overwriteCts=FALSE)

#' calculate ZScores
newfds <- calculateZScores(newfds)

#' calculte P-values
newfds <- calculatePValues(newfds)

#' save final fds object
newfds <- saveFraseRDataSet(newfds)

oriFds <- fds
fds <- newfds
benchsetgr <- grinj
benchset <- inj
