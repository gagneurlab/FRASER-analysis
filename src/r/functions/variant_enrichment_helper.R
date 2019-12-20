load_fraser_enrichment_data <- function(file, internCPUs=3, minCoverage=10,
                    minDeltaPsi=0.1, debug=FALSE){
    curTissue <- basename(strsplit(file, "__")[[1]][1])
    curMethod <- dirname(strsplit(file, "__")[[1]][2])
    
    # get results
    fds <- loadFraseRDataSet(file=file)
    if(debug == TRUE){
        fds <- fds[seqnames(fds) %in% c(1, 2, "chr1", "chr2"), 1:50]
    }
    fds <- annotateRanges(fds)
    mc.cores <- 3
    resls <- mclapply(psiTypes, extractCalls4Enrichment, 
            mc.cores=mc.cores, Ncpus=internCPUs, fds=fds, 
            minCoverage=minCoverage, minDeltaPsi=minDeltaPsi)
    
    res <- rbindlist(resls)[!is.na(hgnc_symbol),.(
        pvalue=min(pvalue),
        pvalue_norm=min(pvalue_norm),
        zscore=zscore[abs(zscore) == max(abs(zscore))][1]),
        by="hgnc_symbol,sampleID"]
    
    ans <- res[,.(subjectID=sampleID, geneID=hgnc_symbol, p=pvalue,
                  z=zscore, p_norm=pvalue_norm, tissue=curTissue)][order(p)]
    
    names(ans) <- c("subjectID", "geneID",
                    paste0(curMethod, c("_p", "_z", "_np")), "tissue")
    
    if(isTRUE(debug)){
        ans <- list(enrich_obj=ans[,-"tissue"], fds=fds, res=res[pvalue<0.1])
    } else {
        ans <- list(enrich_obj=ans[,-"tissue"])
    }
    ans
}


extractCalls4Enrichment <- function(fds, type, minCoverage=10, minDeltaPsi=0.1,
                    Ncpus=1){
    BPPARAM <- SerialParam()
    if(Ncpus > 1){
        BPPARAM <- MulticoreParam(Ncpus)
    }
    
    #' 
    #' get data
    #' 
    index  <- getSiteIndex(fds, type)
    symbol <- mcols(fds, type=type)[,"hgnc_symbol"]
    pv     <- as.matrix(pVals(fds, type))
    dp     <- as.matrix(getDeltaPsi(fds, type))
    zs     <- as.matrix(zScores(fds, type))
    n      <- as.matrix(N(fds, type))
    
    #
    # z score extraction
    # take absolute maximum per gene
    # 
    zs[abs(dp) < minDeltaPsi | n < minCoverage] <- 0
    zs[is.na(zs)] <- 0
    zsBG <- getAbsMaxByGroup(mat=zs, index=symbol, BPPARAM=BPPARAM)
    
    # 
    # p value (betaBinomial) extraction
    # correct with FWER per gene and take minimum
    # 
    pv[abs(dp) < minDeltaPsi | n < minCoverage] <- 1
    pv <- pmin(pmax(pv, 1e-250), 1)
    pvBS <- 10^getAbsMaxByGroup(mat=log10(pv), index=index, BPPARAM=BPPARAM)
    symbolIndex <- symbol[!duplicated(index)]
    pvBG_ls <- bplapply(seq_col(pvBS), adjust_FWER_PValues,
            pvals=pvBS, index=symbolIndex, BPPARAM=BPPARAM)
    pvBG <- do.call(cbind, pvBG_ls)[!duplicated(symbolIndex),]
    rownames(pvBG) <- symbolIndex[!duplicated(symbolIndex)]
    
    # 
    # p value (normal) extraction
    # correct with FWER per gene and take minimum
    # 
    if(!paste0("pvaluesNormal_", type) %in% assayNames(fds)){
        pnvBG <- matrix(1, nrow=sum(!duplicated(symbolIndex)))
        rownames(pnvBG) <- symbolIndex[!duplicated(symbolIndex)]
    } else {
        pnv <- as.matrix(pVals(fds, type, dist="Normal"))
        pnv[abs(dp) < minDeltaPsi | n < minCoverage] <- 1
        pnv <- pmin(pmax(pnv, 1e-250), 1)
        pnvBS <- 10^getAbsMaxByGroup(mat=log10(pnv), index=index, 
                BPPARAM=BPPARAM)
        symbolIndex <- symbol[!duplicated(index)]
        pnvBG_ls <- bplapply(seq_col(pnvBS), adjust_FWER_PValues,
                pvals=pnvBS, index=symbolIndex, BPPARAM=BPPARAM)
        pnvBG <- do.call(cbind, pnvBG_ls)[!duplicated(symbolIndex),]
        rownames(pnvBG) <- symbolIndex[!duplicated(symbolIndex)]
    }
    
    # 
    # merge final result table
    # 
    final_symbols <- na.omit(intersect(rownames(zsBG), rownames(pvBG)))
    resdt <- data.table(
        hgnc_symbol=final_symbols,
        sampleID=rep(colnames(zsBG), each=length(final_symbols)),
        pvalue=c(pvBG[final_symbols,]),
        pvalue_norm=c(pnvBG[final_symbols,]),
        zscore=c(zsBG[final_symbols,]))
    
    return(resdt)
}


calculateEnrichment <- function(dt, cols, cutOff, isZscore=TRUE, na.rm=FALSE){
    dt <- copy(dt)
    if(isTRUE(isZscore)){
        if(isTRUE(na.rm)){
            dt[is.na(get(cols)), c(cols):=0]
        }
        dt[,cutoff:=abs(get(cols)) > cutOff]
    } else {
        if(isTRUE(na.rm)){
            dt[is.na(get(cols)), c(cols):=1]
        }
        dt[,cutoff:=get(cols) < cutOff]
    }

    ans <- dt[, .(nRareEvent=sum(!is.na(simple_conseq)), 
            total=.N, nNA=sum(is.na(get(cols)))), by=cutoff]
    if(!any(ans$cutoff == FALSE, na.rm = TRUE)){
        ans <- rbind(ans, data.table(cutoff=FALSE, nRareEvent=0, total=0, nNA=0))
    }
    if(!any(ans$cutoff == TRUE, na.rm = TRUE)){
        ans <- rbind(ans, data.table(cutoff=TRUE, nRareEvent=0, total=0, nNA=0))
    }
    ans <- ans[,.(nRareEvent, total, fraction=(nRareEvent)/(total), nNA), by=cutoff]
    enrich <- ans[cutoff==TRUE,fraction]/ans[cutoff==FALSE,fraction]

    # get bound as done in Li et al
    out.var      <- ans[cutoff == TRUE,  nRareEvent]
    out.total    <- ans[cutoff == TRUE,  total]
    nonout.var   <- ans[cutoff == FALSE, nRareEvent]
    nonout.total <- ans[cutoff == FALSE, total]
    log.se = sqrt(1/out.var - 1/out.total + 1/nonout.var - 1/nonout.total)
    max.ci = enrich * exp(1.96*log.se)
    min.ci = enrich * exp(-1.96*log.se)

    return(list(dt=ans, enrichment=enrich, max.ci=max.ci, min.ci=min.ci))
}

getEnrichmentForMethod <- function(dt, method, cutoffs=c(0.005, 2)){
    message(method)
    isZscore <- grepl('_z$', method)
    curCutoff <- cutoffs[isZscore + 1]
    curType <- paste0(c('P-value (<', 'Z score (>'), cutoffs, ')')[isZscore + 1]
    ans <- calculateEnrichment(dt, method, curCutoff, isZscore)
    return(data.table(
            Method=gsub('_[pz]$', '', method), Type=curType, Cutoff=curCutoff,
            enrichment=ans$enrichment, max.ci=ans$max.ci, min.ci=ans$min.ci,
            out.var      = ans$dt[cutoff == TRUE,  nRareEvent],
            out.total    = ans$dt[cutoff == TRUE,  total],
            nonout.var   = ans$dt[cutoff == FALSE, nRareEvent],
            nonout.total = ans$dt[cutoff == FALSE, total]))
}

getEnrichmentForTissues <- function(tissue, rdsFiles, cutoffs){
    message(tissue)
    rds <- readRDS(rdsFiles[[tissue]])
    methods <- grep('_[pz]$', colnames(rds), value=TRUE)
    enrichments <- rbindlist(lapply(methods, getEnrichmentForMethod,
            dt=rds, cutoffs=cutoffs))
    enrichments[,Tissue:=tissue]
    return(enrichments)
}

plotEnrichment <- function(dt, numVarOffset=0.8){

    ggplot(dt, aes(x=Method, y=enrichment, col=Type)) +
        geom_point(position=position_dodge(.75)) +
        geom_errorbar(aes(ymin=min.ci, ymax=max.ci, col=Type), width=.2,
                position=position_dodge(.75)) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        geom_hline(yintercept=1, col='orange') +
        geom_point(aes(Method, numVarOffset, col=Type, shape=numVars),
                position=position_dodge(.75)) +
        grids() +
        facet_wrap('Tissue')
}

calculateRecallRank <- function(dt, cols, isPvalue=TRUE){
    dt <- dt[order(abs(get(cols)), decreasing=!isPvalue)]
    dt[,c(paste0(cols, '_recall')):=cumsum(!is.na(simple_conseq))]
    dt[,c(paste0(cols, '_rank')):=seq_len(.N)]

    return(dt)
}

simplifyConsequences <- function(vardt, splicingFirst=FALSE){
    SORanking <- c(
        # HIGH
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        # MODERAT
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        # LOW
        "splice_region_variant",
        "synonymous_variant")
    if(isTRUE(splicingFirst)){
        idx <- startsWith(SORanking, "splice_")
        SORanking <- c(SORanking[idx], SORanking[!idx])
    }

    vardt[,simple_conseq:=NA_character_]
    sapply(SORanking, function(x){
        vardt[is.na(simple_conseq) & grepl(x, Consequence),
                c("rank", "simple_conseq"):=list(which(SORanking==x), x)]
    })
    
    vardt[,simple_conseq:=factor(simple_conseq, 
            levels=unique(c(SORanking, simple_conseq)))]
    return(vardt)
}

plotRecallRankForEnrichment <- function(dt, maxRank, maxPoints, logy=FALSE, logx=FALSE){
    if(!missing(maxRank)){
        dt <- dt[rank<maxRank]
    }
    if(!missing(maxPoints)){
        prob4Samp <- min(1, maxPoints / (nrow(dt)/nrow(unique(dt[,.(Method, Type)]))))
        dt <- dt[
            rank < 1e3 |
            rank > max(rank) - 1e3 |
            sample(c(TRUE, FALSE), .N, replace = TRUE, prob = c(prob4Samp, 1 - prob4Samp))]
    }

    gg <- ggplot(dt[,.(rank=rank, recall=recall, Method, Type)],
            aes(rank, recall, col=Method, linetype=Type)) +
        geom_line() + grids()

    if(isTRUE(logx)){
        gg <- gg + scale_x_log10()
    }
    if(isTRUE(logy)){
        gg <- gg + scale_y_log10()
    }
    gg
}

maf2number <- function(maf, aggregateFun=max, naMAF=FALSE){
    mafls <- strsplit(maf, '&')
    ans <- unlist(bplapply(mafls, agrFun=aggregateFun,
        function(x, agrFun){ agrFun(as.numeric(gsub('[A-Z-]+:', '', x)))}))
    if(isScalarNumeric(naMAF)){
        ans[is.na(ans)] <- naMAF
    }
    return(ans)
}

readVcfParallel <- function(vcfFile, ranges, threads, nPerChunk=25,
                    vcfParam=ScanVcfParam()){
    rangesOfInt <- reduce(ranges)
    nchunks <- max(10, ceiling(length(rangesOfInt)/nPerChunk))
    BPPARAM <- MulticoreParam(threads, nchunks, progressbar=TRUE)
    chunks <- chunk(seq_along(rangesOfInt), n.chunks=nchunks)

    message(date(), ':',
            ' Run with number of chunks: ', length(chunks),
            ' and size of chunks: ', length(chunks[[1]]),
            ' with threads in parallel: ', threads)

    if(isFALSE(bpisup(BPPARAM))){
        bpstart(BPPARAM)
    }
    # read vcf file in chunks
    vcfgtls <- bplapply(chunks, f=vcfFile, ranges=rangesOfInt, param=vcfParam,
            BPPARAM=BPPARAM,
            function(x, ranges, f, param){
                    vcfWhich(vcfParam) <- ranges[x]
                    readVcf(TabixFile(f), param=vcfParam) })
    message(date(), ': done extracting genotypes')
    vcfgt <- do.call(rbind, vcfgtls)

    bpstop(BPPARAM)

    # remove unwanted variants lying in the same region
    vcfgt <- vcfgt[unique(names(ranges))]
    return(vcfgt)
}

