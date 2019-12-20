varHandFile <- "../mitomultiomics/src/r/functions/variant_handling/omim_parser.R"
if(file.exists(varHandFile)){
    source(varHandFile)
}

#'
#' Add additional sample infos to data.table
#'
#' The final result table should look like:
#'   * sample identifiers: RNA_ID, FIBROBLAST_ID
#'   * FraseR restuls: gene, type, range, pvalue, fds, zscore, psivalue
#'   * Found by other algorithms: Leafcutter, JunctionsSeq, Cummings et al
#'   * Annotations: is mito vip gene, number of findings, MAE, variant, Abex
#'
addSampleNGeneInfosToResults <- function(res, fdrCutoff=0.1, addMissing=FALSE){
    if(!is(res, "data.table")){
        res <- as.data.table(res)
    }

    if(!"MITOMAP_DATATABLE" %in% ls(envir=.GlobalEnv)){
        eval(load_mitomap(), envir=globalenv())
    }

    # add RNA_ID if not there yet
    if(!"RNA_ID" %in% colnames(res)){
        res[,RNA_ID:=sampleID]
    }

    # is the patient solved already or has a known mutation?
    if(!"KNOWN_MUTATION" %in% colnames(res)){
        mut2merge <- unique(MITOMAP_DATATABLE[!is.na(RNA_ID),
                .(RNA_ID, KNOWN_MUTATION, SOLVED)])
        res <- merge(res, mut2merge, by="RNA_ID", all.x=TRUE)
    }


    # add fibroblast id if needed
    if(!"FIBROBLAST_ID" %in% colnames(res)){
        sampleAnno <- MITOMAP_DATATABLE[!is.na(RNA_ID) & !is.na(FIBROBLAST_ID)]
        sampleAnno <- unique(sampleAnno[,.(RNA_ID, FIBROBLAST_ID)])
        stopifnot(all(sampleAnno[,!duplicated(RNA_ID)]))

        res <- merge(res, sampleAnno)
    }

    # add FraseR columen and say true since all results coming from FraseR
    if(!"FraseR" %in% colnames(res)){
        res <- res[fdr <= fdrCutoff]
        res[, FraseR:=TRUE]
    }

    # add Leafcutter (those are our old results) from Kremer, Bader et al
    if(!"leafcutter" %in% colnames(res)){
        res <- addLeafcutter(res, fdrCutoff, addMissing)
    }

    # add junctionSeq results
    if(!"junctionSeq" %in% colnames(res)){
        warning("Need to add junctionSeq results")
        res[,junctionSeq:=FALSE]
    }

    # add cummingsEtAl paper pipeline
    if(!"cummingsEtAl" %in% colnames(res)){
        warning("Need to add cummings et al results")
        res[,cummingsEtAl:=FALSE]
    }

    if(!"found_MAE" %in% colnames(res)){
        res <- addMAE(res, fdrCutoff, maeCutoff=MAE_LIMIT)
    }

    if(!"found_AbEx" %in% colnames(res)){
        res <- addAbEx(res)
    }

    if(!"found_variant" %in% colnames(res)){
        res <- addVariants(res)
    }

    if(!"isMitoVIP" %in% colnames(res)){
        if(!"MGSA_GO_FULL" %in% ls()){
            eval(load_mito_paper_data(), envir=globalenv())
        }
        source("../gagneurlab_shared/r/disease/get_vip_info_table.R")
        vip_genes <- get_vip_info_table()[,.(
                hgnc_symbol=gene, isMitoVIP=ifelse(causal,"causal", "vip"))]

        res <- merge(res, vip_genes, all.x=TRUE, "hgnc_symbol")
    }

    if(!"numSamplesSameGene" %in% colnames(res)){
        res[p.adj<=fdrCutoff,
            numSamplesSameGene:=length(unique(RNA_ID)),by="hgnc_symbol"
        ]
    }

    # final order
    res2return <- res[,.(seqnames, start, end, width, strand, type,
            RNA_ID, FIBROBLAST_ID, hgnc_symbol, numSamplesSameGene,
            zscore, psiValue, pvalue, p.adj,
            FraseR, leafcutter, junctionSeq, cummingsEtAl,
            found_AbEx, found_MAE, found_variant, isMitoVIP,
            SOLVED, KNOWN_MUTATION
    )][order(zscore)]

    # report removed entries
    message("Removed the following columns:\n\t", paste(collapse=" ",
        colnames(res)[!colnames(res) %in% colnames(res2return)])
    )

    return(res2return)
}

#'
#' add the leafcutter results to the final result table
#'
addLeafcutter <- function(res, fdrCutoff=0.1, addMissing=FALSE){
    if(!"leafcutter_all" %in% ls(envir=.GlobalEnv)){
        if(!"FILE_splice_all_leafcutter_results" %in% ls(envir=.GlobalEnv)){
            eval(load_mito_paper_data(), envir=globalenv())
        }
        eval(envir=globalenv(),
            leafcutter_all <- readRDS(FILE_splice_all_leafcutter_results))
    }
    leafcutter2merge <- leafcutter_all[fdr <= fdrCutoff]
    leafcutter2merge <- leafcutter2merge[!is.na(hgncid),
            .(RNA_ID=RNA_IDs, hgnc_symbol=hgncid)
    ]
    leafcutter2merge <- unique(leafcutter2merge[,
            .(hgnc_symbol=unlist(strsplit(hgnc_symbol, ",")), leafcutter=TRUE),
            by=RNA_ID
    ])
    stopifnot(any(!grepl("[,;:-]", leafcutter2merge[,RNA_ID])))

    # add it to the result table and set the others to FALSE
    res <- merge(res, leafcutter2merge, all.x=TRUE, all.y=addMissing,
            by=c("RNA_ID", "hgnc_symbol")
    )
    res[is.na(leafcutter),leafcutter:=FALSE]

    return(res)
}

#'
#' add aberrent gene expression results
#'
addAbEx <- function(res, fdrCutoff){
    abEx2merge <- data.table(gather("hgnc_symbol", key="RNA_ID",
        data=data.table(
            hgnc_symbol=rownames(signi_counts_sample), signi_counts_sample
        )
    ))
    colnames(abEx2merge) <- c("hgnc_symbol", "FIBROBLAST_ID", "found_AbEx")

    # add it to the result table and set the others to FALSE
    res <- merge(res, abEx2merge, all.x=TRUE,
            by=c("FIBROBLAST_ID", "hgnc_symbol")
    )
    res[is.na(found_AbEx),found_AbEx:=FALSE]

    return(res)
}

#'
#' add MAE to results
#'
addMAE <- function(res, fdrCutoff, maeCutoff=0.8){
    if(!"mae_all" %in% ls()){
        mae_all <- readRDS(MONO_ALLELIC_EXPRESSION_FILE)
    }

    # filter mae
    mae2merge <- mae_all[padj <= fdrCutoff & alt_freq >= maeCutoff][,
            .(RNA_ID=rna_id, hgnc_symbol=hgncid, found_MAE=TRUE)
    ]

    # add it to the result table and set the others to FALSE
    res <- merge(res, mae2merge, all.x=TRUE, by=c("RNA_ID", "hgnc_symbol"))
    res[is.na(found_MAE),found_MAE:=FALSE]

    return(res)
}

#'
#' add variants to results
#'
addVariants <- function(res){
    variants_all <- rbindlist(fill=TRUE, sapply(unique(res[,FIBROBLAST_ID]),
        function(fib){
            ans <- get_patient_report_result(fib)
            if(length(ans) > 1 && !is.null(ans[["VARIANTS"]])){
                ans <- ans[["VARIANTS"]]
                ans[,FIBROBLAST_ID:=fib]
            } else {
                ans <- data.table()
            }
            ans
        }
    ))

    if(nrow(variants_all) == 0){
        res[,found_variant:=character(0)]
        return(res)
    }

    var2merge <- filter_data(filter_exome(
        filter_prot_affect(filter_rare(variants_all))
    ))
    var2merge <- var2merge[,
            .(found_variant=paste(mstype, gt, collapse=";"))
            ,by="FIBROBLAST_ID,hgncid"
    ]
    setnames(var2merge, "hgncid", "hgnc_symbol")

    # add it to the result table and set the others to FALSE
    res <- merge(res, var2merge,
            by=c("FIBROBLAST_ID", "hgnc_symbol"), all.x=TRUE
    )
    # leave na if no match was found

    return(res)
}

#'
#' load the mito paper data functions and variables
#'
load_mito_paper_data <- function(){
    currWD <- getwd()
    try({
        setwd("../mitomultiomics/")
        source("src/r/paper_natgen/load_paper_data.R")
    })
    setwd(currWD)
}

swapPsiSiteValue <- function(res){

    # swap zscore to - x
    res[res$type=="psiSite"]$zscore <- - res[res$type=="psiSite"]$zscore

    # swap psi value to 1 - x
    res[res$type=="psiSite"]$psiValue <- 1 - res[res$type=="psiSite"]$psiValue

    return(res)
}

getDTTable <- function(dt, colOrder=TRUE, doRound=3){
    dt <- copy(dt)
    colHGNC <- ifelse("hgnc_symbol" %in% colnames(dt), "hgnc_symbol", "hgncid")
    colChr  <- ifelse("seqnames" %in% colnames(dt), "seqnames", "chr")
    if(!"start" %in% colnames(dt)){
        dt[,start:=pos]
    }
    if(!"end" %in% colnames(dt)){
        dt[,end:=start]
    }
    if(!"seqnames" %in% colnames(dt)){
        dt[,end:=start]
    }
    # get links
    dt[,genecards:=get_html_link(get(colHGNC), website="genecards", TRUE)]
    dt[,hgnc:=get_html_link(get(colHGNC), website="hgnc", TRUE)]
    dt[,entrez:=get_html_link(get(colHGNC), website="entrez", TRUE)]
    dt[,locus:=get_html_link(paste0(get(colChr), ":", start, "-", end), website="locus", TRUE)]

    if("GMIM" %in% colnames(dt)){
        dt[,omim:=get_html_link(GMIM, website="omim", TRUE)]
    }

    #'
    #' Round numbers
    #'
    if(isTRUE(doRound)){
        doRound <- 3
    }
    if(is.numeric(doRound)){
        for(i in which(sapply(dt, is.numeric))){
            dt[,c(colnames(dt)[i]):=list(signif(get(colnames(dt)[i]), doRound))]
        }
    }

    # set correct order
    if(isTRUE(colOrder)){
        colOrder <- unique(c("sampleID", "genecards", "p.adj", "deltaPsi",
                "type", "numSamplesPerGene", "numEventsPerGene",
                "numSamplesPerJunc", "isMitoVIP", "omim", "PMIM", "PINH",
                "locus", "hgnc", colOrder))
    } else if(!is.character(colOrder)){
        colOrder <- colnames(dt)
    }

    colOrder <- colOrder[colOrder %in% colnames(dt)]
    DT::datatable(dt[,..colOrder], options=list(scrollX=TRUE), escape=FALSE)
}
