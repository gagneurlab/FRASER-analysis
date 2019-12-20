#'
#'
#'
get_idx_for_region <- function(dataset, readType, chr, start, end, filter = 5){
    pos <- which(
        start(slot(dataset, readType))    >= start   &
        end(slot(dataset, readType))      <= end     &
        seqnames(slot(dataset, readType)) == chr
    )

    filt_pos <- which(
        apply(
            assays(slot(dataset, readType)[pos])$rawCounts,
            1,
            max
        ) > filter
    )

    return(pos[filt_pos])
}

get_table_for_idx <- function(dataset, readType, assays=NULL, mcols=NULL, idx){

    if(!is.null(assays)){
        x <- assays[[1]]
        a <- sapply(assays, function(x){
            assays(slot(dataset, readType))[[x]][idx,]
        })
    }
    if(!is.null(mcols)){
        x <- mcols[[1]]
        m <- sapply(mcols, function(x){
            mcols(slot(dataset, readType))[[x]][idx]
        })
    }
    list(a=a,m=m)
}


#'
#' function to update the FraseR documentation and reload the package
#'
install_fraser <- function(){
    unloadNamespace(basename(PKG_ROOT))
    document(PKG_ROOT)
    install(PKG_ROOT)
    library(FraseR)
}


#'
#' run BiocCheck on FraseR
#'
test_fraser <- function(){
    document(PKG_ROOT)
    BiocCheck(PKG_ROOT)
}


#'
#' build the full vignette for FraseR
#'
buildVignetteFraseR <- function(repackFraseR=FALSE){
    # install FraseR
    if(repackFraseR){
        install_fraser()
    }

    # set correct texlive path
    texdir="/opt/modules/i12g/texlive/2016/"
    Sys.setenv(PATH=paste(sep=":",
            file.path(texdir, "texlive/bin/x86_64-linux"),
            file.path(texdir, "bin/x86_64-linux"),
            Sys.getenv("PATH")
    ))

    # remove old function and variables
    rm(list=grep(
            "buildVignetteFraseR|PKG_ROOT|install_fraser",
            ls(), value=TRUE, invert=TRUE
    ))

    # go to build folder
    buildDir <- "FraseR/vignettes/build"
    srcDir <- getwd()
    if(!dir.exists(buildDir)){
        dir.create(buildDir, recursive=TRUE)
    }
    setwd(buildDir)

    # run knitr and pdflatex
    tryCatch({
        # Sweave2knitr(file.path(srcDir, "FraseR/vignettes/FraseR.Rnw"))
        knit(file.path(srcDir, "FraseR/vignettes/FraseR.Rnw"))
        system("pdflatex FraseR && pdflatex FraseR")
    })

    # go back to src folder
    setwd(srcDir)
}


#'
#' load the mitoMultiOmics project
#'
load_mitomap <- function(){
    if(!dir.exists("../mitomultiomics")){
        message("Can't find mitomultiomics project. If you don't need it, ignore this message!")
        return(FALSE)
    }

    # mitomultiomics functions
    wd <- getwd()
    setwd("../mitomultiomics/")
    tryCatch({ source("./src/r/config.R") },
        #warning = function(w) { print(w, "\n"); source("./src/r/mitomap/load_mitomap.R") },
        error = function(e) { print(e, "\n"); source("./src/r/mitomap/load_mitomap.R") }
    )

    # new sample annotation
    setwd(wd)
    setwd("../sample_annotation/")
    tryCatch({ source("./src/r/config.R") },
            error = function(e) { print(e, "\n"); }
    )
    setwd(wd)
}

#'
#'
#'
getMitoMapTestSet <- function(nSamples=30){
    rnaIDsOfInterest <<- c(
        "MUC1344", "MUC1365", "MUC1350", "MUC1361", "76626",
        "MUC1369", "76631", "MUC1406", "MUC1398", "MUC1391",
        "MUC1367", "MUC1385", "MUC1360", "MUC1373"
    )
    grOfInterest <<- GRanges(
        seqnames = c("chr3", "chr19", "chr22", "chr5", "chr6", "chr19", "chr2"),
        ranges = IRanges(
            start = c(118996051, 7126380, 29883155, 71032670, 158109504, 6361452, 32946972),
            end   = c(119812193, 8056527, 30030868, 71067689, 158168280, 6368908, 33399509)
        ),
        mcols = DataFrame(
            gene_symbol = c("TIMMDC1", "MCOLN1", "MTMR3", "GFT2H2", "SERAC1", "CLPP", "LTBP1")
        )
    )

    # get good rna ids
    rna_ids <- MITOMAP_DATATABLE[
        !is.na(RNA_ID) &
            TISSUE == "FIBROBLAST" &
            is.na(TRANSDUCED_GENE) &
            GROWTH_MEDIUM == "GLU" &
            (is.na(KNOWN_MUTATION) | KNOWN_MUTATION != "NHDF"),
        ][,unique(RNA_ID)]

    # cleanup rna ids
    rnaIDsOfInterest <- rnaIDsOfInterest[rnaIDsOfInterest %in%rna_ids]
    rna_ids <- rna_ids[!rna_ids %in% rnaIDsOfInterest]

    # calculate missing number of extra samples
    nExtraSamplesNeeded <- 1:min(
            max(0, nSamples - length(rnaIDsOfInterest)),
            length(rna_ids)
    )
    if(any(nExtraSamplesNeeded<1)){
        nExtraSamplesNeeded <- integer(0)
    }

    # fill up to samples
    rna_ids_to_take <- c(rnaIDsOfInterest, rna_ids[nExtraSamplesNeeded])

    # get bam files
    bamfiles <- unlist(BamFileList(
        sapply(rna_ids_to_take, get_helmholtz_file, source = "RNA", type = "BAM.STAR")
    ))

    # create setting options
    settings <- FraseRDataSet(colData = data.table(
            sampleID=rna_ids_to_take,
            bamFile=bamfiles,
            group=1:length(rna_ids_to_take
    )))

    # set to interesting regions only
    bamWhich(settings) <- grOfInterest

    # update multicore params
    parallel(settings) <- MulticoreParam(30, progressbar = T)

    return(settings)
}

getFraseRRunTime <- function(startTime, progress){
    if(is.character(startTime)){
        startTime <- as.POSIXct(strptime(startTime, format = "%c"))
    }

    # get diff
    currDiff <- difftime(as.POSIXct(Sys.time()), startTime, units="hours")

    # print current runtime
    message(date(), ":\n\tCurr Runtime in hours: ", round(currDiff,2))

    # print full runtime
    message("\tFull Runtime in hours: ", round(100*currDiff/progress, 2))
}

#'
#' Returns a colData object for the FraseR package based on the given samples
#' Containing at least following columns
#'   * sampleID
#'   * bamFile
#'   * condition
#'
getColDataForFraseR <- function(sampleIDs){

    # get only unique samples
    sampleIDs <- unique(sampleIDs)

    # default col data
    colData <- data.table(sampleID=character(0), bamFile=character(0))

    # check mitomap for sampleIDs
    if(!"MITOMAP_DATATABLE" %in% ls(envir=.GlobalEnv)){
        eval(load_mitomap(), envir = globalenv())
    }
    if(any(MITOMAP_DATATABLE[,RNA_ID %in% sampleIDs])){
        colData2add <- copy(MITOMAP_DATATABLE[RNA_ID %in% sampleIDs])
        colData2add[,sampleID:=RNA_ID]
        colData2add[,bamFile:=sapply(sampleID,
                get_helmholtz_file, source="RNA", type="BAM.STAR"
        )]

        message("REMOVE MUC2214 because of the BAM file problem ... take care later about it!!!!")
        colData2add <- colData2add[RNA_ID != "MUC2214"]

        colData <- merge(colData, colData2add, all=TRUE)
    }

    missingSamples <- sampleIDs[!colData[,sampleIDs %in% sampleID]]
    # get Aligne out new samples
    if(length(missingSamples) > 0){
        colData2add <- data.table(sampleID=missingSamples, bamFile=file.path(
                RAWDATADIR, missingSamples, "RNAout/paired-endout/",
                paste0(missingSamples, "Aligned.out.sort.bam")))
        colData2add <- colData2add[file.exists(bamFile)]
        colData <- merge(colData, colData2add, fill=TRUE, all=TRUE)
    }

    missingSamples <- sampleIDs[!colData[,sampleIDs %in% sampleID]]
    if(length(missingSamples) > 0){
        stop("Please edit this if a sample is still missing!")
    }

    # sanity check
    stopifnot(nrow(colData) == length(sampleIDs))
    stopifnot(colData[,sampleIDs %in% sampleID])
    stopifnot(colData[,!is.na(bamFile)])
    stopifnot(colData[,file.exists(bamFile)])

    if(!"condition" %in% colnames(colData)){
        colData[,condition:=1:.N]
    }

    colData <- DataFrame(colData)
    rownames(colData) <- colData$sampleID

    return(colData)
}


#'
#' render a rmakrdown file
#'
CURR_REPORT_FILE=NULL
renderReport <- function(file=CURR_REPORT_FILE){
    if(is.null(file)){
        message("Please set the CURR_REPORT_FILE parameter",
            " or provide it as a parameter.")
        return()
    }
    library(knitr)
    library(rmarkdown)
    options(width=200)
    opts_knit$set(root.dir = getwd())
    render(file)
}

#'
#' Nice data table with download button
#'
downloadableDT <- function(data, tablename="dt-table", scrollX=TRUE){
    if(!is.data.table(data)){
        data <- as.data.table(data)
    }

    DT::datatable(data, extensions='Buttons', rownames=FALSE,
            options=list(scrollX=scrollX, dom='lBfrtip', server=FALSE,
                    buttons=list('colvis', 'copy', list(
                            extend='collection',
                            buttons = list(
                                    list(extend='csv',   filename = tablename),
                                    list(extend='excel', filename = tablename),
                                    list(extend='pdf',   filename = tablename)),
                            text='Download'))))
}


#'
#' write tsv file
#'
write_tsv <- function(x, file, row.names = FALSE, ...){
    write.table(x=x, file=file, quote=FALSE, sep='\t', row.names= row.names, ...)
}


#'
#' GTEx tissue names for plotting
#'
dName4plot <- function(x, removeMethod=FALSE){
    if(isTRUE(removeMethod)){
        x <- gsub("__.*", "", x)
    }
    ans <- toTitleCase(gsub("^raw-", "", gsub("_", " ", x)))
    ans <- gsub(" BA24| BA9| c 1|", "", ans, perl=TRUE)
    
    ans 
}

#'
#' Method names to nice plotting names
#' 
mName4Plot <- function(x, removeTest=TRUE, AE_Name){
    if(isTRUE(removeTest)){
        x <- gsub("_n?[pz]$", "", x)
    }
    sapply(x, function(i){
        suffix <- ""
        if(endsWith(i, "_p")){
            suffix <- "_p"
        } else if(endsWith(i, "_z")){
            suffix <- "_z"
        } else if(endsWith(i, "_np")){
            suffix <- "_np"
        }
        ans <- switch(gsub("_[pz]$", "", i),
            `AE_Name` = "FRASER",
            BB = "BetaBinom",
            Leafcutter = "Kremer et al.",
            PCA = "PCA",
            `PCA-BB-Decoder` = "BB-AE-weight",
            `PCA-BB-Decoder-no-weights` = "BB-AE",
            stop("Could not find method:", i)
        )
        paste0(ans, suffix)
    })
}

#'
#' Correct Extraction protocols
#'
SMNABTCHT4plot <- function(x){
    correctName <- function(y){
        sapply(y, function(yy) switch (yy,
            "RNA Extraction from Paxgene-derived Lysate Plate Based" = "PAX-lysate",
            "RNA isolation_PAXgene Blood RNA (Manual)"               = "PAX-Blood",
            "RNA isolation_PAXgene Tissue miRNA"                     = "PAX-miRNA",
            "RNA isolation_QIAGEN miRNeasy"                          = "QIAGEN",
            "RNA isolation_Trizol Manual (Cell Pellet)"              = "Trizol",
            yy
        ))
    }
    if(any(class(x) %in% c("data.frame", "data.table", "DataFrame"))){
        x$SMNABTCHT <- correctName(x$SMNABTCHT)
    } else {
        x <- correctName(x)
    }
    x
}


#'
#' Correct gender in coldata
#'
GENDER4plot <- function(x){
    correctName <- function(y){
        sapply(y, function(yy) switch(yy,
                "1" = "male",
                "2" = "female",
                yy
        ))
    }
    if(any(class(x) %in% c("data.frame", "data.table", "DataFrame"))){
        x$GENDER <- correctName(x$GENDER)
    } else {
        x <- correctName(x)
    }
    x
}


#'
#' Correct the RIN number for plotting
#'
SMRIN4plot <- function(x){
    if(is.numeric(x$SMRIN)){
        cuts <- cut(x$SMRIN, 5:10)
        levels(cuts) <- gsub(",", "-", gsub("\\(|\\]", "", levels(cuts)))
        x$SMRIN <- cuts
    }
    x
}


#'
#' Get max cpus per run
#'
bpMaxWorkers <- function(max=bpworkers(), ...){
    ans <- min(unlist(list(max, ...)))

    # local limits (cpus/slurm)
    ans <- min(ans, bpworkers())
    ans <- min(ans, as.numeric(Sys.getenv("SLURM_NPROCS")), na.rm=TRUE)

    ans
}

#'
#' gtable_remove_grobs bugfix!!!
#'
#' TODO submit bugfix report to gtable package
#'
gtable_remove_grobs <- function (table, names, ...){
    kept_names <- table$layout$name[!(table$layout$name %in% names)]
    gtable::gtable_filter(table, paste0("^(",
            paste(kept_names, sep = "", collapse = "|"), ")$"), ...)
}

