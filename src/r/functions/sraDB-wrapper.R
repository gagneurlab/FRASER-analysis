#'
#' load SRAdb and needed packages
#'
suppressPackageStartupMessages({
    library(SRAdb)
    library(stringr)
    library(data.table)
    library(parallel)
    library(BBmisc)
})

#'
#' get SraDB file (download it if needed)
#'
getSraDBConnection <- function(sra_con="Data/filemapping/SRAmetadb.sqlite"){
    # run it with default if its missing or NULL
    if(missing(sra_con) || is.null(sra_con)){
        return(getSraDBConnection(formals(getSraDBConnection)$sra_con))
    }

    # check what kind of input we got
    if(isScalarCharacter(sra_con)){
        # check if file exists or if we need to download it
        if(!file.exists(sra_con)){
            message("Need to download and unzip the SraDB file.")
            message("This can take up to 20 min or longer.")
            sra_con <- getSRAdbFile(destdir = dirname(sra_con))
            message("Finished with getting the SraDB file.")
        }

        # get connection
        real_sra_con <- dbConnect(SQLite(),sra_con)
    } else if(class(sra_con) == "SQLiteConnection"){
        real_sra_con <- sra_con
    } else {
        stop("Can not detect a SraDB object within the given object: ",
                str(sra_con)
        )
    }

    return(real_sra_con)
}

#'
#' converting a string which contains "key: value" pairs separated by " || "
#' into a data.table.
#' The keys are the colnames and the values are the row entries.
#'
getDataTableFromString <- function(str){
    strsp <- unlist(strsplit(str, " *\\|\\| *"))
    strkv <- str_split(strsp, ": ", 2)
    value <- lapply(strkv, "[", 2)
    names(value) <- gsub(" ", "_", sapply(strkv, "[", 1))
    as.data.table(value)
}

#'
#' extract for the given SRR ID all relevant information
#' and return it as a data.table
#'
getSampleDataTable <- function(id, sra_con=NULL, sraData=NULL){
    if(is.null(sraData)){
        sra_con <- getSraDBConnection(sra_con)
        neededCols <- "sample_attribute,experiment_attribute,run_attribute"
        ids <- as.data.table(sraConvert(id, sra_con=sra_con))
        stopifnot(dim(ids)[1] == 1)
        rssra <- dbGetQuery(sra_con, paste0("select ", neededCols,
                " from sra where run_accession == '", ids[,"run"], "'"
        ))
        dbDisconnect(sra_con)
    } else {
        stopifnot(class(sraData) == "data.frame")
        rssra <- sraData
        ids <- rssra[,c("study", "submission", "sample", "experiment", "run")]
    }

    sampledt <- getDataTableFromString(rssra$sample_attribute)
    expdt    <- getDataTableFromString(rssra$experiment_attribute)
    rundt    <- getDataTableFromString(rssra$run_attribute)

    finaldt  <- do.call(cbind, list(ids, sampledt, expdt, rundt))

    # remove duplicated columns
    if(any(duplicated(colnames(finaldt)))){
        finaldt  <- finaldt[,!duplicated(colnames(finaldt))]
    }

    return(finaldt)
}

#'
#' returns all samples and its information in a data.table for a given
#' study_accession number from the SRA database
#'
getSRAProjectTable <- function(study_accession="SRP012682", sra_con=NULL,
            mc.cores=10){
    sra_con <- getSraDBConnection(sra_con)

    # get data for full study
    neededCols <- paste(sep=",",
            "run_accession", "sample_attribute",
            "experiment_attribute", "run_attribute"
    )
    fullStudydf <- dbGetQuery(sra_con, paste0("select ", neededCols,
            " from sra where study_accession == '", study_accession, "'"
    ))

    # get all ids for study
    ids <- sraConvert(study_accession, sra_con=sra_con)
    dbDisconnect(sra_con)

    datatablels <- mclapply(1:dim(fullStudydf)[1], mc.cores=mc.cores,
        FUN=function(idx){
            data <- cbind(fullStudydf[idx,], ids[idx,])
            getSampleDataTable(id=fullStudydf[idx,1], sraData=data)
        }
    )

    # unlist result and return it
    datatable <- rbindlist(datatablels, fill=TRUE)
    return(datatable)
}
