#'
#' MAIN R CONFIGURATION FOR FRASER PIPELINE
#'
load_wbuild <- function(){
    library(reticulate)
    if(file.exists("/opt/modules/i12g/anaconda/envs/fraser-analysis/bin/python")){
        use_condaenv(
                condaenv = "/opt/modules/i12g/anaconda/envs/fraser-analysis/", 
                conda="/opt/modules/i12g/anaconda/3-2019.10/condabin/conda")
    } else if(file.exists("/home/mertes/miniconda3/bin/python")){
        use_python("/home/mertes/miniconda3/bin/python", required=TRUE)
        Sys.setenv(PATH=paste0("/home/mertes/miniconda3/bin:", Sys.getenv("PATH")))
    }
    py_config()
    if(system('which snakemake') != 0){
        warning("Could not load snakemake/conda environment.")
    }
}

print_comment <- function(...){
    border <- paste0(rep("#", 80), collapse="")
    message("\n", border, "\n#\n# ", ..., "\n#\n", border, "\n")
}
print_comment("Start loading the 'src/r/config.R' ...")

loadConfig <- function(confFile="./wbuild.yaml"){
    if(file.exists(confFile)){
        CONFIG <<- read_yaml(confFile)
    }
}

loadFraseR <- function(){
    load_all(CONFIG$PACKAGE_DIR)
}

##--------------------------------------------
## Hack to bypass the parallel function
## rewriting in the full analysis code
##--------------------------------------------
`parallel<-` <- function(fds, value){
    stopifnot(is(value, "BiocParallelParam"))
    message("Please rewrite the code to use register(...)")
    register(value)
    fds
}

parallel <- function(fds){
    message("Please rewrite the code to use register(...)")
    bpparam()
}

##--------------------------------------------
## required packages
message("Load needed packages")
suppressPackageStartupMessages({
    library(markdown)
    library(knitr)
    library(cowplot)
    library(BBmisc)
    library(DelayedMatrixStats)
    library(devtools)
    library(GenomicAlignments)
    library(ggplot2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(gtable)
    library(plotly)
    library(png)
    library(RColorBrewer)
    library(tidyr)
    library(tools)
    library(yaml)
})

suppressPackageStartupMessages({
    loadConfig()
    loadFraseR()
})

##--------------------------------------------
## global defined functions
sourceFolder <- function(dir){
    DEV_NULL <- sapply(
        list.files(dir, pattern="\\.R$", full.names=TRUE),
        function(f) tryCatch({source(f)}, error=function(e) message(e, "in: ", f))
    )
}

# source extra functions
sourceFolder("src/r/functions")

##--------------------------------------------

if("snakemake" %in% ls()){
    if(!dir.exists("tmp")){
        dir.create("tmp")
    }
    saveRDS(snakemake, "tmp/snakemake.RDS")
}
getSnakemakeTmp <- function(){
    if(file.exists("./tmp/snakemake.RDS")){
        message("Load snakemake to snakemake variable")
        snakemake <<- readRDS("./tmp/snakemake.RDS")
    } else {
        message("snakemake file not found!")
    }
    invisible(TRUE)
}

if(dir.exists("Data/paperPipeline/datasets/savedObjects")){
    if(!rhdf5::h5testFileLocking("Data/paperPipeline/datasets/savedObjects/test.locking.h5")){
        warning("Disabling HDF5 file locking, because it is not supported on this system!")
        rhdf5::h5disableFileLocking()
    }
}


print_comment("Finished with loading the .Rprofile.\n#\n# Now you can start curing diseases! :)")




if(!any(grepl("i12g/gcc/", Sys.getenv("LDPATH")))){
    gcc_version <- "7.4.0"
    Sys.setenv(LDPATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib:",            Sys.getenv("LDPATH")))
    Sys.setenv(LDPATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib64:",          Sys.getenv("LDPATH")))
    Sys.setenv(LD_LIBRARY_PATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib:",   Sys.getenv("LD_LIBRARY_PATH")))
    Sys.setenv(LD_LIBRARY_PATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib64:", Sys.getenv("LD_LIBRARY_PATH")))
}
