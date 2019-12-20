#'
#' MAIN R CONFIGURATION FOR FRASER PIPELINE
#'
library(reticulate)
if(file.exists("/home/mertes/miniconda3/bin/python")){
    use_python("/home/mertes/miniconda3/bin/python", required=TRUE)
    Sys.setenv(PATH=paste0("/home/mertes/miniconda3/bin:", Sys.getenv("PATH")))
} else {
    use_python("/opt/modules/i12g/anaconda/3-5.0.1/envs/omicsOUTRIDER/bin/python", required=TRUE)
    Sys.setenv(PATH=paste0("/opt/modules/i12g/anaconda/3-5.0.1/bin:", Sys.getenv("PATH")))
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

print_comment("Finished with loading the .Rprofile.\n#\n# Now you can start curing diseases! :)")




if(!any(grepl("i12g/gcc/", Sys.getenv("LDPATH")))){
    gcc_version <- "7.4.0"
    Sys.setenv(LDPATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib:",            Sys.getenv("LDPATH")))
    Sys.setenv(LDPATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib64:",          Sys.getenv("LDPATH")))
    Sys.setenv(LD_LIBRARY_PATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib:",   Sys.getenv("LD_LIBRARY_PATH")))
    Sys.setenv(LD_LIBRARY_PATH=paste0("/opt/modules/i12g/gcc/", gcc_version, "/lib64:", Sys.getenv("LD_LIBRARY_PATH")))
}
