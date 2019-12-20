#'---
#' title: Create paper plots
#' author: Ines Scheller
#' wb:
#'  input:
#'   - allOut: '`sm expand("Output/html/FraseR/{dataset}/{{AE_impl}}_autoencoder_fit.html", dataset=config["datasets"])`'
#'  output:
#'   - outPlots: '`sm expand(config["FIGDATADIR"] + "/resultsOverview/{dataset}__{{AE_impl}}_plots.RDS", dataset=config["datasets"])`'
#'   - wBhtml: '`sm config["FIGDATADIR"] + "/resultsOverview/paper_plots_{AE_impl}.html"`'
#'  type: noindex
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=15, fig.height=7)

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(AE_impl="PCA-BB-Decoder")
    parseWBHeader2("./Scripts/DataAnalysis/paperPlots.R",
            wildcards=wildcards, rerun=TRUE)
}

#+ input
datasets          <- snakemake@config$datasets
AE_implementation <- snakemake@wildcards$AE_impl
workingDir  <- file.path(CONFIG$DATADIR, "datasets")

#+ output
outPlots <- snakemake@output$outPlots

#'
#' # Get plots for each dataset
datasets

#+ plot sample rank and qq-plot for all psi types and datasets
plots <- lapply(datasets, function(dataset){

    print(dataset)

    fds <- loadFraseRDataSet(dir=workingDir, name=paste0(dataset, "__", AE_implementation))

    message(date(), ": Creating plots for dataset ", dataset)
    plot_list=list()
    plot_list[["SampleRank"]] <- plotAberrantPerSample(fds)
    plot_list[["GlobalQQ"]]   <- plotGlobalQQPerGene(fds)
    plot_list[["GlobalQQ"]]   <- plotQQ(fds, global=TRUE)
    
    saveRDS(plot_list, file=outPlots[which(grepl(pattern=dataset, x=outPlots))])

    g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
    g
})

