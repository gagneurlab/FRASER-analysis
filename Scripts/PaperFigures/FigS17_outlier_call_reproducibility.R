#'---
#' title: Paper figure S17 reproduibilty
#' author: Christian Mertes
#' wb:
#'  threads: 10
#'  input:
#'   - data:   '`sm config["DATADIR"] + "/GTEx_variant_enrichment/rareSplicing__0.0_reproducability.RDS"`'
#'   - table:  '`sm config["DATADIR"] + "/GTEx_variant_enrichment/rareSplicing__0.0_reproducability.tsv.gz"`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/FigureS17_outlier_call_reproducibility.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/FigureS17_outlier_call_reproducibility.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    options <- c("--configfile", "wbuild.yaml")
    parseWBHeader2("Scripts/PaperFigures/FigS17_outlier_call_reproducibility.R",
        rerun=TRUE, options=options)
}

#+ input
data_file  <- snakemake@input$data
table_file <- snakemake@input$table
AE_NAME    <- CONFIG$AE_IMPLEMENTATION
outPng     <- snakemake@output$outPng
outPdf     <- snakemake@output$outPdf


data_file
table_file
AE_NAME
outPng


#' 
#' Read in data
#' 
data <- readRDS(data_file)
res <- readRDS(table_file)

names(data$datatables) <- c("dt2p", "dt2p9p", "dt2p7p", "dt2p5p")
attach(data$datatables)


#' 
#' rename methods
renameTheMethod <- function(dt){
    if(is.factor(dt$Method) && grepl("FRASER", levels(dt$Method))){
        return(dt)
    }
    dt[,Method:=factor(Method, levels=c("PCA_p", "BB_p", "LeafcutterMD_p", "SPOT_p"))]
    dt <- dt[!is.na(Method)]
    levels(dt$Method) <- mName4Plot(levels(dt$Method), removeTest=TRUE, AE_Name=AE_NAME)
    dt
}
dt2p <- renameTheMethod(dt2p)
dt2p5p <- renameTheMethod(dt2p5p)
dt2p7p <- renameTheMethod(dt2p7p)
dt2p9p <- renameTheMethod(dt2p9p)

#'
#'
#'
gt1 <- ggplot(res[,.(Method, total)], aes(total, fill=Method)) + 
    geom_histogram(position="dodge") + 
    scale_fill_brewer(palette="Dark2") +
    ylab("Number of tested events") + 
    xlab(bquote("Number of tissues event is tested")) + 
    theme_bw() + 
    grids() + 
    scale_y_log10()
gt1

g1 <- ggplot(dt2p[hits5 != 0], aes(x=hits5, fill=Method)) + 
    geom_bar(position="dodge") + 
    theme_bw() + 
    facet_wrap(~spliceVariant) + 
    scale_fill_brewer(palette="Dark2") + 
    ylab("Number of events") + 
    xlab(bquote("Number of tissues outlier is present (" ~ italic(P) < 10^-5 ~ ")")) + 
    grids() + 
    scale_y_log10()
g1

g2 <- ggplot(dt2p[hits7 != 0], aes(x=hits7, fill=Method)) + 
    geom_bar(position="dodge") + 
    theme_bw() + 
    facet_wrap(~spliceVariant) + 
    scale_fill_brewer(palette="Dark2") + 
    ylab("Number of events") + 
    xlab("Number of tissues outlier is present")
g2


#' 
#' Plot percentages
#' 
plotPercentage <- function(dt, value){
    ggplot(dt, aes(y=freq*100, x=names, fill=Method)) + 
        geom_bar(stat="identity", position="dodge") + 
        labs(x=bquote("Number of tissues outlier is present (" ~ italic(P) < 10^-.(value) ~ ")"),
             y="Percentage\nwithin Method") + 
        scale_fill_brewer(palette="Dark2") + 
        theme_bw() + 
        grids()
}

gg9p <- plotPercentage(dt2p9p, "9")
gg9p

gg7p <- plotPercentage(dt2p7p, "7")
gg7p

gg5p <- plotPercentage(dt2p5p, "5")
gg5p



#'
#' Arrange the plots
#'
g <- ggarrange(labels=letters[1:4], ncol=1, common.legend=TRUE, legend="bottom",
    g1,
    gg5p,
    gg7p,
    gg9p)
g


#+ save figure
factor <- 0.65
outPng
ggsave(outPng, g, width = 13*factor, height = 17*factor)
ggsave(outPdf, g, width = 13*factor, height = 17*factor)

