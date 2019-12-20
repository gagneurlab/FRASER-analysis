#'---
#' title: Paper figure (TAZ)
#' author: Ines Scheller
#' wb:
#'  input:
#'   - sashimi: 'Data/figures/TAZ.png'
#'   - fdsin:   '`sm config["DATADIR"] + "/datasets/savedObjects/Kremer__" + config["AE_IMPLEMENTATION"] + "/pajdBetaBinomial_psiSite.h5"`'
#'  output:
#'   - outPng: '`sm config["FIGDIR"] + "/Figure7_TAZ.png"`'
#'   - outPdf: '`sm config["FIGDIR"] + "/Figure7_TAZ.pdf"`'
#' output:
#'  html_document
#'---

#+ echo=FALSE
source("./src/r/config.R")

# Debug data
if(FALSE){
    source(".wBuild/wBuildParser2.R")
    parseWBHeader2("Scripts/PaperFigures/Fig7_TAZ.R", rerun=TRUE)
}

#+ input
sashimiPlot <- readPNG(snakemake@input$sashimi)
fdsin       <- snakemake@input$fdsin
FraseR_implementation <- snakemake@config$AE_IMPLEMENTATION

FraseR_implementation
gene <- "TAZ"
sample <- "MUC1398"
fibID  <- "74116"

#+ output
outPdf <- snakemake@output$outPdf
outPng <- snakemake@output$outPng

fdsin
outPng

#+ load fds
fds <- loadFraseRDataSet(file=fdsin)
fds <- annotateRanges(fds)

#+ get pvalues per gene
sampleIndex  <- which(samples(fds) == sample)
grOfInterest <- GRanges(seqnames="chrX", strand="+", type="psi3",
        ranges=IRanges(start=153641905, end=153642437))
type         <- grOfInterest$type

dt <- getPlottingDT(fds, axis="col", idx=sampleIndex, aggregate=TRUE, type=type)
dt[featureID == "TAZ"]

#'
#' Assemble TAZ figure
#'
#+ assemble TAZ figure
g1tmp <- plotVolcano(fds, sample, type=type, aggregate=TRUE)
g1 <- g1tmp + annotate("text", label="TAZ", 
            y=g1tmp$data[featureID == "TAZ"][1, -log10(pval)], 
            x=g1tmp$data[featureID == "TAZ"][1, deltaPsi + 0.14]) + 
    xlab(expression(paste(Delta, psi[3]))) + 
    theme(plot.title=element_blank())
g1


g2 <- plotExpression(fds, result=grOfInterest) + 
    theme(plot.title=element_blank()) +
    ylab("Split reads (K)\nfrom exon 4 to exon 5") +
    xlab("All split reads (N)\nat exon 5 acceptor") + 
    scale_x_log10(breaks=xticks + 1, labels=xticks) + 
    scale_y_log10(breaks=yticks + 1, labels=yticks)
g2


g3 <- plotQQ(fds, result=grOfInterest) + 
    theme(plot.title=element_blank())
g3


g4 <- plotExpectedVsObservedPsi(fds, result=grOfInterest) + 
    scale_color_manual(values=c("firebrick", "gray70")) + 
    theme(plot.title=element_blank())
g4


e_factor <- 3.5
g5 <- ggplot() +
    background_image(sashimiPlot) +
    theme_nothing() +
    theme(plot.margin=margin(
            t=.2*e_factor, l=.2, r=-.5, b=.2*e_factor, unit="cm"))
g5


#'
#' Assemble figures
#'
#+ assemble figure for psi5
g <- ggarrange(ncol=2, widths=c(1.2,1),
    ggarrange(nrow=2, ncol=2, labels=LETTERS[1:4], align="hv",
        g1, g2, g4, g3),
    ggarrange(g5, labels=LETTERS[5]))
g


#+ save figure
factor <- 0.7
ggsave(outPng, g, width = 16*factor, height = 7.5*factor)
ggsave(outPdf, g, width = 16*factor, height = 7*factor)

