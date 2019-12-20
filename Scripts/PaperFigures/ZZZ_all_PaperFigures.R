#'---
#' title: Merge figures
#' author: Ines Scheller, Christian Mertes
#' wb:
#'  input:
#'   - PaperNums: '`sm htmlOutputPath + "/Scripts_PaperFigures_PaperNumbers.html"`'
#'   - Figure1:   '`sm config["FIGDIR"] + "/Figure1_scheme.pdf"`'
#'   - Figure2:   '`sm config["FIGDIR"] + "/Figure2_heatmap_psi5.pdf"`'
#'   - Figure3:   '`sm config["FIGDIR"] + "/Figure3_resultOverview.pdf"`'
#'   - Figure4:   '`sm config["FIGDIR"] + "/Figure4_precRec_clean/Skin_Not_Sun_Exposed_Suprapubic/uniformDistr/precRec_byJunctionGroup_psi5_clean.pdf"`'
#'   - Figure5:   '`sm config["FIGDIR"] + "/Figure5_GTEx_enrichment.pdf"`'
#'   - Figure6:   '`sm config["FIGDIR"] + "/Figure6_Kremer_results.pdf"`'
#'   - Figure7:   '`sm config["FIGDIR"] + "/Figure7_TAZ.pdf"`'
#'   - FigureS01: '`sm config["FIGDIR"] + "/FigureS1_filtering.pdf"`'
#'   - FigureS02: '`sm config["FIGDIR"] + "/FigureS2_heatmap_psi3.pdf"`'
#'   - FigureS03: '`sm config["FIGDIR"] + "/FigureS2_heatmap_psiSite.pdf"`'
#'   - FigureS04: '`sm config["FIGDIR"] + "/FigureS4_finding_q.pdf"`'
#'   - FigureS05: '`sm config["FIGDIR"] + "/FigureS5_qq_plots.pdf"`'
#'   - FigureS06: '`sm config["FIGDIR"] + "/FigureS6_outlier_numbers.pdf"`'
#'   - FigureS07: '`sm config["FIGDIR"] + "/Figure4_precRec/Skin_Not_Sun_Exposed_Suprapubic/uniformDistr/precRec_byJunctionGroup_psi3.pdf"`'
#'   - FigureS08: '`sm config["FIGDIR"] + "/Figure4_precRec/Skin_Not_Sun_Exposed_Suprapubic/uniformDistr/precRec_byJunctionGroup_psiSite.pdf"`'
#'   - FigureS09: '`sm config["FIGDIR"] + "/FigureS9_injectedVsFittedDpsi_PCA-BB-Decoder.pdf"`'
#'   - FigureS10: '`sm config["FIGDIR"] + "/FigureS10_GTEx_enrichment_pval.pdf"`'
#'   - FigureS11: '`sm config["FIGDIR"] + "/FigureS11_IntronRetention.pdf"`'
#'  output:
#'  - figure:  '`sm config["FIGDIR"] + "/Figure_all_main.pdf"`'
#'  - figureS: '`sm config["FIGDIR"] + "/Figure_all_sup.pdf"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---
source("./src/r/config.R")


# source config and load package
mainFigures <- unlist(unique(snakemake@input[grepl("Figure[0-9]",  snakemake@input)]))
supFigures  <- unlist(unique(snakemake@input[grepl("FigureS[0-9]", snakemake@input)]))

message(mainFigures)


message(date(), ': Start with main figures ... ')
system(paste0(
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dNumRenderingThreads=10 -dNOGC ',
    ' -dBandBufferSpace=5000000000 -dBufferSpace=10000000000 -sBandListStorage=memory ',
    '-sOutputFile=', snakemake@output$figure, ' ',
    paste(mainFigures, collapse = ' ')))
system(paste0(
    'gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE ',
    '-dQUIET -dBATCH -sOutputFile=',
    gsub(".pdf$", "_reduced.pdf", snakemake@output$figure), ' ', snakemake@output$figure))



message(supFigures)


message(date(), ': Start with supplement figures ... ')
system(paste0(
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dNumRenderingThreads=10 -dNOGC ',
    ' -dBandBufferSpace=5000000000 -dBufferSpace=10000000000 -sBandListStorage=memory ',
    '-sOutputFile=', snakemake@output$figureS, ' ',
    paste(supFigures, collapse = ' ')))
system(paste0(
    'gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE ',
    '-dQUIET -dBATCH -sOutputFile=',
    gsub(".pdf$", "_reduced.pdf", snakemake@output$figureS), ' ', snakemake@output$figureS))

