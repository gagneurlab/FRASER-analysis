#'---
#' title: Convert figure1 to pdf and png
#' author: Ines Scheller
#' wb:
#'  input:
#'   - Figure1: 'Data/paper/figures_svg/Fig1-FraseR-overview.svg'
#'  output:
#'  - png: '`sm config["FIGDIR"] + "/Figure1_scheme.png"`'
#'  - pdf: '`sm config["FIGDIR"] + "/Figure1_scheme.pdf"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

source('./src/r/config.R')

system(paste0("inkscape -z -b 'white' -d 600 -e ",
              snakemake@output$png, ' ',
              snakemake@input$Figure1))

system(paste0("inkscape -z -d 600 -A ",
              snakemake@output$pdf, ' ',
              snakemake@input$Figure1))
