#'---
#' title: Run all Benchmarks
#' author: Christian Mertes
#' wb:
#'  input:
#'   - benchmarks: '`sm expand("Output/html/Benchmark/inSilicoBenchmark/efreq{efreq}_n{nsamples}_{dataset}_final_benchmark.html", efreq="0.0001", nsamples="200", dataset="gtex-skin-all")`'
#'---

paste(snakemake@input[[1]])
