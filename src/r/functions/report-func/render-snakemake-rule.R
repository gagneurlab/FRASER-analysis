## # /usr/bin/env Rscript --vanilla
#'
#' author: render snakemake rule
#'
if(!"snakemake" %in% ls()){
    stop("this scipt renders a snakemake call based on the snakemake object.")
}
stopifnot("params" %in% slotNames(snakemake))
stopifnot("script" %in% names(snakemake@params))
stopifnot(file.exists(snakemake@params[['script']]))

#library(devtools)
#load_all("../fraser")

library(knitr)
library(rmarkdown)

renderScript  <- snakemake@params[['script']]
wildcardNames <- grep("^$", names(snakemake@wildcards), invert=TRUE, value=TRUE)
paramNames    <- grep("^(script)?$", names(snakemake@params),
        perl=TRUE, value=TRUE, invert=TRUE)
params        <- snakemake@params[paramNames]
params        <- append(params, c(threads=snakemake@threads))
wd            <- getwd()

renderObj <- list(renderScript=renderScript, wildcardNames=wildcardNames,
        params=params, wd=wd, snakemake=snakemake)

#'
#' save it for debuging
#'
rdsFile <- file.path("./tmp/snakemake/", renderScript, paste0(paste(
        wildcardNames, unlist(attributes(snakemake)$wildcards[wildcardNames])
        , collapse="_", sep=":"), ".RDS"))
dir.create(dirname(rdsFile), recursive=TRUE)
saveRDS(renderObj, rdsFile)

opts_knit$set(root.dir = wd)
render(renderScript, output_file = unlist(snakemake@output), params = params)


