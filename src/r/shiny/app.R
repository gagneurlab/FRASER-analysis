#'
#' Run the FraseR shiny server
#'
#'
#' author: mertes
#'
#' @examples
#'   options(shiny.reactlog=TRUE)
#'   shiny::runApp("src/r/shiny/", port=3843, host="0.0.0.0")
#'

#' this code is only for our inhouse adaptation
setwd(gsub("/src/r/shiny/?$", "", getwd()))
if(file.exists("src/r/config.R")){
    #'
    #' source our configs
    source("src/r/config.R")

    #'
    #' check that we have access to the hard drives
    message("Check kerberos ticket")
    kerberos.renew()
    if(!kerberos.check(4)){
        kerberos.init(
            minLifeTime=10,
            keyFile="/etc/shiny-server/keytabs/krb5.nfs.ouga03.keytab",
            user="nfs/ouga03.i12g.informatik.tu-muenchen.de")
    }
    stopifnot(kerberos.check(4))

    load_mitomap()
    load_mito_paper_data()

    library(devtools)
    load_all(PKG_ROOT)
}

FRASER_DATA_DIR <- "/s/project/fraser/analysis/prokisch-full"
FRASER_DATA_NAME <- "Full Prokisch Analysis n 2 Batch"
FRASER_DATA_NAME <- "Full Prokisch Analysis"


#'
#' Source needed packages
message("Source packages")
suppressPackageStartupMessages({
    require(FraseR)
    require(plotly)
    require(DT)
    require(shiny)
})

#'
#' load data
message("Load FraseR results")
if(!"fdsShiny" %in% ls()){
    fdsShiny    <<- loadFraseRDataSet(FRASER_DATA_DIR, FRASER_DATA_NAME)
    fdsShinyRes <<- results(fdsShiny)

    if("MITOMAP_DATATABLE" %in% ls()){
        fdsShinyRes <- swapPsiSiteValue(fdsShinyRes)
        resDT <- addSampleNGeneInfosToResults(fdsShinyRes)
    }
}

#'
#' get shiny object
message("Load shiny app")
shinyObj <- FraseRShinyApp(fdsShiny, fdsShinyRes, server=TRUE)


#'
#' the shiny ui object
message("Start shiny app")
shinyApp(ui=shinyObj$shinyObj$uiMain, server=shinyObj$shinyObj$serverMain)

