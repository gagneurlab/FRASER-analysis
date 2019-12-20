devtools::load_all("../FRASER")

allTissues <- sort(gsub("raw-", "", list.files("./Data/paperPipeline/datasets/savedObjects/", pattern="raw-")))
allTissues <- grep("^(Cells_Leukemia_cell_line_CML|MLL|Kremer|DROP_example|example|SCZ|ProkischFib|SimulationBB|ProkischFull|SimulationDM)$", allTissues, value=TRUE, invert=TRUE, perl=TRUE)
length(allTissues)
allTissues

# outDir <- "./Data/paperPipeline/GTEx-Splice-Map"
outDir <- "/s/public_webshare/project/fraser/GTEx-splicing-map"
fdsDir <- "./Data/paperPipeline/datasets/savedObjects"

tissue <- allTissues[1]
for(tissue in allTissues){
    message(date(), " : work on tissue: ", tissue)
    
    fdsFile <- paste0(fdsDir, "/raw-", tissue, "/fds-object.RDS")
    outJunctionFile   <- paste0(outDir, "/", tissue, "_junction_cts.tsv.gz")
    outSpliceSiteFile <- paste0(outDir, "/", tissue, "_splice_site_cts.tsv.gz")
    
    tissue
    fdsFile
    outJunctionFile
    outSpliceSiteFile
    
    # read file
    fdsRaw <- loadFraseRDataSet(file=fdsFile)
    fdsRaw
    
    # extract counts
    countdt <- as.data.table(K(fdsRaw, type="psi5"))
    annodt <- as.data.table(rowRanges(fdsRaw, type="psi5"))[,.(seqnames, start, end, strand, startID, endID, maxCount)]
    outJunctionTab <- cbind(annodt, countdt)
    outJunctionTab
    
    countdt <- as.data.table(K(fdsRaw, type="psiSite"))
    annodt <- as.data.table(rowRanges(fdsRaw, type="psiSite"))[,.(seqnames, start, end, strand, spliceSiteID, type)]
    outSpliceSiteTab <- cbind(annodt, countdt)
    outSpliceSiteTab
    
    
    # write file
    if(!dir.exists(dirname(outFile))){
        dir.create(dirname(outFile), recursive=TRUE)
    }
    message(date(), ": Writing counts with ", paste(dim(outJunctionTab), collapse=" x "), " to file: ", outJunctionFile)
    fwrite(outJunctionTab, file=outJunctionFile, sep="\t")
    message(date(), ": Writing counts with ", paste(dim(outSpliceSiteTab), collapse=" x "), " to file: ", outSpliceSiteFile)
    fwrite(outSpliceSiteTab, file=outSpliceSiteFile, sep="\t")
}


