import os
import os.path
import pandas as pd

LeafcutterAnnoDict = {}

def getLeafcutterAnno(dataset):
    annoFile = (config["DATADIR"] + "/processedData/leafcutter/" + 
            dataset + "/rlds_anno.tsv")
    if annoFile in LeafcutterAnnoDict:
        return(LeafcutterAnnoDict[annoFile])
    
    if os.path.isfile(annoFile) is False:
        print("Please rerun the pipeline after generating the anno file" +
                "for dataset: " + dataset)
        anno = pd.DataFrame(columns=['sampleID', 'condition'])
    else:
        #print("Load anno file: " + annoFile)
        anno = pd.read_csv(annoFile, sep="\t")
    
    LeafcutterAnnoDict[annoFile] = anno
    return(anno)
    

def getInputLeafcutterDS(wildcards):
    if hasattr(wildcards, "dataset") is False:
        print("Did not find any wildcards. It should contain dataset!")
        return(list())
        
    dataset = wildcards.dataset
    anno = getLeafcutterAnno(dataset)
    
    ans = [config["DATADIR"] + "/processedData/leafcutter/" + dataset + 
        "/ds_testing/" + x + "_versus_rest/leafcutter_ds_cluster_significance.txt" 
        for x in anno["condition"].unique()]
    
    return(list(ans))
    

def getInputLeafcutterCounts(wildcards):
    dataset = wildcards.dataset
    anno = getLeafcutterAnno(dataset)
    
    ans = [config["DATADIR"] + "/processedData/leafcutter/" + dataset +
        "/leafcutter_processed/junction_files/" + x + ".junc" 
        for x in anno["sampleID"].unique()]
    
    return(list(ans))
    

def getInputLeafcutterBams(wildcards):
    dataset = wildcards.dataset
    anno = getLeafcutterAnno(dataset)
    
    ans  = anno[anno.sampleID == wildcards.sampleID]["bamFile"]
    
    return(list(ans))
    


