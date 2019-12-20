import pandas as pd
import pathlib
import wbuild

configfile: "wbuild.yaml"
htmlOutputPath = config['htmlOutputPath']

# sub-part of enrichment pipeline
include: "Scripts/Leafcutter/leafcutter_snakemake_functions.py"
include: "Scripts/GTExEnrichment/enrichment.snake"

# parse main wbuild script
include: ".wBuild/wBuild.snakefile"

rule all:
	input: rules.Index.output, htmlOutputPath + "/readme.html"
	output: touch("Output/all.done")


rule prepub:
    shell: """
        rsync -Ort {config[htmlOutputPath]} {config[webDir]}
        chmod -R uo+rX {config[webDir]}
    """

rule defineDatasets:
    input: 
        expand(config["DATADIR"] + "/annotations/{dataset}.tsv", dataset=config["datasets"] + config["EnrichmentTissues"]),
        expand(config["DATADIR"] + "/processedData/leafcutter/{dataset}/rlds_anno.tsv", dataset=config["datasets"] + config["EnrichmentTissues"])

rule supplement_pdf:
    input:
        tex='src/latex/supplement_tex/main_supplement.tex',
        gradient='src/latex/supplement_tex/gradient.tex',
        figures=htmlOutputPath + "/Scripts_PaperFigures_ZZZ_all_PaperFigures.html"
    output: config["FIGDIR"] + '/supplement_final.pdf'
    shell: """
        set -x
        cd src/latex/supplement_tex
        pdflatex main_supplement
        bibtex main_supplement
        pdflatex main_supplement
        pdflatex main_supplement
        cp main_supplement.pdf '../../../{output}'
        pwd
        gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE \
            -dQUIET -dBATCH -sOutputFile='../../../{output}.reduced.pdf' '../../../{output}'
    """
