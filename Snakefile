import pandas as pd
import pathlib
import wbuild
config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)

configfile: "wbuild.yaml"

htmlOutputPath = config['htmlOutputPath']

# sub-part of enrichment pipeline
include: "Scripts/Leafcutter/leafcutter_snakemake_functions.py"
include: "Scripts/GTExEnrichment/enrichment.snake"

# parse main wbuild script
include: config['wBuildPath'] + "/wBuild.snakefile"

rule all:
	input: rules.Index.output, htmlOutputPath + "/readme.html"
	output: touch("Output/all.done")


rule prepub:
    shell: """
        rsync -Ort {config[htmlOutputPath]} {config[webDir]}
        chmod -R uo+rX {config[webDir]}
    """

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

rule exportData:
    input:
        expand(config["DATADIR"] + "/export/fds_objects/{dataset}__FRASER_FDS.xz.RDS", dataset=list(set(config["EnrichmentTissues"] + config["datasets"])))
    shell: "echo Export done"


