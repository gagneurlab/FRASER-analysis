# FRASER-analysis

This is the accompanying analysis repository of the paper: `Detection of aberrant splicing events in RNA-seq data with FRASER`.

The paper can be found [here](https://www.biorxiv.org/content/10.1101/2019.12.18.866830v1).

This repository contains the full pipeline and code to reproduce the results published in the paper using [snakemake](https://snakemake.readthedocs.io/en/stable/) and [wBuild](https://github.com/gagneurlab/wBuild). 

## Project structure

This project is setup as a [wBuild workflow](https://github.com/gagneurlab/wBuild). This is an automatic build tool for R reports based on [snakemake](https://snakemake.readthedocs.io/en/stable/).

* The `wbuild.yaml` is the main configuration file to setup up the workflow
* The `Scripts` folder contains scripts which will be rendered as HTML reports
* The `src` folder contains additional helper functions and scripts
* The `Output` folder will contain all files produced in the analysis pipeline
    * `Output/data` has all raw RDS output files
    * `Output/html` contains the final HTML report
    * `Output/paper_figures` has all paper figures

## Dependencies

This project depends on the python package `wBuild` and the R package `FRASER`.
The pipeline starts with the raw aligned GTEx samples V7P and their genotype calls, which can be downloaded from [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v7.p2). Since the data are not publicly shareable one has to apply for the data at [dbGaP]( https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v6.p1). 

## Repository setup

First download the repo and its dependencies:

```
# R package used throughout the workflow
git clone https://github.com/gagneurlab/FRASER FRASER
# analysis code
git clone https://github.com/gagneurlab/FRASER-analysis FRASER-analysis
cd FRASER-analysis
```

and install wbuild using pip by running.

```
pip install wBuild
wBuild init
```

Since `wBuild init` will reset the current `Snakefile`, ` readme.md`, and `wbuild.yaml` we have to revert them again with git.

```
git checkout Snakefile
git checkout wbuild.yaml
git checkout readme.md
```

To make sure all packages needed in the analysis are installed source the following file in R

```
Rscript ./src/r/install_dependencies.R
```

## Run the full pipeline

To run the full pipeline, execute the following command with 10 cores in parallel:

```
snakemake -j 10
```

or to run it on the cluster with SLUM installed: 

```
snakemake -k --restart-times 2 --cluster "sbatch -N 1 -n 10 --mem 80G" --jobs 20
```
