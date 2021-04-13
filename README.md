[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4542753.svg)](https://doi.org/10.5281/zenodo.4542753)

# Analysis code for Ing-Simmons et al. 2021

This repository contains analysis code associated with Ing-Simmons et al., *Nature Genetics*, 2021, '[Independence of chromatin conformation and gene regulation during Drosophila dorsoventral patterning](https://www.nature.com/articles/s41588-021-00799-x)'.

This code is provided as an enhancement to the Methods described in the paper, in order to increase the reproducibility of our work and help interested readers understand how the analysis was carried out. Please note that while we aim for easy reproducibility, this repository contains the actual code that was developed over time and used for the analysis and may not be a fully self-contained workflow to reproduce all results in the paper. 

## Overview
The analysis is set up as several separate [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows, for the initial ChIP-seq analysis, for the scRNA-seq analysis, and one for all other analysis and integration of different datasets. 

### Data

This analysis uses data from the following ArrayExpress accessions:
* ChIP-seq: E-MTAB-9303
* scRNA-seq: E-MTAB-9304
* Hi-C: E-MTAB-9306
* Micro-C: E-MTAB-9784

We also used publicly available ChIP-seq and RNA-seq data from [Koenecke et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1057-2), with GEO accession GSE68983, as well as ChIP-seq data from [modENCODE](http://www.modencode.org/).

### Software requirements

Python packages and versions used can be found in `requirements.txt`. R packages and versions used can be found in `renv.lock` (for more information on using renv, see [here](https://rstudio.github.io/renv/articles/renv.html)).

### ChIP-seq and RNA-seq analysis

This can be found in the `external_data/` directory. Each individual data source has its own subdirectory with its own `Snakefile`. These are set up so that the relationship between files and samples is specified in a `metadata.txt` file. Only the input fastq files and the `metadata.txt` files should be needed to run the workflow to produce `.bam` and `.bw` files. 

### scRNA-seq analysis

This analysis can be found in the `scrna-seq/` directory. The data was analysed using [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) as shown in the `Snakefile`, with downstream analysis and visualisation in R as shown in the `.Rmd` file in the `scripts/` directory. 

### Hi-C analysis and integrative analysis

This analysis is split over multiple `.smk` files for better organisation. These can be found in the `workflow/` directory. This workflow includes identification of putative tissue-specific enhancers and their downstream analysis, as well as all Hi-C and Micro-C analysis, integrative analysis, and visualisation.

Hi-C analysis was carried out using [FAN-C](https://fan-c.readthedocs.io/en/latest/). Both Python and R were used for analysis; individual analysis and plotting scripts used in different steps of the workflow can be found in the `scripts/` directory. 

The workflow assumes that .fastq files for Hi-C and Micro-C are present in `data/fastq/`. The relationship between samples and files is described in `metadata.txt` (and `micro-c_metadata.txt`).

