# Analysis code for Ing-Simmons et al. 2020

This repository contains analysis code associated with Ing-Simmons et al. 2020, '[Independence of 3D chromatin conformation and gene regulation during Drosophila dorsoventral patterning](https://www.biorxiv.org/content/10.1101/2020.07.07.186791v1)'.

## Overview

This analysis uses data from the following ArrayExpress accessions (which will be publicly available upon publication):
* ChIP-seq: E-MTAB-9303
* scRNA-seq: E-MTAB-9304
* Hi-C: E-MTAB-9306
* Micro-C: E-MTAB-9784

We also used publicly available ChIP-seq and RNA-seq data from [Koenecke et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1057-2), with GEO accession GSE68983, as well as ChIP-seq data from [modENCODE](http://www.modencode.org/).

The analysis is set up as two separate [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows, one for the scRNA-seq analysis and one for all other analysis. 

### ChIP-seq and Hi-C analysis

This analysis is split over multiple `.smk` files for better organisation. These can be found in the `workflow/` directory. Hi-C analysis was carried out using [FAN-C](https://fan-c.readthedocs.io/en/latest/). Both Python and R were used for analysis; individual analysis and plotting scripts used in different steps of the workflow can be found in the `scripts/` directory. 

The workflow assumes that .fastq files for Hi-C and Micro-C are present in `data/fastq/`. The relationship between samples and files is described in `metadata.txt` (and `micro-c_metadata.txt`.

### scRNA-seq analysis

This analysis can be found in the `scrna-seq/` directory. The data was analysed using [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) as shown in the `Snakefile`, with downstream analysis and visualisation in R as shown in the `.Rmd` file in the `scripts/` directory. 
 
