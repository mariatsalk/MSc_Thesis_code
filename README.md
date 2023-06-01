# M.Sc. Thesis: Host and graft glial composition in a stem cell transplantation model of Parkinson's disease in rats.

This repository consists of a collection of R scripts that were used for analysing the glial composition 

As some of the scripts are re-used with different data they need user input (e.g. resolution for UMAP), so running them from the command line is not recommended.

## Research Background
Thesis abstract



## Tools/Software

All the analysis was performed using R(v4.2.2) and the IDE of choice was Rstudio(v2022.07.02+576). The R packages used were:

- devtools
- BiocManager(v1.30.20)
- Seurat(v4.3.0)
- dplyr(v1.1.0)
- patchwork(v1.1.2)
- harmony(v0.0.1)
- rliger
- ggplot2(v3.4.2)
- stringr(v1.5.0)
- openxlsx(v4.2.5.2)
- tidyverse(v2.0.0)
- presto(v1.0.0)
- msigdbr()
- fgsea(v1.24.0)
- enrichplotv(v1.18.4)


## Workflow
The file paths used are specific to my machine so they need to be edited for the code to run properly.

There is a number of R scripts that are not run as part of the workflow but are sourced in the R scripts:

- Directories_Packages.R - The main paths and variables are defined in the script and should be modified there.
- snRNA_Functions.R - Contains functions used in the workflow.
- Silhouette scores.R - Calculates the mathematically optimal value for resolution for seurat clustering.



### Pre-process of data
As this is a meta-analysis of snRNA-seq data the datasets were provided in different formats, either as a seurat object or a gene matrix so pre-processing was necessary in order to merge the datasets correctly.


The scripts should be run in the following order:
1. Setup_SeuratObj.R
2. Merge_SeuratObj.R


### Analysis of host and graft environment

After merging the data, two new datasets were created; one for the host environment and one for the graft environment.

Initially the data were filtered for low quality cells and then normalised and integrated with HARMONY. Next, the distinct cell types were identified by comparing the highly expressed genes per cluster to well-known canonical markers for cell types in substantia nigra and striatum. After that, the cell types were assigned and some plots and information were extracted. Finally, the glial cell populations of each environment were subset to idnividual seurat object for downstream analysis.

The analysis followed the same steps for the host and graft environment, which are: 

1. QC_Filtering.R
2. NormaliseIntegrate.R
3. IdentifyCellTypes.R
4. AssignCellTypes.R
5. Plots_Stats.R
6. SplitObjPerCluster.R


### High resolution analysis for glial cell populations

Each glial population identified in the host and graft environments is further analysed with the following sub-workflow. This analysis re-uses some of the previous scripts.

First, the data needed to be normalised and integrated with HARMONY again. Similarly to the host and graft environment analysis, the identification and assignment of cell types withing the glial populations followed after normalisation. Next, the differentially expressed (DE) genes were extracted and gene set enrichment analysis (GSEA) was performed to investigate the functionality of the glial sub-populations.

The scripts used for this should be run as:

1. NormaliseIntegrate.R
2. IdentifyCellTypes.R
3. AssignCellTypes.R
2. DEanalysis.R
3. GO_enrichment_analysis.R
