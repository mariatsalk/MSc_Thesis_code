
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("Biocmanager")
}

remotes::install_github("mojaveazure/seurat-disk")

devtools::install_github('satijalab/seurat-data')

install_github('immunogenomics/presto')
BiocManager::install("fgsea")