## snRNA analysis - Preproccess seurat object befor CellOracle GRN construction
## Author: Maria Tsalkitzidou
## Date: 24/03/2023
## Updated: 

setwd("~/Desktop/scRNAseq_analysis/4.Scripts/")
#----------------------------------------------------------------------------------------------
## Define the neseccary directories, files and objects
source("Directories_Packages.R")
library(SeuratData)
library(SeuratDisk)


## User defined variables
#obj.dir = "SubstantiaNigra&Striatum_RatCells_RemovedLowQualitySamples&HumanGenes_29112022.rds" #The object to be loaded in the script
obj.dir = "Rat_Astrocyte_subclusters_VarFeat_Scaled_HARMONY_15PCs_res005_obj.rds" #The object to be loaded in the script

#Astrocyte_SN_STR_v2_HARMONY_normalizedscaled_annotated_obj.rds
#Astrocyte_subclusters_VarFeat_Scaled_HARMONY_obj.rds
dataset = "Human" #Name of the dataset to be loaded in the script (in this case: Human, Rat, Merged)

## Read the object
s.obj <- readRDS(file = paste0(seurat.dir, obj.dir))
DimPlot(s.obj)

#----------------------------------------------------------------------------------------------
## SEURAT TO ANNDATA WITH SEURATDISK (https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)

SaveH5Seurat(s.obj, filename = "SaveH5Seurat_Human_Astrocyte_subclusters_VarFeat_Scaled_HARMONY_obj.h5Seurat")
Convert("SaveH5Seurat_Human_Astrocyte_subclusters_VarFeat_Scaled_HARMONY_obj.h5Seurat", dest = "h5ad")




#----------------------------------------------------------------------------------------------
## SEURAT TO ANNDATA WITH PETTER'S CODE


remotes::install_github("pmbio/MuDataSeurat")

# step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
s.obj.slim = DietSeurat(
  s.obj,
  counts = TRUE, # so, raw counts save to adata.layers['counts']
  data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
  features = rownames(s.obj), # export all genes, not just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

MuDataSeurat::WriteH5AD(s.obj.slim, "DietSeurat_Human_Astrocyte_subclusters_VarFeat_Scaled_HARMONY_obj.h5ad", assay="RNA")




#----------------------------------------------------------------------------------------------
## SEURAT TO ANNDATA WITH PETTER'S CODE + my modifications: get the top 2-3K variable genes only

DimPlot(s.obj)

# step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
s.obj.slim2 = DietSeurat(
  s.obj,
  counts = TRUE, # so, raw counts save to adata.layers['counts']
  data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = TRUE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
  features = s.obj@assays$RNA@var.features, # export just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)


# step 2: factor to character, or else your factor will be number in adata 
i <- sapply(s.obj.slim2@meta.data, is.factor)
s.obj.slim2@meta.data[i] <- lapply(s.obj.slim2@meta.data[i], as.character)


# optional step: downsample the seurat object
s.obj.extraslim <- s.obj.slim2[, sample(colnames(s.obj.slim2), size = 20000, replace=F)]
s.obj.extraslim <- s.obj.extraslim[, sample(s.obj.extraslim@assays$RNA@var.features, size = 2000, replace=F)]
ncol(s.obj.extraslim)

# step 3: Save in h5ad format.
MuDataSeurat::WriteH5AD(s.obj.extraslim, "DietSeurat_WriteH5AD_20Kcells_Scaled_Human_SN_STR_v2_HARMONY_normalizedscaled_annotated_obj.h5ad", assay="RNA", overwrite = TRUE)



SaveH5Seurat(s.obj.slim2, filename = "DietSeurat_SaveH5Seurat_Human_Astrocyte_subclusters_20PCs_res005_VarFeat_Scaled_HARMONY_obj.h5Seurat", overwrite = TRUE)
Convert("DietSeurat_SaveH5Seurat_Human_Astrocyte_subclusters_20PCs_res005_VarFeat_Scaled_HARMONY_obj.h5Seurat", dest = "h5ad", overwrite = TRUE)



#----------------------------------------------------------------------------------------------
## SEURAT TO ANNDATA WITH EDO'S CODE

DimPlot(s.obj)

# step 0: Remove prefix and store celltype information as metadata (optional)
# Remove prefix
species_prefix <- "premRNARnor6--"

#Extract the needed information from the "old" seurat object
new.meta.data <- s.obj@meta.data
new.active.ident <- s.obj@active.ident
new.var.features <- s.obj@assays$RNA@var.features

# Extract the counts
counts <- GetAssayData(s.obj, slot = "counts", assay = "RNA")

# Remove the prefix from the counts and the variable features
rownames(counts) <- gsub(species_prefix, "", rownames(counts))
new.var.features <- gsub(species_prefix, "", new.var.features)

#Create the new seurat object and add the extracted information
seurat.obj <- CreateSeuratObject(counts, meta.data = new.meta.data, assay = "RNA")
Idents(seurat.obj) <- new.active.ident
seurat.obj@assays$RNA@var.features <- new.var.features

rownames(seurat.obj) #check that the prefix is removed

#Add the cluster names as metadata information
seurat.obj$celltype <- Idents(seurat.obj)


# step 1: factor to character, or else your factor will be number in adata 
seurat.obj$orig.ident <- as.character(seurat.obj$orig.ident)
seurat.obj$location <- as.character(seurat.obj$location)
seurat.obj$Age <- as.character(seurat.obj$Age)
seurat.obj$celltype <- as.character(seurat.obj$celltype)


# step 2: Slim down a Seurat object. So you get raw counts, lognorm counts
s.obj.slim2 = DietSeurat(
  seurat.obj,
  counts = TRUE, # so, raw counts save to adata.layers['counts']
  data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = TRUE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
  features = seurat.obj@assays$RNA@var.features, # export just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)



# step 3: Save in h5ad format.
SaveH5Seurat(s.obj.slim2, filename = "../CellOracleAnalysis/DietSeurat_SaveH5Seurat_Rat_Astrocyte_HARMONY_normalizedscaled_annotated_obj.h5Seurat", overwrite = TRUE)
Convert("../CellOracleAnalysis/DietSeurat_SaveH5Seurat_Rat_Astrocyte_HARMONY_normalizedscaled_annotated_obj.h5Seurat", dest = "h5ad", overwrite = TRUE)







