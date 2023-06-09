"
## snRNA analysis - Assign Cell Types
## Author: Maria Tsalkitzidou
## Created: 08/09/2022
## Updated: 31/05/2023
"
############################################################################################### 
#### Load the necessary packages and user defined variables ####

rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts") ## Set the working directory


#-------------------------------------------------------------------------------------
## Load the directories and necessary packages
source("Directories_Packages.R")

#-------------------------------------------------------------------------------------
# User defined variables.

obj.dir = "" #The object to be loaded in the script
species = "" #Name of the dataset to be loaded in the script (in this case: Human, Rat, Merged)

colour.palette <- c() #colours for the UMAP per cell type

celltype.names <- c() #The cell types identified


###############################################################################################
#### Load the Seurat object ####
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))

#DimPlot(seurat.obj) #checkpoint

###############################################################################################
#### Rename the clusters by their celltype ####

names(celltype.names) <- levels(seurat.obj)
annotated.seurat.obj <- RenameIdents(seurat.obj, celltype.names)
annotated.seurat.obj$celltype <- Idents(annotated.seurat.obj)



###############################################################################################
#### Produce UMAP plots ####

#UMAP per brain region
UMAP_BrainRegion <- DimPlot(annotated.seurat.obj, reduction = "umap", group.by = "location", label.size = 5, shuffle = T) 
UMAP_BrainRegion


#UMAP per Cell Type
UMAP_CellType <- DimPlot(annotated.seurat.obj, reduction = "umap", label = F, label.size = 5, cols = coulour.palette)
UMAP_CellType


#UMAP combined
UMAP_Combined <- UMAP_BrainRegion + UMAP_CellType
UMAP_Combined



#save the UMAPs in png images
png(filename=paste0(plot.dir, dataset, "_UMAP_ManuallyAnnotated_PerCellType.png"), width = 1200, height = 800)
print(UMAP_CellType)
dev.off()

png(filename=paste0(plot.dir, dataset, "_UMAP_annotated_PerSpecies.png"), width = 1200, height = 800)
print(UMAP_BrainRegion)
dev.off()

png(filename=paste0(plot.dir, dataset, "_UMAP_annotated_perSpecies_Combined.png"), width = 1200, height = 800)
print(UMAP_Combined)
dev.off()



###############################################################################################
#### Save the annotated seurat obj ####

saveRDS(annotated.seurat.obj, file = paste0(seurat.dir, dataset, "_annotated_obj.rds"))


