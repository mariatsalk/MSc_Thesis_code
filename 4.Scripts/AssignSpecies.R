"
## snRNA analysis - Assign Species
## Author: Maria Tsalkitzidou
## Created: 08/09/2022
## Updated: 24/10/2022

Description:
  The script takes as input graft samples from scRNA (or snRNA) sequencing,. It creates a seuray object, filters the seurat object and assigns the species (rat or human cells) as meta data information in the object and lastly splits the object in two new objects based on the species.


Procedure:



Limitations:
  1) Doesn't take input from the terminal
  2) Isn't 100% generic
  

"

#### STEP 1: Load the necessary packages and user defined variables ####

## Set the working directory
rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Desktop/glial_snRNAseq_analysis/4.Scripts")

## Load the directories and necessary packages
source("Directories_Packages.R")
source("snRNA_Functions.R")

#-------------------------------------------------------------------------------------
# User defined variables. 
dataset = "vMB_vFB_merged.seurat.all.20210920" #The Seurat object will start with this phrase/word

obj.dir = "vMB_vFB_merged.seurat.all.20210920_raw_merged_object.rds" #Seurat object to load

species = "rat" #rat or human


#### STEP 2: Load the Seurat object ####

merged.obj <- readRDS(paste0(seurat.dir, obj.dir))



#### STEP 3: Add the species information as metadata to the seurat object ####


merged.obj <- species_assignment(merged.obj, species.1 = )

merged.obj$Human <- PercentageFeatureSet(merged.obj, pattern = "^premRNAGRCH38-.+")
merged.obj$Rat <- PercentageFeatureSet(merged.obj, pattern = "^premRNARnor6-.+")
merged.obj$Species <- NA

for (i in 1:length(merged.obj$Species)){
  if (merged.obj$Human[i] > 80 && merged.obj$Rat[i] < 10){
    merged.obj$Species[i] <- "human"
  } else if (merged.obj$Rat[i] > 80 && merged.obj$Human[i] < 10){
    merged.obj$Species[i] <- "rat"
  }
  else{
    merged.obj$Species[i] <- "Undefined"
  }
}



#### STEP 4: Subset the dataset to extract the cells of the desired species ####
seurat.obj <- subset(merged.obj, subset = Species == species)




#### Convert to SCE and remove prefix from gene names ####
seurat.sce <- as.SingleCellExperiment(seurat.obj, assay = "RNA")
rownames(seurat.sce) <- gsub("premRNAGRCH38-", "", rownames(seurat.sce))

# Convert back to seurat object
mod.seurat.obj <- as.Seurat(seurat.sce)





