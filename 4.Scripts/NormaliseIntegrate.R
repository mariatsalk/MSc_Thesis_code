"
## snRNA analysis - Normalisation & Integration
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
#### Load the necessary packages and user defined variables ####

## Set the working directory
rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts")

## Load the directories and necessary packages
source("Directories_Packages.R")

#-------------------------------------------------------------------------------------
# User defined variables. 
dataset = "SN&STR_" #The Seurat object will start with this phrase/word

obj.dir = "human_SN&STR_v2_RemovedLowQualitySamples&prefix_03032023.rds" #Seurat object to load

species = "rat" #rat or human (all lowercase!!)
#species_prefix = "premRNAGRCH38-" #premRNARnor6-- or premRNAGRCH38- for rat and human respectively


###############################################################################################
#### Load the Seurat object ####

s.obj <- readRDS(paste0(seurat.dir, obj.dir))


###############################################################################################
#### Normalize and Scale the seurat object ####

s.obj <- s.obj%>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = counts@var.genes, npcs = 50, verbose = FALSE) %>%
  JackStraw(num.replicate = 100) %>%
  ScoreJackStraw(dims=1:20)

#Determine the dimensionality of the dataset
jackstrawplot <- JackStrawPlot(s.obj)
elbowplot <- ElbowPlot(s.obj)

#Save the JackStrawPlot and Elbow plot
png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_JackStrawplot.png"), width = 1000, height = 700)
print(jackstrawplot)
dev.off()

png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_Elbowplot.png"), width = 1000, height = 700)
print(elbowplot)
dev.off()

###############################################################################################
#### Find Neighbors and Clusters with no integration ####

s.obj1 <- s.obj %>% 
  RunUMAP(reduction = "pca", dims = 1:50) %>% 
  FindNeighbors(reduction = "pca", dims = 1:50) %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()

## Plot the UMAP
dimplot.obj1 <- DimPlot(s.obj1, label = F)
dimplot.obj1.location <- DimPlot(s.obj1, group.by = "location", shuffle = T)
dimplot.obj1.age <- DimPlot(s.obj1, group.by = "Age", shuffle = T)

#Save the UMAP
png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_UMAP.png"), width = 1000, height = 700)
print(dimplot.obj1)
dev.off()

png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_perBrainArea_UMAP.png"), width = 1000, height = 700)
print(dimplot.obj1.location)
dev.off()

png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_perAge_UMAP.png"), width = 1000, height = 700)
print(dimplot.obj1.age)
dev.off()


###############################################################################################
#### Integrate with Harmony and find Neighbors and Clusters ####

s.obj2 <- RunHarmony(s.obj, group.by.vars = "orig.ident")

## Repeat the necessary steps with harmony embeddings
s.obj2 <- s.obj2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.05) %>% 
  identity()

## Plot the UMAP with Harmony integration
dimplot.obj2 <- DimPlot(s.obj2, label = F)
dimplot.obj2.location <- DimPlot(s.obj2, group.by = "location", shuffle = T)
dimplot.obj2.age <- DimPlot(s.obj2, group.by = "Age", shuffle = T)

#Save the UMAPs with HARMONY integration
png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_HARMONY_UMAP.png"), width = 1000, height = 700)
print(dimplot.obj2)
dev.off()

png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_HARMONY_perBrainArea_UMAP.png"), width = 1000, height = 700)
print(dimplot.obj2.location)
dev.off()

png(filename=paste0(plot.dir, species, "/", norm.dir, dataset, "_HARMONY_perAge_UMAP.png"), width = 1000, height = 700)
print(dimplot.obj2.age)
dev.off()



#----------------------------------------------------------------------------------------------
## Create dataframes with the cell number per cluster
s.obj1.table <- table(s.obj1@active.ident)
s.obj1.df <- as.data.frame(s.obj1.table)

s.obj2.table <- table(s.obj2@active.ident)
s.obj2.df <- as.data.frame(s.obj2.table)


#Save the dataframes with the cell number per cluster as excel files
write.xlsx(s.obj1.df, file = paste0(plot.dir, species, "/", norm.dir, dataset, "_CellsPerCluster_NoIntegration.xlsx"), rowNames = FALSE, sheetName = "No Integration")
write.xlsx(s.obj2.df, file = paste0(plot.dir, species, "/", norm.dir, dataset, "_CellsPerCluster_HARMONY.xlsx"), rowNames = FALSE, sheetName = "HARMONY", append = TRUE)



###############################################################################################
#### Save the objects ####

saveRDS(s.obj1, file=paste0(seurat.dir, species, "_", dataset, "_normalized&scaled_obj.rds"))
saveRDS(s.obj2, file=paste0(seurat.dir, species, "_", dataset, "_HARMONY_normalized&scaled_obj.rds"))

