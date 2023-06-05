"
## snRNA analysis - Individual seurat objects for each cluster
## Author: Maria Tsalkitzidou
## Date: 13/01/2023
## Updated: 17/01/2023
"

###############################################################################################
#### Load the necessary packages and user defined variables ####

rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts") ## Set the working directory

## Load the directories and necessary packages
source("Directories_Packages.R")

#-------------------------------------------------------------------------------------
# User defined variables. 
dataset = "" #The Seurat object will start with this phrase/word after the species information

obj.dir = "" #Seurat object to load

species = "" #rat or human (all lowercase!!)

celltypes.to.subset <- c() #provide the name of the clusters to subset as individual objects e.g. c("Astrocyte", "GPC")

###############################################################################################
#### Load the Seurat object ####

seurat.obj <- readRDS(paste0(seurat.dir, obj.dir))


###############################################################################################
#### Subset the oblect ####

for (cluster in celltypes.to.subset){
  
  #subset
  s.obj <- subset(seurat.obj, idents = cluster)
  
  ## Save the seurat object
  saveRDS(s.obj, file = paste0(seurat.dir, species, cluster, dataset, "_obj.rds"))
}

