"
## scRNA analysis - Merge Seurat Objects
## Author: Maria Tsalkitzidou
## Created: 09/11/2022
## Updated: 14/02/2023
"

##########################################################################################
#### Load the necessary packages and user defined variables ####

rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Desktop/glial_snRNAseq_analysis/4.Scripts") ## Set the working directory

## Load the directories and necessary packages
source("Directories_Packages.R")


##########################################################################################
#### Load the Seurat objects to merge ####

seurat.obj1 <- readRDS(paste0(seurat.dir, "merged.seurat.all.20210920.rds"))
seurat.obj2 <- readRDS(paste0(seurat.dir, "vMB_hVM_raw_merged_object.rds"))
seurat.obj3 <- readRDS(paste0(seurat.dir, "vMB_vFB_raw_merged_object.rds"))

#-------------------------------------------------------------------------------------
## Merge the objects

rat.obj <- merge(x = seurat.obj1, y = seurat.obj3) #the dataset for the host environment analysis

human.obj <- merge(x = seurat.obj1, y = seurat.obj2) #the dataset for the graft cells analysis



##########################################################################################
#### Get some basic info from the each dataset ####

## Number of Cells per brain region per time point
nCell.rat.obj <- table(rat.obj$location, rat.obj$Age) %>%
  as.data.frame()


nCell.human.obj <- table(human.obj$location, human.obj$Age) %>%
  as.data.frame()


write.xlsx(nCell.rat.obj, file = paste0(plot.dir, qc.dir, "nCells_perBrainArea_perAge_vMB_vFB_merged.seurat.all.20210920.xlsx"))

write.xlsx(nCell.human.obj, file = paste0(plot.dir, qc.dir, "nCells_perBrainArea_perAge_vMB_hVM_merged.seurat.all.20210920.xlsx"))


##########################################################################################
#### Save the raw merged object ####
saveRDS(rat.obj, file = paste0(seurat.dir, "vMB_vFB_merged.seurat.all.20210920_raw_merged_object.rds"))

saveRDS(human.obj, file = paste0(seurat.dir, "vMB_hVM_merged.seurat.all.20210920_raw_merged_object.rds"))
