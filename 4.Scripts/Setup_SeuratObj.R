"
## scRNA analysis - 10x reads to Seurat Objects
## Author: Maria Tsalkitzidou
## Created: 09/11/2022
## Updated: 14/02/2023

Description:



Procedure:



Limitations:
  

"

##########################################################################################
# STEP 1: Load the necessary packages and user defined variables


#-------------------------------------------------------------------------------------
## Set the working directory
rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Desktop/glial_snRNAseq_analysis/4.Scripts")

## Load the directories and necessary packages
source("Directories_Packages.R")

#-------------------------------------------------------------------------------------
# User defined variables. 
dataset = "vMB_hVM" #The Seurat object will start with this phrase/word

#The names of the samples that will become a merged Seurat object
samples =  c("hVM_high_16C_1","hVM_high_32C_1","hVM_low_10A_1","hVM_low_30A_1","AF030_vMB_only_30_core_barcoded","AF030_vMB_only_30_edge_1_barcoded","AF030_vMB_only_30_edge_2_barcoded")

#c("vMB_high_16_A_1", "vMB_high_16_A_2", "vMB_vFB_7_C_1", "vMB_vFB_7_C_2", "vMB_vFB_7_C_3", "vMB_vFB_7_C_4", "vMB_vFB_7_C_5", "vMB_vFB_7_C_6")

#c("hVM_high_16C_1","hVM_high_32C_1","hVM_low_10A_1","hVM_low_30A_1","AF030_vMB_only_30_core_barcoded","AF030_vMB_only_30_edge_1_barcoded","AF030_vMB_only_30_edge_2_barcoded") 


##########################################################################################
# STEP 2: Load the 10x reads and convert them to Seurat objects

#-------------------------------------------------------------------------------------
## Load the samples
for (file in samples){
  seurat.data <- Read10X(data.dir = paste0(data.dire, file, "/"))
  seurat.obj <- CreateSeuratObject(counts = seurat.data,
                                   min.cells = 3,
                                   min.features = 200, project = file)
  
  assign(file, seurat.obj)
}


#-------------------------------------------------------------------------------------
## Create one merged Seurat object
merged.obj <- merge(x = hVM_high_16C_1, y = c(hVM_high_32C_1, hVM_low_10A_1,hVM_low_30A_1,AF030_vMB_only_30_core_barcoded,AF030_vMB_only_30_edge_1_barcoded,AF030_vMB_only_30_edge_2_barcoded),add.cell.id = c("16C1", "32C1", "10A1", "30A1", "AF030", "AF0301", "AF0302"))

table(merged.obj$orig.ident)
#merge(x = vMB_high_16_A_1,y = c(vMB_high_16_A_2, vMB_vFB_7_C_1, vMB_vFB_7_C_2, vMB_vFB_7_C_3, vMB_vFB_7_C_4, vMB_vFB_7_C_5, vMB_vFB_7_C_6), add.cell.ids = c("16A1", "16A2", "7C1", "7C2", "7C3", "7C4", "7C5", "7C6"))

#merge(x = hVM_high_16C_1, y = c(hVM_high_32C_1, hVM_low_10A_1,hVM_low_30A_1,AF030_vMB_only_30_core_barcoded,AF030_vMB_only_30_edge_1_barcoded,AF030_vMB_only_30_edge_2_barcoded),add.cell.id = c("16C1", "32C1", "10A1", "30A1", "AF030", "AF0301", "AF0302"))

## Add the missing information (species, brain region, age of graft) as metadata to the seurat object
merged.obj$location <- "Striatum"
merged.obj$Age <- '6m'


merged.obj$Human <- PercentageFeatureSet(merged.obj, pattern = "^(premRNAGRCH38-)?[A-Z]+")
merged.obj$Rat <- PercentageFeatureSet(merged.obj, pattern = "^(premRNARnor6--)?[A-Z]{1}[a-z]+")
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


# Extract basic information
nCell <- table(merged.obj$location, merged.obj$Age) %>%
  as.data.frame()


write.xlsx(nCell, file = paste0(plot.dir, qc.dir, dataset, "_nCells.xlsx"))

#-------------------------------------------------------------------------------------
#save the raw merged object
saveRDS(merged.obj, file = paste0(seurat.dir, dataset, "_raw_merged_object.rds"))


