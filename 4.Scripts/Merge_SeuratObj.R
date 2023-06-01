"
## scRNA analysis - Merge Seurat Objects
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


##########################################################################################
# STEP 2: Load the Seurat objects to merge

seurat.obj1 <- readRDS(paste0(seurat.dir, "merged.seurat.all.20210920.rds"))
seurat.obj2 <- readRDS(paste0(seurat.dir, "vMB_hVM_raw_merged_object.rds"))
seurat.obj3 <- readRDS(paste0(seurat.dir, "vMB_vFB_raw_merged_object.rds"))

#-------------------------------------------------------------------------------------
## Merge the objects

rat.obj <- merge(x = seurat.obj1, y = seurat.obj3) #the dataset for the host environment analysis

human.obj <- merge(x = seurat.obj1, y = seurat.obj2) #the dataset for the graft cells analysis



##########################################################################################
# STEP 3: Get some basic info from the each dataset

## Number of Cells per brain region per time point
nCell.rat.obj <- table(rat.obj$location, rat.obj$Age) %>%
  as.data.frame()


nCell.human.obj <- table(human.obj$location, human.obj$Age) %>%
  as.data.frame()


write.xlsx(nCell.rat.obj, file = paste0(plot.dir, qc.dir, "nCells_perBrainArea_perAge_vMB_vFB_merged.seurat.all.20210920.xlsx"))

write.xlsx(nCell.human.obj, file = paste0(plot.dir, qc.dir, "nCells_perBrainArea_perAge_vMB_hVM_merged.seurat.all.20210920.xlsx"))


#-------------------------------------------------------------------------------------
#save the raw merged object
saveRDS(rat.obj, file = paste0(seurat.dir, "vMB_vFB_merged.seurat.all.20210920_raw_merged_object.rds"))

saveRDS(human.obj, file = paste0(seurat.dir, "vMB_hVM_merged.seurat.all.20210920_raw_merged_object.rds"))




rownames(seurat.obj1)
raw_seurat <- as.data.frame(table(seurat.obj1$Species, seurat.obj1$orig.ident))
species.table <- as.data.frame(table(seurat.obj1$Species, seurat.obj1$rat_counts, seurat.obj1$human_counts))


rat.obj <- subset(seurat.obj1, subset = Species == 'rat')
table(rat.obj$Species, rat.obj$orig.ident)
rownames(rat.obj)

merged.obj <- seurat.obj1
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

renamed_seurat <- as.data.frame(table(merged.obj$Species, merged.obj$orig.ident))


test.obj <- subset(merged.obj, subset = Species == 'Undefined')
tail(rownames(test.obj))

grep("^(premRNARnor6--)?[A-Z]{1}[a-z]+", rownames(merged.obj), value = T, ignore.case = F)
tail(grep("^(premRNAGRCH38-)?[A-Z]+", rownames(merged.obj), value = T, ignore.case = F))

saveRDS(merged.obj, file = paste0(seurat.dir, "RENAMED_merged.seurat.all.20210920.rds"))
