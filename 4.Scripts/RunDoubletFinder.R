"
## snRNA analysis - Run DoubletFinder
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
###############################################################################################
#### STEP 1: Load the necessary packages and user defined variables ####

## Set the working directory
rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts")

## Load the directories and necessary packages
source("Directories_Packages.R")
library(DoubletFinder)
#-------------------------------------------------------------------------------------
# User defined variables. 
dataset = "Rat_SN&STR_" #The Seurat object will start with this phrase/word

obj.dir = "Rat_SN&STR_Filtered&RemovedLowQualitySamples&HumanGenes_normalized&scaled_obj.rds" #Seurat object to load

species = "rat" #rat or human (all lowercase!!)
species_prefix = "premRNARnor6--" #premRNARnor6-- or premRNAGRCH38- for rat and human respectively


###############################################################################################
#### STEP 2: Load the Seurat object ####

sn.str.obj <- readRDS(paste0(seurat.dir, obj.dir))




# Rename the active.ident with the sample names
Idents(sn.str.obj) <- "orig.ident"
################################################################################################
################################################################################################
################################ IGNORE DOESNT WORK #######################################
################################################################################################
################################################################################################
#-----
## Split the merged object into objects per sample to be able to apply doublet finder
obj.list <- c()
for (sample in levels(sn.str.obj@active.ident)){
  #Subset the merged seurat object per orig.ident 
  seurat.obj <- subset(sn.str.obj, idents = sample)
  #assign the name of the file (aka sample) as the name of the object
  assign(sample, seurat.obj)
  
  # Keep only the samples that have more than 10 cells
  if(length(colnames(seurat.obj)) > 10){
    #add the names of the seurat objects in the obj.list
    obj.list <- append(obj.list, seurat.obj)
  }
  
}


#-----
## Pre-process function with DoubletFinder
# Pre-process seurat object (one sample per time)
for (sample.obj in obj.list){
  
  #extract the name of the sample/seurat object
  sample_name <- levels(sample.obj@active.ident)
  
  #assign the name of the file (aka sample) as the name of the object
  assign(sample_name, sample.obj)
  
  sample.obj <- sample.obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = 'vst', nfeatures = 4000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:10)
  
  
  ## pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(sample.obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn.sample <- find.pK(sweep.stats)
  
  #Plot the pK and BCmetric values of each sample
  pKplot <- ggplot(bcmvn.sample, aes(pK, BCmetric, group=1)) +
    geom_point() +
    geom_line()
  
  #store the plot as png image
  png(filename = paste0(plot.dir, qc.dir, sample_name, "_pK_plot.png"), width = 1000, height = 700)
  print(pKplot)
  dev.off()
  
  #Select the optimal pK that corresponds to max bcmvn to optimize doublet detection
  pK.optimal <- bcmvn.sample %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK.optimal <- as.numeric(as.character(pK.optimal[[1]]))
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- sample.obj@meta.data$seurat.clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(sample.obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  
  ## Run DoubletFinder
  sample.obj <- doubletFinder_v3(sample.obj,
                                 PCs = 1:10,
                                 pN = 0.25,
                                 pK = pK.optimal,
                                 nExp = nExp_poi.adj,
                                 reuse.pANN = FALSE,
                                 sct = FALSE)
  
  #extract the DF.classifications column name from the seurat object (the column is named differently for every sample since the name includes the pN and pK value)
  DF.classifications = colnames(sample.obj@meta.data)[grepl(paste0("DF.classifications_0.25_", pK.optimal, '_', nExp_poi.adj), colnames(sample.obj@meta.data))]
  
  #visualise doublets
  doublet.dimplot <- DimPlot(sample.obj, reduction = 'umap', group.by = DF.classifications)
  doublet.dimplot
  png(filename = paste0(plot.dir, qc.dir, sample_name, "_SN&STR_doublet_dimplot.png"), width = 1000, height = 700)
  print(doublet.dimplot)
  dev.off()
  
  # number of singlets and doublets
  doublet.table <- table(sample.obj@meta.data[DF.classifications])
  doublet.df <- as.data.frame(doublet.table)
  
  #Write the output in an excel file
  write.xlsx(doublet.df, file = paste0(plot.dir, qc.dir, sample_name, "_SN&STR_N_Singlets_Doublets.xlsx"), rowNames = FALSE, overwrite = TRUE)
  
  #Subset the object and keep only the singlets
  sample.obj <- sample.obj[, sample.obj@meta.data[, DF.classifications] == "Singlet"]
}

# Merge the filtered samples
doublefree.obj <- merge(x=obj.list[1], y=obj.list)


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
## DOUBLE FINDER WITHOUT SPLITTING THE OBJECT

## Pre-process function with DoubletFinder
# Pre-process seurat object (one sample per time)

  

## pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(sn.str.obj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn.sample <- find.pK(sweep.stats)
  
#Plot the pK and BCmetric values of each sample
pKplot <- ggplot(bcmvn.sample, aes(pK, BCmetric, group=1)) +
  geom_point() +
  geom_line()

#store the plot as png image
png(filename = paste0(plot.dir, qc.dir, dataset, "_pK_plot.png"), width = 1000, height = 700)
print(pKplot)
dev.off()
  
#Select the optimal pK that corresponds to max bcmvn to optimize doublet detection
pK.optimal <- bcmvn.sample %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK.optimal <- as.numeric(as.character(pK.optimal[[1]]))
  
## Homotypic Doublet Proportion Estimate
annotations <- sn.str.obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(sn.str.obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  

## Run DoubletFinder
sn.str.obj <- doubletFinder_v3(sn.str.obj,
                               PCs = 1:10,
                               pN = 0.25,
                               pK = pK.optimal,
                               nExp = nExp_poi.adj,
                               reuse.pANN = FALSE,
                               sct = FALSE)
  
#extract the DF.classifications column name from the seurat object (the column is named differently for every sample since the name includes the pN and pK value)
DF.classifications = colnames(sn.str.obj@meta.data)[grepl(paste0("DF.classifications_0.25_", pK.optimal, '_', nExp_poi.adj), colnames(sn.str.obj@meta.data))]

#visualise doublets
doublet.dimplot <- DimPlot(sn.str.obj, reduction = 'umap', group.by = DF.classifications)
doublet.dimplot
png(filename = paste0(plot.dir, qc.dir, dataset, "_SN&STR_doublet_dimplot.png"), width = 1000, height = 700)
print(doublet.dimplot)
dev.off()

# number of singlets and doublets
doublet.table <- table(sn.str.obj@meta.data[DF.classifications])
doublet.df <- as.data.frame(doublet.table)

#Write the output in an excel file
write.xlsx(doublet.df, file = paste0(plot.dir, qc.dir, dataset, "_SN&STR_N_Singlets_Doublets.xlsx"), rowNames = FALSE, overwrite = TRUE)

#Subset the object and keep only the singlets
doubletfree.obj <- sn.str.obj[, sn.str.obj@meta.data[, DF.classifications] == "Singlet"]

# Save the doublet free object
saveRDS(doubletfree.obj, file = paste0(seurat.dir, dataset, "_DoubletFree_obj.rds"))
