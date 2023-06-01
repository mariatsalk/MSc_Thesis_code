"
## snRNA analysis - QC analysis
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
#### Load the necessary packages and user defined variables ####

## Set the working directory
rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts")

## Load the directories and necessary packages
source("Directories_Packages.R")

#-------------------------------------------------------------------------------------
# User defined variables. 
dataset = "vMB_hVM_merged.seurat.all.20210920" #The Seurat object will start with this phrase/word

obj.dir = "vMB_hVM_merged.seurat.all.20210920_raw_merged_object.rds" #Seurat object to load

species = "human" #rat or human (all lowercase!!)
species_prefix = "premRNAGRCH38-" #premRNARnor6-- or premRNAGRCH38- for rat and human respectively


###############################################################################################
#### Load the Seurat object ####

merged.obj <- readRDS(paste0(seurat.dir, obj.dir))

###############################################################################################
#### STEP 3: Subset the dataset to extract the cells of the desired species and remove prefix ####
seurat.obj <- subset(merged.obj, subset = Species == species)

rownames(seurat.obj)

## Remove the human genes from the dataset
# Get the counts from the matrix ans store it an a variable called counts
counts <- GetAssayData(seurat.obj, slot = "counts")

# Filter the counts and keep only the rat genes. Set value=TRUE, if FALSE a vector containing the (integer) indices of the matches determined by grep is returned, and if TRUE, a vector containing the matching elements themselves is returned.
counts <- counts[grep("^(premRNAGRCH38-)?[A-Z]+", rownames(counts), value = TRUE),]

# Remove prefix from gene names if it exists
for (gene in rownames(counts)){
  if (startsWith(species_prefix, gene)){
    gene <- gsub(species_prefix, "", gene)
  }
}


rownames(counts) #Checkpoint

# Extract the meta data information and the active.ident information from the initial object and store it in a variable
new.meta.data <- seurat.obj@meta.data
new.active.ident <- seurat.obj@active.ident

# Create a new seurat object with the new counts and the meta data information from the initial object and rename the active.idents to match how they were named in the previous object
sn.str.obj <- CreateSeuratObject(counts, meta.data = new.meta.data, assay = "RNA")
Idents(sn.str.obj) <- new.active.ident

rownames(sn.str.obj) #check that the prefix is removed

# Save the modified object
saveRDS(sn.str.obj, file = paste0(seurat.dir, species, "_raw_", dataset, ".rds"))


################################################################################################

#### STEP 4: Extraction of basic information before filtering ####
nCell <- table(sn.str.obj$location, sn.str.obj$Age) %>%
  as.data.frame()


write.xlsx(nCell, file = paste0(plot.dir, qc.dir, species, "_", dataset, "_nCells.xlsx")) 

################################################################################################
#### STEP 5: QC and filtering ####

#calculate the mitochondrial and ribosomal QC metrics for the dataset
sn.str.obj$percent.mt <- PercentageFeatureSet(sn.str.obj, pattern = "M[T|t]-.+")
sn.str.obj$percent.rb <- PercentageFeatureSet(sn.str.obj, pattern = "^R[P|p][SL|sl]")
#tail(grep("^(premRNARnor6--)?[A-Z]{1}[a-z]+", rownames(seurat.obj), value = T, ignore.case = F)) #checkpoint
#visualize the QC metrics with a violin plot
violinplotFeatCount <- VlnPlot(sn.str.obj, features = c("nFeature_RNA", "nCount_RNA"), pt.size=0)
vplotFeatCount <- list()
for (i in seq_along(violinplotFeatCount)){
  vplotFeatCount[[i]] = violinplotFeatCount[[i]] + theme(axis.text.x = element_text(size = 5))
}

violin.plotFeatCount <- plot_grid(plotlist=vplotFeatCount, ncol=2)

violinplotMtRb <- VlnPlot(sn.str.obj, features = c("percent.mt", "percent.rb"), pt.size=0)
vplotMtRb <- list()
for (i in seq_along(violinplotMtRb)){
  vplotMtRb[[i]] = violinplotMtRb[[i]] + theme(axis.text.x = element_text(size = 5))
}

violin.plotMtRb <- plot_grid(plotlist=vplotMtRb, ncol=2)

png(filename=paste0(plot.dir, qc.dir, species, "_", dataset, "_QC_BeforeFiltering_Feat_Count_ViolinPlot.png"), width = 1400, height = 1000)
print(violin.plotFeatCount)
dev.off()

png(filename=paste0(plot.dir, qc.dir, species, "_", dataset, "_QC_BeforeFiltering_MT_RB_ViolinPlot.png"), width = 1400, height = 1000)
print(violin.plotMtRb)
dev.off()

#visualize feature-feature relationships
featurescatter.plot1 <- FeatureScatter(sn.str.obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
featurescatter.plot2 <- FeatureScatter(sn.str.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() 
featurescatter.plot3 <- FeatureScatter(sn.str.obj, feature1 = "nCount_RNA", feature2 = "percent.rb") + theme(legend.text = element_text(size=10))
featurescatter.plot <- featurescatter.plot1 + featurescatter.plot2 + featurescatter.plot3
featurescatter.plot

png(filename=paste0(plot.dir, qc.dir, species, "_", dataset, "_QC_BeforeFiltering_FeatureScatterPlot.png"), width = 1700, height = 1000)
print(featurescatter.plot)
dev.off()


#filter based on the QC metrics
sn.str.obj <- subset(sn.str.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 500 & percent.mt < 3 & percent.rb <3)
table(sn.str.obj$orig.ident)

## Remove samples with less than 100 cells
nCells_table <- table(sn.str.obj$orig.ident)
nCells_df <- as.data.frame(nCells_table)
levels(sn.str.obj@active.ident)

samples_to_remove <- c()
for (i in 1:length(levels(sn.str.obj@active.ident))){
  
  if (nCells_df$Freq[i] < 100){
    samples_to_remove <- append(samples_to_remove, as.character(nCells_df$Var1[i]))
    
  }
}

#Subset the seurat object and remove the samples with less than 100 cells
for (s in samples_to_remove){
  sn.str.obj <- subset(sn.str.obj, subset = orig.ident != s)
}

# Check that the low quality samples are removed
table(sn.str.obj@active.ident)
table(sn.str.obj$orig.ident)

#save the merged object so far
saveRDS(sn.str.obj, file = paste0(seurat.dir, species, "_", dataset, "_QCfiltered.rds"))
#sn.str.obj <- readRDS(file = paste0(seurat.dir, dataset, "_QCfiltered.rds"))





## QC plots after filtering
#visualize the QC metrics with a violin plot
violinplotFeatCount <- VlnPlot(sn.str.obj, features = c("nFeature_RNA", "nCount_RNA"), pt.size=0)
vplotFeatCount <- list()
for (i in seq_along(violinplotFeatCount)){
  vplotFeatCount[[i]] = violinplotFeatCount[[i]] + theme(axis.text.x = element_text(size = 5))
}

violin.plotFeatCount <- plot_grid(plotlist=vplotFeatCount, ncol=2)

violinplotMtRb <- VlnPlot(sn.str.obj, features = c("percent.mt", "percent.rb"), pt.size=0)
vplotMtRb <- list()
for (i in seq_along(violinplotMtRb)){
  vplotMtRb[[i]] = violinplotMtRb[[i]] + theme(axis.text.x = element_text(size = 5))
}

violin.plotMtRb <- plot_grid(plotlist=vplotMtRb, ncol=2)

png(filename=paste0(plot.dir, qc.dir, species, "_", dataset, "_QC_AfterFiltering_Feat_Count_ViolinPlot.png"), width = 1400, height = 1000)
print(violin.plotFeatCount)
dev.off()

png(filename=paste0(plot.dir, qc.dir, species, "_", dataset, "_QC_AfterFiltering_MT_RB_ViolinPlot.png"), width = 1400, height = 1000)
print(violin.plotMtRb)
dev.off()

#visualize feature-feature relationships
featurescatter.plot1 <- FeatureScatter(sn.str.obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
featurescatter.plot2 <- FeatureScatter(sn.str.obj, feature1 = "nCount_RNA", feature2 = "percent.rb") + NoLegend()
featurescatter.plot3 <- FeatureScatter(sn.str.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.text = element_text(size=10))
featurescatter.plot <- featurescatter.plot1 + featurescatter.plot2 + featurescatter.plot3
featurescatter.plot

png(filename=paste0(plot.dir, qc.dir, species, dataset, "_QC_AfterFiltering_FeatureScatterPlot.png"), width = 1700, height = 1000)
print(featurescatter.plot)
dev.off()

#Extract basic information
filtered.nCell <- table(sn.str.obj$location, sn.str.obj$Age) %>%
  as.data.frame()


write.xlsx(filtered.nCell, file = paste0(plot.dir, qc.dir, species, "_", dataset, "_filtered_nCells.xlsx"))


