 "
## snRNA analysis - Identify cell types by signature markers
## Author: Maria Tsalkitzidou
## Date: 11/10/2022
## Updated: 15/02/2023

Description:
  The script takes as input a seurat object and an R script with cell signature markers as vectors and plots the expression of the markers in violin and feature plots per marker.
  As the file that contains the markers might not be applicable to our organism there is function that customizes the markers to our dataset (Human: premRNAGRCH38-, Rat: premRNARnor6--)
  This script is customized to check for specific markers from the signature markers file and produce feature plots of the highest expressed markers per cell type


Procedure:



Limitations:
  1) Doesn't take input from the terminal
  2) Isn't 100% generic
  3) If the merged dataset is loaded (both rat and human cells) it still finds and plots differentially expressed features for one of the two species

"
 
#### Load the necessary packages and user defined variables ####
 
## Set the working directory
rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts")
 
## Load the directories and necessary packages
source("Directories_Packages.R")
source("Markers.R")
source("snRNA_functions.R")

#-------------------------------------------------------------------------------------------------------
## User defined variables
obj.dir = "Human_SN&STR_v2_HARMONY_normalized&scaled_obj.rds" #The object to be loaded in the script
species = "Human" #Name of the dataset to be loaded in the script (in this case: Human, Rat, Merged)
prefix = NULL #prefix in front of the marker (in this case: "premRNAGRCH38-", "premRNARnor6--" or NULL if the prefix is removed)



#--------------------------------------------------------------------------------------------------------
## Load the data
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))

markers <- read_excel("../Markers.xlsx", sheet = 1)

#--------------------------------------------------------------------------------------------------------
## Produce the Violin and Feature Plots

# Iterate throught the markers dataframe and produce the violin and feature plots for each celltype in the excel file
for (celltype in 1:length(markers)){
  
  #split the markers in individual character strings in vector
  marker <- strsplit(as.character(markers[celltype,2]), ", ") %>%
    unlist()
  
  #Produce the Violin and Feature plots for every cell type
  violinplot <- VlnPlot(seurat.obj, features = convert_markers(marker, species, prefix), pt.size = 0)
  featplot <- FeaturePlot(seurat.obj, features = convert_markers(marker, species, prefix))
  
  #Save the plots
  png(filename=paste0(plot.dir, species, "/3.CellTypeIdentification/", as.character(markers[celltype,1]), "_Violinplot.png"), width = 1000, height = 700)
  print(violinplot)
  dev.off()
  
  png(filename=paste0(plot.dir, species, "/3.CellTypeIdentification/", as.character(markers[celltype,1]), "_Featureplot.png"), width = 1000, height = 700)
  print(featplot)
  dev.off()
}


