 "
## snRNA analysis - Identify cell types by signature markers
## Author: Maria Tsalkitzidou
## Date: 11/10/2022
## Updated: 15/02/2023
"
############################################################################################### 
#### Load the necessary packages and user defined variables ####

rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts") ## Set the working directory
 
## Load the directories and necessary packages
source("Directories_Packages.R")
source("snRNA_functions.R")

#-------------------------------------------------------------------------------------------------------
## User defined variables
obj.dir = "" #The object to be loaded in the script
species = "" #Name of the dataset to be loaded in the script (in this case: Human, Rat, Merged)
prefix = NULL #prefix in front of the marker (in this case: "premRNAGRCH38-", "premRNARnor6--" or NULL if the prefix is removed)
markers.file = "" #path for excel file with markers to be tested. The first column should include the cell type and the second column the markers.


###############################################################################################
#### Load the data ####
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))

markers <- read_excel(markers.file, sheet = 1)

###############################################################################################
#### Produce the Violin and Feature Plots ####

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


