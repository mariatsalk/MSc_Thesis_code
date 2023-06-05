"
## snRNA analysis - Extract information for the glial subpopulations and get an excel file with the differentially exressed (DE) genes for every cluster and location
## Author: Maria Tsalkitzidou
## Date: 27/04/2023
## Updated: 28/05/2023
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
obj.dir = "" #The object to be loaded in the script
celltype = "" #cell type for which we want to run this script e.g. Astrocyte
species = "" #rat or human (all lowercase!!)




############################################################################################
## Load the data
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))
DimPlot(seurat.obj) #checkpoint



############################################################################################
#### QC check and information extraction ####

# QC data for the cluster
cluster.qc <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Save the QC plot
png(filename = paste0(plot.dir, species, "/4.DE/QC/", str_to_title(species), "_", celltype, "_QCplot.png"))
print(cluster.qc)
dev.off()

#-----------------------------------------------------------------------------------------------
# Get the number of cells per cluster per brain region
nCells.df <- as.data.frame(table(Idents(seurat.obj), seurat.obj$location))

#calculate the percentage of cells in each cluster
nCells.df$percentage <- round(100*(nCells.df$Freq/sum(nCells.df$Freq)),2)

#Calculate the total number of cells per cluster 
CellNumberperCluster <- nCells.df %>% group_by(Var1) %>% summarize(total = sum(Freq))


# Extract the table in an excel file
write.xlsx(nCells.df, file = paste0(plot.dir, species, "/4.DE/QC/", str_to_title(species), celltype, "_nCells_perCluster_perBrainRegion.xlsx"), rowNames=F)



############################################################################################
#### DE per cluster ####

#Create an excel workbook to store the names of the DE genes from each cluster
clusters.wb = createWorkbook()

#Iterate through the clusters to calculate the DE genes per cluster (against the rest of the clusters) and save the table to a sheet in the excel workbook
for (cluster in levels(seurat.obj$seurat_clusters)){
  
  #Find the markers
  cluster.genelist<- FindMarkers(seurat.obj, ident.1 = cluster, ident.2 = NULL)
  
  #Arrange the genes in a descending order based to avg_log2FC and keep only the genes with avg_log2FC > 1.0 and p_val_adj < 0.05
  desc.cluster.genelist <- cluster.genelist %>%
    arrange(desc(avg_log2FC)) %>%
    dpylr::filter(avg_log2FC > 1.0 & p_val_adj < 0.05)
  
  
  #Create a separate variable for each genelist in case we want to investigate further
  assign(paste0("genelist.cluster", cluster), desc.cluster.genelist)
  
  #Create a worksheet in the excel workbook for this cluster and its DE genes
  seurat.cluster = addWorksheet(clusters.wb, paste0("Cluster ", cluster))
  
  #Write the data in the worksheet
  writeData(clusters.wb, x= desc.cluster.genelist, sheet=seurat.cluster, rowNames=T)
  
}


#Save the workbook
saveWorkbook(clusters.wb, file = paste0(plot.dir, species, "/4.DE/", str_to_title(species), celltype, "_DifferentiatedMarkers_perCluster.xlsx"), overwrite = T)

############################################################################################
#### DE per brain area ####

#Create an excel workbook to store the names of the DE genes from each cluster
location.wb = createWorkbook()


temp.seurat.obj <- seurat.obj

Idents(temp.seurat.obj) <- "location"

#Iterate through the clusters to calculate the DE genes per cluster (against the rest of the clusters) and save the table to a sheet in the excel workbook
for (region in levels(Idents(temp.seurat.obj))){
  
  
  #Find the markers
  region.genelist<- FindMarkers(temp.seurat.obj, ident.1 = region, ident.2 = NULL)
  
  #Arrange the genes in a descending order based to avg_log2FC and keep only the genes with avg_log2FC > 1.0
  desc.region.genelist <- region.genelist %>%
    arrange(desc(avg_log2FC)) %>%
    dpylr::filter(avg_log2FC > 1.0 & p_val_adj < 0.05)
  

  
  #Create a separate variable for each genelist in case we want to investigate further
  assign(paste0("genelist.region", region), desc.region.genelist)
  
  #Create a worksheet in the excel workbook for this cluster and its DE genes
  seurat.region = addWorksheet(location.wb, paste0("Region ", region))
  
  #Write the data in the worksheet
  writeData(location.wb, x= desc.region.genelist, sheet=seurat.region, rowNames=T)
  
}


#Save the workbook
saveWorkbook(location.wb, file = paste0(plot.dir, species, "/4.DE/", str_to_title(species), celltype, "_DifferentiatedMarkers_perBrainRegion.xlsx"), overwrite = T)




