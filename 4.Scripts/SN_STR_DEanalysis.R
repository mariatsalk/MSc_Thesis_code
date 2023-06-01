## snRNA analysis - DE analysis for each cluster
## Author: Maria Tsalkitzidou
## Date: 13/01/2023
## Updated: 31/01/2023

setwd("/Users/Maria/Desktop/scRNAseq_analysis/4.Scripts")
#----------------------------------------------------------------------------------------------
## Load the necessary packages
source("Directories_Packages.R")
source("snRNA_Functions.R")
source("Markers.R")

## Colour palettes
human.astrocyte.palette <- c("Blue3", "deepskyblue", "deepskyblue3", "cadetblue3")
human.opc.palette <- c("#FDB863","#C9D12E", "#DC8315", "yellow4")


#-----------------------------------------------------------------------------------------------

source("AstrocyteSubtypesMarkers.R")
prefix="premRNARnor6--"
org = "Human"

obj.dir = "GPC_Human_obj.rds" #The object to be loaded in the script
dataset = "OPC" #Name of the dataset to be loaded in the script (in this case: Human, Rat, Merged)
species = "human"
chosen.palette <- human.opc.palette

#-----------------------------------------------------------------------------------------------
## Load the data
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))
DimPlot(seurat.obj)


#-----------------------------------------------------------------------------------------------
## Re-cluster the object

## Repeat the necessary steps with harmony embeddings
seurat.obj <- seurat.obj %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = counts@var.genes, npcs = 50, verbose = FALSE) #%>%
  #JackStraw(num.replicate = 100) %>%
  #ScoreJackStraw(dims=1:20)

JackStrawPlot(seurat.obj)
ElbowPlot(seurat.obj)

seurat.obj <- seurat.obj %>%
  RunHarmony(group.by.vars = "orig.ident") %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.05) %>% 
  identity()

DimPlot(s.obj)

s.obj <- readRDS(file = paste0(seurat.dir, "Human_OPC_subclusters_10PCs _res005_VarFeat_Scaled_HARMONY_obj.rds"))
  
  


source("Silhouette scores.R")
compute_silhouette_scores(s.obj, res.from = 0.1, res.to = 1.5, by= 0.01, plot=F)

#Produce the UMAPs per different 
dimplot.general <- DimPlot(seurat.obj, cols = chosen.palette)
dimplot.location <- DimPlot(seurat.obj, group.by = 'location')
dimplot.age <- DimPlot(seurat.obj, group.by = 'Age')

## Save the plots
png(filename = paste0(plot.dir, dataset, "_subclusters_10PCs_res005_UMAP_perCluster.png"), width = 900, height = 700)
print(dimplot.general)
dev.off()

png(filename = paste0(plot.dir, dataset, "_subclusters_10PCs_res005_UMAP_perBrainArea.png"), width = 900, height = 700)
print(dimplot.location)
dev.off()

png(filename = paste0(plot.dir, dataset, "_subclusters_10PCs_res005_UMAP_perAge.png"), width = 900, height = 700)
print(dimplot.age)
dev.off()

## Save the seurat object
saveRDS(seurat.obj, file = paste0(seurat.dir, dataset, "_subclusters_20PCs _res005_VarFeat_Scaled_HARMONY_obj.rds"))

seurat.obj <- readRDS(file = paste0(seurat.dir, org, "_", dataset, "_subclusters_10PCs_res005_VarFeat_Scaled_HARMONY_obj.rds"))

seurat.obj <- readRDS(file = paste0(seurat.dir, "Human_OPC_subclusters_10PCs_res005_VarFeat_Scaled_HARMONY_obj.rds"))
DimPlot(seurat.obj)
## Extract some basic information

nCell <- table(s.obj$seurat_clusters) %>%
  as.data.frame()

write.xlsx(nCell, file = paste0(plot.dir, qc.dir, species, "_", dataset, "_perCluster_nCells.xlsx"))
# Number of cells per brain area and age
nCell.location.age <- table(s.obj$location, s.obj$Age) %>%
  as.data.frame()


write.xlsx(nCell.location.age, file = paste0(plot.dir, qc.dir, species, "_", dataset, "_perBrainArea_perAge_nCells.xlsx"))


opc.obj <- readRDS(file = paste0(seurat.dir, "OPC", "_subclusters_VarFeat_Scaled_HARMONY_obj.rds"))

#-----------------------------------------------------------------------------------------------
## DE per cluster

#Create an excel workbook to store the names of the DE genes from each cluster
clusters.wb = createWorkbook()

#Iterate through the clusters to calculate the DE genes per cluster (against the rest of the clusters) and save the table to a sheet in the excel workbook
for (cluster in levels(seurat.obj$seurat_clusters)){
  
  #Find the markers
  cluster.genelist<- FindMarkers(seurat.obj, ident.1 = cluster, ident.2 = NULL)
  
  #Arrange the genes in a descending order based to avg_log2FC and keep only the genes with avg_log2FC > 1.0
  desc.cluster.genelist <- cluster.genelist %>%
    arrange(desc(avg_log2FC)) %>%
    filter(avg_log2FC > 1.0)
  
  #Get the names of the genes without the prefix
  desc.cluster.genelist$genename <- rownames(desc.cluster.genelist)
  desc.cluster.genelist$genename <- gsub(prefix, "", rownames(desc.cluster.genelist))
  
  #Create a separate variable for each genelist in case we want to investigate further
  assign(paste0("genelist.cluster", cluster), desc.cluster.genelist)
  
  #Create a worksheet in the excel workbook for this cluster and its DE genes
  seurat.cluster = addWorksheet(clusters.wb, paste0("Cluster ", cluster))
  
  #Write the data in the worksheet
  writeData(clusters.wb, x= desc.cluster.genelist, sheet=seurat.cluster, rowNames=F)
  
}


#Save the workbook
saveWorkbook(clusters.wb, file = paste0(plot.dir, dataset, "_10PCs_res005_DifferentiatedMarkers_perCluster.xlsx"), overwrite = T)


## DE per brain area

seurat.obj <- s.obj
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
    filter(avg_log2FC > 1.0)
  
  #Get the names of the genes without the prefix
  desc.region.genelist$genename <- rownames(desc.region.genelist)
  desc.region.genelist$genename <- gsub(prefix, "", rownames(desc.region.genelist))
  
  #Create a separate variable for each genelist in case we want to investigate further
  assign(paste0("genelist.region", region), desc.region.genelist)
  
  #Create a worksheet in the excel workbook for this cluster and its DE genes
  seurat.region = addWorksheet(location.wb, paste0("Region ", region))
  
  #Write the data in the worksheet
  writeData(location.wb, x= desc.region.genelist, sheet=seurat.region, rowNames=T)
  
}


#Save the workbook
saveWorkbook(location.wb, file = paste0(plot.dir, dataset, "_DifferentiatedMarkers_perBrainRegion.xlsx"), overwrite = T)


####################3

head(AverageExpression(seurat.obj, group.by = "location"))

#######################


temp.seurat.obj <- seurat.obj

Idents(temp.seurat.obj) <- "location"

genelist.location<- FindMarkers(temp.seurat.obj, ident.1 = "Nigra")

desc.genelist.location <- genelist.location %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 1.0)


desc.genelist.location$genename <- rownames(desc.genelist.location)
desc.genelist.location$genename <- gsub(prefix, "", rownames(desc.genelist.location))

write.xlsx(desc.genelist.location, file = paste0(plot.dir, dataset, "_DifferentiatedMarkers_perBrainArea.xlsx"), rowNames=T)




#-----------------------------------------------------------------------------------------------
## Check the expression of specific markers 

DimPlot(seurat.obj)
DimPlot(seurat.obj, group.by = "location")
DimPlot(seurat.obj, group.by = "Age")

# ASTROCYTE
#astro specific
VlnPlot(seurat.obj, features = c("GFAP", "AQP4", "S100B", "GJA1", "SLC1A2", "SLC1A3"), pt.size = 0, cols = chosen.palette)
VlnPlot(seurat.obj, features = c("AGT", "KCNC2", "GFAP", "SLC1A2", "GJA1", "SLC1A3"), pt.size = 0, group.by = "location")

#reactive astrocytes
VlnPlot(seurat.obj, features = c("CHI3L1", "C3", "CRYAB", "MAOB", "NFAT5", "HSPB1"), pt.size = 0, cols = chosen.palette)

#synapse function
VlnPlot(seurat.obj, features = c("AGT", "SLC7A10", "GLUL", "SLC6A11", "SLC6A1", "SLC6A9"), pt.size = 0, cols = chosen.palette)

#GPC

#opc and cop
VlnPlot(seurat.obj, features = str_to_upper(c("PDGFRA", "Cspg4", "Sox10", "Neu4", "Bmp4", "Tns3")), pt.size = 0, cols = chosen.palette)

#odc and cycling
VlnPlot(seurat.obj, features = str_to_upper(c("Opalin", "Mog", "Mobp", "Top2a", "Cdk1", "Cenpa")), pt.size = 0, cols = chosen.palette)




VlnPlot(seurat.obj, Astrocytes.M, pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(Astrocytes.M, prefix, org))

VlnPlot(seurat.obj, str_to_upper(astrocyte.development), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.development, prefix, org))

VlnPlot(seurat.obj, str_to_upper(astrocyte.iontransportation), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.iontransportation, prefix, org))

VlnPlot(seurat.obj, str_to_upper(astrocyte.immunefunction), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.immunefunction, prefix, org))

VlnPlot(seurat.obj, convert_markers(astrocyte.intermediateprogenitor, prefix, org), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.intermediateprogenitor, prefix, org))

VlnPlot(seurat.obj, str_to_upper(astrocyte.mature), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.mature, prefix, org))

VlnPlot(seurat.obj, convert_markers(astrocyte.phagocytosis, prefix, org), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.phagocytosis, prefix, org))

VlnPlot(seurat.obj, str_to_upper(astrocyte.progenitors), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.progenitors, prefix, org))

VlnPlot(seurat.obj, str_to_upper(astrocyte.synapsefunction), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.synapsefunction, prefix, org))

VlnPlot(seurat.obj, str_to_upper(astrocyte.synaptogenesis), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.synaptogenesis, prefix, org))

VlnPlot(seurat.obj, convert_markers(astrocyte.neuraltissuedevelopment, prefix, org), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.neuraltissuedevelopment, prefix, org))



VlnPlot(seurat.obj, str_to_upper(c("Slc1a3", "Slc1a2", "Crym", "sparc", "Gfap", "Slc6a1", "Gjb6", "Olr1", "Gins3")), group.by = "location", pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(c("Slc1a3", "Slc1a2", "Crym", "sparc", "Gfap", "Gjb6"), prefix, org))

VlnPlot(seurat.obj, str_to_upper(c("Agt", "Kcnc2", "Gfap", "Slc1a3", "Slc1a2", "Gja1")), group.by = "location", pt.size = 0)


VlnPlot(seurat.obj, convert_markers(astrocyte.synapsecommunication, prefix, org), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.synapsecommunication, prefix, org))

VlnPlot(seurat.obj, convert_markers(astrocyte.synapsesynthesis, prefix, org), pt.size = 0)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.synaptogenesis, prefix, org))






VlnPlot(seurat.obj, convert_markers(), prefix, org), pt.size = 0, cols = chosen.palette)
FeaturePlot(seurat.obj, features = convert_markers(astrocyte.synaptogenesis, prefix, org))