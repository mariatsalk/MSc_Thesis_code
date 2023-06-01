## snRNA analysis - Gene Set Enrichment Analysis for Human DA neurons
## Author: Maria Tsalkitzidou
## Date: 11/4/2023
## Updated: 


#----------------------------------------------------------------------------------------
## Define the neseccary directories, files and objects
setwd("/Users/Maria/Desktop/scRNAseq_analysis/4.Scripts/")
source("Directories_Packages.R")
source("snRNA_Functions.R")

## User defined variables
#obj.dir = "SubstantiaNigra&Striatum_RatCells_RemovedLowQualitySamples&HumanGenes_29112022.rds" #The object to be loaded in the script
obj.dir = "Human_SN_STR_v2_HARMONY_normalizedscaled_annotated_obj.rds" #The object to be loaded in the script
dataset = "Human_OPC" #Name of the dataset to be loaded in the script (in this case: Human, Rat, Merged)
prefix = "premGRCH38-"


#----------------------------------------------------------------------------------------

fgsea_analysis <- function(seurat.object, group_type, group_name, species_name){
  # Install presto which performs a fast Wilcoxon rank sum test
  #install_github('immunogenomics/presto')
  #BiocManager::install("fgsea")
  library(devtools)
  library(presto)
  library(msigdbr)
  library(fgsea)
  
  print("Processing...")
  
  print(seurat.object)
  print(group_type)
  print(group_name)
  
  #Run the wilcoxon test
  seurat.genes <- wilcoxauc(seurat.object, group_type)
  
  print("DE genes per group type found..")
  
  #Remove the prefix in the gene symbol
  #seurat.genes$feature <- gsub(prefix, "", seurat.genes$feature)
  head(seurat.genes)
  
  #Check out the available species and collections in the database and choose the wanted species and category
  msigdbr_species()
  print(msigdbr_collections(), n=30)
  m_df<- msigdbr(species = species_name, category = "C5")
  
  m_df <- subset(m_df, m_df$gs_subcat != "HPO")
  
  #tail(m_df)
  
  # Create a gene set where genes are grouped per GO function
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  
  
  # select only the feature and auc columns for fgsea, which statistics to use is an open question
  groupname.genes <- seurat.genes %>%
    dplyr::filter(group == group_name) %>%
    arrange(desc(logFC), desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  
  ranks.groupname <- deframe(groupname.genes)
  
  head(ranks.groupname)
  
  #Run GSEA 
  fgseaRes.groupname <- fgsea(fgsea_sets,
                              stats = ranks.groupname,
                              eps=0.0,
                              minSize = 10,
                              maxSize = length(ranks.groupname)-1)
  
  
  # Return the GSEA result
  return(fgseaRes.groupname) 
  
  
}


#----------------------------------------------------------------------------------------




## Read the object
s.obj <- readRDS(file = paste0(seurat.dir, obj.dir))




#----------------------------------------------------------------------------------------------
## DA neurons for all clusters
seurat.obj <- readRDS(file = paste0(seurat.dir, dataset, "_subclusters_20PCs_VarFeat_Scaled_HARMONY_obj.rds"))

seurat.obj <- s.obj
DimPlot(seurat.obj)
## NIGRA

fgseaRes.nigra <- fgsea_analysis(seurat.object=seurat.obj, group_type = "seurat_clusters", group_name = "0", species_name = "Homo sapiens")

# Tidy up the data
fgseaResTidy.nigra <- fgseaRes.nigra %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy.nigra %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()


#Plot a barplot with the normalized GO enrichment score (the first 20)
ggplot(fgseaResTidy.nigra %>% dplyr::filter(padj < 0.05) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title= paste0(dataset, " Cluster 0 GO NES from GSEA")) + 
  theme_minimal()


# Get only the Biological Process GO terms
BP_fgseaResTidy.nigra <- fgseaResTidy.nigra[grep("GOBP_", fgseaResTidy.nigra$pathway),]

ggplot(BP_fgseaResTidy.nigra %>% dplyr::filter(padj < 0.05) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title=paste0(dataset, " Cluster 1 GO NES from GSEA")) + 
  theme_minimal()




## STRIATUM

fgseaRes.striatum <- fgsea_analysis(seurat.object=seurat.obj, group_type = "location", group_name = "Striatum", species_name = "Homo sapiens")

# Tidy up the data
fgseaResTidy.striatum<- fgseaRes.striatum %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy.striatum %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()


#Plot a barplot with the normalized GO enrichment score (the first 20)
ggplot(fgseaResTidy.striatum %>% dplyr::filter(padj < 0.05) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title= paste0(dataset, " Striatum GO NES from GSEA")) + 
  theme_minimal()


# Get only the Biological Process GO terms
BP_fgseaResTidy.striatum <- fgseaResTidy.striatum[grep("GOBP_", fgseaResTidy.striatum$pathway),]

ggplot(BP_fgseaResTidy.striatum %>% dplyr::filter(padj < 0.05) %>% head(n= 10), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title=paste0(dataset, " Striatum GO NES from GSEA")) + 
  theme_minimal()


#----------------------------------------------------------------------------------------------
## DA neurons for TH+ clusters

VlnPlot(seurat.obj, features = "TH")
FeaturePlot(seurat.obj, features = "TH")

th.obj <- subset(seurat.obj, subset = TH > 0)




## NIGRA

fgseaRes.nigraTH <- fgsea_analysis(seurat.object=th.obj, group_type = "location", group_name = "Nigra")

# Tidy up the data
fgseaResTidy.nigraTH <- fgseaRes.nigraTH %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy.nigraTH %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()


#Plot a barplot with the normalized GO enrichment score (the first 20)
ggplot(fgseaResTidy.nigraTH %>% dplyr::filter(padj < 0.05) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title= paste0(dataset, " Nigra TH+ GO NES from GSEA")) + 
  theme_minimal()


# Get only the Biological Process GO terms
BP_fgseaResTidy.nigraTH <- fgseaResTidy.nigraTH[grep("GOBP_", fgseaResTidy.nigraTH$pathway),]

ggplot(BP_fgseaResTidy.nigraTH %>% dplyr::filter(padj < 0.05) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title=paste0(dataset, " Nigra TH+ GO NES from GSEA")) + 
  theme_minimal()




## STRIATUM

fgseaRes.striatumTH <- fgsea_analysis(seurat.object=th.obj, group_type = "location", group_name = "Striatum")

# Tidy up the data
fgseaResTidy.striatumTH <- fgseaRes.striatumTH %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy.striatumTH %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()


#Plot a barplot with the normalized GO enrichment score (the first 20)
ggplot(fgseaResTidy.striatumTH %>% dplyr::filter(padj < 0.05) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title= paste0(dataset, " Striatum TH+ GO NES from GSEA")) + 
  theme_minimal()


# Get only the Biological Process GO terms
BP_fgseaResTidy.striatumTH <- fgseaResTidy.striatumTH[grep("GOBP_", fgseaResTidy.striatumTH$pathway),]

ggplot(BP_fgseaResTidy.striatumTH %>% dplyr::filter(padj < 0.05) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="GO term", y="Normalized Enrichment Score",
       title=paste0(dataset, " Striatum TH+ GO NES from GSEA")) + 
  theme_minimal()


