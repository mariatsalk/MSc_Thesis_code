"
## snRNA analysis - GO enrichment analysis
## Author: Maria Tsalkitzidou
## Date: 14/03/2023
## Updated: 28/05/2023

I followed the tutorial from Harvard Chan Bioinformatics Core (HBC) learning hub:
https://hbctraining.github.io/DGE_workshop/lessons/09_functional_analysis.html 
"
###############################################################################################
#### Load the necessary packages and user defined variables ####

rm(list = ls()) #Remove (possible preloaded) objects from the environment to avoid conflicts
setwd("/Users/Maria/Dropbox (DNPL)/snRNA_Maria_Final/glial_snRNAseq_analysis/4.Scripts") ## Set the working directory

## Load the directories and necessary packages
source("Directories_Packages.R")
source("snRNA_Functions.R")


#----------------------------------------------------------------------------------------------
## User defined variables

obj.dir = "" #The object to be loaded in the script
species = "" #rat or human (all lowercase!!)
prefix = NULL #prefix in front of the marker (in this case: "premRNAGRCH38-", "premRNARnor6--" or NULL if the prefix is removed)
celltype = "" #cell type for which we want to run GSEA e.g. Astrocyte




###############################################################################################
#### Load the data ####
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))




###############################################################################################
#### GSEA per brain area ####
for (area in levels(as.factor(seurat.obj$location))){
  
  fgseaRes.location <- fgsea_analysis(seurat.object=seurat.obj, group_type = "location", group_name = area, species_name = species)
  
  
  # Tidy up the data
  fgseaResTidy.location <- fgseaRes.location %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  
  fgseaResTidy.location %>% 
    dplyr::select(-leadingEdge, -ES) %>% 
    arrange(padj) %>% 
    head()
  
  
  # Get only the Biological Process GO terms
  BP_fgseaResTidy.location <- fgseaResTidy.location[grep("GOBP_", fgseaResTidy.location$pathway),]
  
  GSEA.plot.location <- ggplot(BP_fgseaResTidy.location %>% dplyr::filter(padj < 0.05) %>% head(n= 10), aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj)) +
    coord_flip() +
    labs(x="GO term", y="Normalized Enrichment Score",
         title=paste0(str_to_title(species), " ", celltype, " ", area, " GO NES from GSEA")) + 
    theme_minimal() +
    theme(axis.text = element_text(size=12))
  
  
  png(filename=paste0(plot.dir, species, "/4.DE/", celltype, "_", area, "_GSEA_BP_Barplot.png"), width = 1000, height = 400)
  print(GSEA.plot.location)
  dev.off()
  
}



###############################################################################################
#### per cluster ####
for (cluster in levels(seurat.obj$seurat_clusters)){
  fgseaRes.cluster <- fgsea_analysis(seurat.object=seurat.obj, group_type = "seurat_clusters", group_name = as.character(cluster), species_name = species)
  
  # Tidy up the data
  fgseaResTidy.cluster <- fgseaRes.cluster %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  
  fgseaResTidy.cluster %>% 
    dplyr::select(-leadingEdge, -ES) %>% 
    arrange(padj) %>% 
    head()
  
  
  # Get only the Biological Process GO terms
  BP_fgseaResTidy.cluster <- fgseaResTidy.cluster[grep("GOBP_", fgseaResTidy.cluster$pathway),]
  
  GSEA.plot.cluster <- ggplot(BP_fgseaResTidy.cluster %>% dplyr::filter(padj < 0.05) %>% head(n= 10), aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj)) +
    coord_flip() +
    labs(x="GO term", y="Normalized Enrichment Score",
         title=paste0(str_to_title(species), " ", celltype, " Cluster ", cluster, " GO NES from GSEA")) +
    theme_minimal() +
    theme(axis.text = element_text(size=12))
  
  png(filename=paste0(plot.dir, species, "/4.DE/", celltype, "_Cluster", cluster, "_GSEA_BP_Barplot.png"), width = 1000, height = 400)
  print(GSEA.plot.cluster)
  dev.off()
  
}
