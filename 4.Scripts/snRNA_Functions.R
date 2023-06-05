"
## snRNA analysis - Functions
## Author: Maria Tsalkitzidou
## Created: 08/09/2022
## Updated: 24/10/2022
"

####Load the necessary packages ####
library(Seurat)
library(dplyr)
library(patchwork)
library(readxl)




################################################################################################################
#### Function for automated species assignment ####

# The function takes as input a seurat object, one or two species as string, one or two regex patterns (corresponding to the species) that will be used to separate the species, an upper threshold for the PercentageFeatureSet() and one lower threshold for PercentageFeatureSet()
species_assignment <- function(seurat.obj = data.obj,
                               species.1 = "Species1",
                               species.2 = NULL,
                               pattern.species.1 = "pattern1",
                               pattern.species.2 = NULL,
                               upper.thresh = 80,
                               lower.thresh = 10){
  
  # Check if one or two species were provided
  if(!is.null(species.2)){
    #Check that the user also provided a pattern, if not print error message
    if(is.null(pattern.species.2)){
      stop("Please provide a pattern for all the species")
    }
    
    #Add the species information as metadata to the seurat object
    seurat.obj$species.1 <- PercentageFeatureSet(seurat.obj, pattern = pattern.species.1)
    seurat.obj$species.2 <- PercentageFeatureSet(seurat.obj, pattern = pattern.species.2)
    seurat.obj$Species <- NA
    
    for (i in 1:length(seurat.obj$Species)){
      if (seurat.obj$species.1[i] > upper.thresh && seurat.obj$species.2[i] < lower.thresh){
        seurat.obj$Species[i] <- species.1
      } else if (seurat.obj$species.2[i] > upper.thresh && seurat.obj$species.1[i] < lower.thresh){
        seurat.obj$Species[i] <- species.2
      }
      else{
        seurat.obj$Species[i] <- "Undefined"
      }
    }
  } else {
    #Add the species information as metadata to the seurat object
    seurat.obj$species.1 <- PercentageFeatureSet(seurat.obj, pattern = pattern.species.1)
    seurat.obj$Species <- NA
    
    for (i in 1:length(seurat.obj$Species)){
      if (seurat.obj$species.1[i] > upper.thresh){
        seurat.obj$Species[i] <- species.1
      } else{
        seurat.obj$Species[i] <- "Undefined"
      }
    }
    
  }
  return(seurat.obj)  
}




################################################################################################################
#### Function for converting the case of markers based on species #### 
# This function takes as input a vector or list of markers and the organism we want to convert the marker for. The output is a vector with the converted markers
convert_markers <- function(markers, organism, organism.prefix=NULL){
  library(stringr)
  converted.markers <- c()
  
  # Change the user input to lowercase to avoid case-sensitivity
  organism <- str_to_lower(organism) 
  
  #If the user wants to convert to human markers:
  if (organism == 'human'){
    #go through the provided marker vector or list
    for (marker in markers){
      # Change the marker format to the appropriate one (for the human model all the letters should be in uppercase)
      converted.markers <- append(converted.markers, str_to_upper(marker))
    }
  } else if (organism == 'rat'){
    for (marker in markers){
      #for the rat model the first letter should be capital and everything else in lowercase
      converted.markers <- append(converted.markers, str_to_title(marker))
    }
  } else{
    stop('No available conversion for this organism')
  }
  
  # Add the prefix of the marker if the user provides one
  if (!is.na(organism.prefix)){
    converted.markers <- paste0(prefix, converted.markers)
  }
  
  
  # Return the converted marker vector for our dataset
  return(converted.markers)
}





################################################################################################################
#### Function for GSEA analysis with fgsea ####
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
  
  #Rename the species with the proper name
  if (str_to_lower(species_name) == 'human') species_name == "Homo sapiens"
  if (str_to_lower(species_name) == 'rat') species_name == "Rattus norvegicus"
  
  print(species_name)
  
  #Run the wilcoxon test
  seurat.genes <- wilcoxauc(seurat.object, group_type)
  
  #Remove the prefix in the gene symbol
  seurat.genes$feature <- gsub(prefix, "", seurat.genes$feature)
  
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
    dplyr::filter(logFC > 1.0) %>%
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

