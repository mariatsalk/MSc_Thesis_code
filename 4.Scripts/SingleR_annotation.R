## Define the neseccary directories, files and objects
source("Directories_Packages.R")
library(SingleR)
library(celldex)
library(scRNAseq)
library(scater)
library(Seurat)

#-----
## Define the necessary directories, files and objects
## User defined variables
obj.dir = "Human_SN&STR_v2_HARMONY_normalized&scaled_obj.rds" #The object to be loaded in the script
ref.dir = "hvm_linnarsson.singler_annotation.20200403.rds"

brewerPallete2 <- c("deepskyblue4", "chocolate3", "cadetblue3", "slateblue3", "#FDB863")
#-----

## Load our dataset adnreference
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))
ref.dat <- readRDS(file = paste0(seurat.dir, ref.dir))
rownames(ref.dat)
levels(ref.dat$X2)


## Remove the uneccessary cell types from the reference to make the computations faster
unwanted.celltypes <- c("hEndo","hGaba","hMgl","hNbGaba" ,"hNbM","hNbML1","hNbML5","hOMTN", "hPeric","hProgBP", "hProgFPL", "hProgFPM", "hProgM", "hRN","hSert", "Unk")
#copy the reference dataset
ref <- ref.dat
#remove the unwanted cell types
for (ctype in unwanted.celltypes){
  ref <- ref[,ref$X2 != ctype]
}

#check that the unwated types are removed
ref$X2
levels(ref$X2 == "hEndo")

## Compute the log-expression values for the reference dataset for use in marker detection
ref.log <- logNormCounts(ref)

## Examine the distribution of labels in this reference
table(ref.dat.log$X2)
rownames(ref.dat.log)

# Run singleR
sn.str.singleRlabels1 <- SingleR(test = as.SingleCellExperiment(seurat.obj), ref = ref.log, labels = ref.log$X2)

# Add the singleR labels as metadata information to the seurat obj
seurat.obj$singleR_labels2 <- sn.str.singleRlabels1$labels

DimPlot(seurat.obj, group.by = "singleR_labels2")











ref <- ref.dat[,ref.dat$X2 != 'hGaba']
ref <- ref[,ref$X2 != 'hNb.+']
ref <- ref[,ref$X2 != 'hOMTN']

c("hDA0","hDA1","hDA2","hOPC","hNProg", "hRgl1","hRgl2a","hRgl2b","hRgl2c","hRgl3")
ref <- ref.dat[,ref.dat$X2 == 'hDA0']
ref <- ref[,ref$X2 == 'hDA1']
ref <- ref[,ref$X2 == 'hDA2']
ref <- ref[,ref$X2 == 'hOPC']
ref <- ref[,ref$X2 == 'hNProg']
ref <- ref[,ref$X2 == 'hRgl1']
ref <- ref[,ref$X2 == 'hRgl1']
ref$X2
levels(ref$X2)
