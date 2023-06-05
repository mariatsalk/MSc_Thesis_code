"
## scRNA analysis - R packages and Directories Script
## Author: Maria Tsalkitzidou
## Date: 19/12/2022

##This script contains all the directories used in this analysis and also loads all the packages needed.
"
##########################################################################################
## Load the necessary packages that are used in all scripts
library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(rliger)
library(ggplot2)
library(stringr)
library(openxlsx)
library(tidyverse)
library(RColorBrewer)



#-----------------------------------------------------------------------------------------------
## Define directories
data.dire = "../1.RawData/"
plot.dir = "../3.Plots_Stats/"
seurat.dir = "../2.SeuratObjects/"
qc.dir = "1.QC/"
norm.dir = "2.NormalizationScaling/"

if(!dir.exists(plot.dir)) dir.create(plot.dir)

#-----------------------------------------------------------------------------------------------
## Color palettes
celltypes.human <- c("deepskyblue4", "chocolate3", "cadetblue3", "slateblue3", "#FDB863")
human.astrocyte.palette <- c("Blue3", "deepskyblue", "deepskyblue3", "cadetblue3")
human.opc.palette <- c("#FDB863","#C9D12E", "#DC8315", "yellow4")

rat.microglia.palette = c("#cb4335", "#ec7063", "#FF5733", "#C70039", "#d35400", "#7b241c")
rat.astrocyte.palette = c("Blue3", "deepskyblue3", "cadetblue3", "deepskyblue", "Blue")
rat.oligo.palette = c("#72d693", "#117733", "green2")
rat.opc.palette = c("#E69F00", "darkorange1", "coral2", "yellow3", "#CB910F")