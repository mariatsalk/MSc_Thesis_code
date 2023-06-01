## scRNA analysis - R packages and Directories Script
## Author: Maria Tsalkitzidou
## Date: 19/12/2022

##This script contains all the directories used in this analysis and also loads all the packages needed.

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




#-----------------------------------------------------------------------------------------------
## Define directories
data.dire = "../1.RawData/"
plot.dir = "../3.Plots_Stats/"
seurat.dir = "../2.SeuratObjects/"
qc.dir = "1.QC/"
norm.dir = "2.NormalizationScaling/"

if(!dir.exists(plot.dir)) dir.create(plot.dir)
if(!dir.exists(qc.dir)) dir.create(qc.dir)
if(!dir.exists(norm.dir)) dir.create(norm.dir)
#-----------------------------------------------------------------------------------------------
## Color palettes
cbPalette5_human <- c("deepskyblue4", "chocolate3", "cadetblue3", "slateblue3", "#FDB863")
cbPalette5a_human <- c("deepskyblue4", "#FDB863", "chocolate3", "cadetblue3", "slateblue3")
cbPalette10_rat <- c("#E69F00", "#56B4E9", "#117733", "#F0E442", "#0072B2", "#D55E00", "#CC6677", "#AA4499", "#44AA99", "#999933")