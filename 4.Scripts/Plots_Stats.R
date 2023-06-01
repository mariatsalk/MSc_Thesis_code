"
## snRNA analysis - Dataset Statistics
## Author: Maria Tsalkitzidou
## Date: 31/10/2022
## Updated: 13/01/2022

Description:
  The script takes as input a seurat object and produces several stacked bar plots and tables in csv format. The stacked bar plots produced are:
  1) Types of cells per species (with %)
  2) Types of cells per samples (with %)
  3) Number of cells per species (with %)
  4) Number of cells per samples (with %)
  5) Types of cells per species & per samples (with %)
The tables contain information about the number of cells (with %) per dataset (Human, Rat, Merged dataset)


Procedure:



Limitations:
  1) Doesn't take input from the terminal
  2) Isn't 100% generic

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
dataset = "vMB_hVM_merged.seurat.all.20210920" #The Seurat object will start with this phrase/word

obj.dir = "vMB_hVM_merged.seurat.all.20210920_raw_merged_object.rds" #Seurat object to load

species = "human" #rat or human (all lowercase!!)

colour.palette <- c()


###############################################################################################
#### Load the Seurat object ####

seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))



###############################################################################################
#### CELL TYPE PER SAMPLE STACKED BAR PLOT ####


#Create a table that consists of the clusters, the sample corresponding to each cluster and the number of cells in each cluster
CellTypePerSample.table <- table(Idents(seurat.obj), seurat.obj$orig.ident, seurat.obj$location)
CellTypePerSample.table <- as.data.frame(CellTypePerSample.table)


#calculate the percentage of cells in each cluster
CellTypePerSample.table$percentage <- round(100*(CellTypePerSample.table$Freq/sum(CellTypePerSample.table$Freq)),2)

CellTypePerSample.table <- CellTypePerSample.table %>%
  subset(Freq!=0) %>%
  group_by(Var3)
  #mutate(Var1 = factor(Var1, levels = celltype.order))
  


#Calculate the total number of cells per sample 
CellNumberperSample <- CellTypePerSample.table %>% group_by(Var3, Var2) %>% summarize(total = sum(Freq))




#Plot the stacked bar graph % cell type per sample
CellTypePerSample.StackedBarPlot <- ggplot(CellTypePerSample.table, aes(x = Var2, y = percentage)) +
  facet_grid(~ Var3, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 15) +
  geom_col(aes(fill = Var1), position="fill") +
  geom_text(data = CellNumberperSample, aes(x=Var2, y=1, label= total, fill=NULL), position = position_fill(vjust = 1.02) , size=4) +
  scale_fill_manual(values = colour.palette) +
  ggtitle("Cell Type per Samples") +
  xlab("Samples") +
  ylab("Abundance (%)") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        plot.margin = margin(10, 10, 10, 110),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 30),
        axis.text.x = element_blank(), #element_text(angle = 45, size = 11, hjust = 1),
        strip.text = element_text(size = 20, face = "bold"),
        strip.background = element_blank())





#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_POSTER_LEGEND_ON_BOTTOM_CellTypePerSample_StackedBarPlot.png"), width = 1700, height = 900)
print(CellTypePerSample.StackedBarPlot)
dev.off()








#---------------------------------------------------------------------------------------------
#### CELL TYPE PER BRAIN AREA ####

#Create a table that consists of the clusters, the species information for each cluster and the number of cells in each cluster
CellTypePerBrainArea.table <- table(Idents(seurat.obj), seurat.obj$location)
CellTypePerBrainArea.table <- as.data.frame(CellTypePerBrainArea.table)

#calculate the percentage of cells in each cluster
CellTypePerBrainArea.table$percentage <- round(100*(CellTypePerBrainArea.table$Freq/sum(CellTypePerBrainArea.table$Freq)),2)

#Transform the name of the clusters in characters
CellTypePerBrainArea.table$Var1 <- as.character(CellTypePerBrainArea.table$Var1)

#Calculate the total number of cells per sample 
CellNumberperBrainArea <- CellTypePerBrainArea.table %>% group_by(Var2) %>% summarize(total = sum(Freq))



#Plot the stacked bar graph % cell type per species
CellTypePerBrainArea.StackedBarPlot <- ggplot(CellTypePerBrainArea.table, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_text(data = CellNumberperBrainArea, aes(x=Var2, y=1, label= total, fill=NULL), position = position_fill(vjust = 1.02) , size=4) +
  scale_fill_manual(values = colour.palette) +
  ggtitle("Cell Type per Brain Area (%)") +
  xlab("Samples") +
  ylab("Proportion") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1.5))

#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_CellTypePerBrainArea_StackedBarPlot.png"), width = 800, height = 800)
print(CellTypePerBrainArea.StackedBarPlot)
dev.off()





#---------------------------------------------------------------------------------------------
#### CELL TYPE PER AGE ####

#Create a table that consists of the clusters, the species information for each cluster and the number of cells in each cluster
CellTypePerAge.table <- table(Idents(seurat.obj), seurat.obj$Age, seurat.obj$location)
CellTypePerAge.table <- as.data.frame(CellTypePerAge.table)

#calculate the percentage of cells in each cluster
CellTypePerAge.table$percentage <- round(100*(CellTypePerAge.table$Freq/sum(CellTypePerAge.table$Freq)),2)

CellTypePerAge.table$Var2 <- factor(CellTypePerAge.table$Var2, levels = c("3m", "6m", "9m", "12m"))


CellTypePerAge.table <- CellTypePerAge.table %>%
  subset(Freq!=0) %>%
  group_by(Var3)

#Calculate the total number of cells per sample 
CellNumberperAge <- CellTypePerAge.table %>% group_by(Var3, Var2) %>% summarize(total = sum(Freq))

#Plot the stacked bar graph % cell type per species
CellTypePerAge.StackedBarPlot <- ggplot(CellTypePerAge.table, aes(x = Var2, y = Freq, fill =Var1)) +
  facet_grid(~ Var3, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_text(data = CellNumberperAge, aes(x=Var2, y=1, label= total, fill=NULL), position = position_fill(vjust = 1.02) , size=4) +
  scale_fill_manual(values = colour.palette) +
  ggtitle("Cell Type per Age") +
  xlab("Samples") +
  ylab("Proportion") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1.5),
        strip.text = element_text(size = 20, face = "bold"),
        strip.background = element_blank())


#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_POSTER_CellTypePerAge_StackedBarPlot.png"), width = 1200, height = 700)
print(CellTypePerAge.StackedBarPlot)
dev.off()


#-----------------------------------------------------------------------------------------------
#### NUMBER OF CELLS PER SAMPLE ####
#Create a table that consists of the species and the number of cells per sample
CellsPerSample.table <- table(seurat.obj$orig.ident)
CellsPerSample.table <- as.data.frame(CellsPerSample.table)

#calculate the percentage of cells in each sample
CellsPerSample.table$percentage <- round(100*(CellsPerSample.table$Freq/sum(CellsPerSample.table$Freq)),2)

#Transform the name of the sample in characters
CellsPerSample.table$Var1 <- as.character(CellsPerSample.table$Var1)

#geom_text(aes(label = sprintf("%0.2f", percentage)), position = position_fill(vjust=0.5), size = 4) +
#Plot the stacked bar graph % cells per sample
CellsPerSample.StackedBarPlot <- ggplot(CellsPerSample.table, aes(x = " ", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_bar(stat = "identity", color = "white") +
  ggtitle("Number of Cells per Samples (%)") +
  xlab("Samples") +
  ylab("Proportion") +
  theme(legend.text = element_text(size=5),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1.5))

#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_CellsPerSample_StackedBarPlot.png"), width = 1000, height = 900)
print(CellsPerSample.StackedBarPlot)
dev.off()



#-----------------------------------------------------------------------------------------------
#### NUMBER OF CELLS PER BRAIN AREA ####

#Create a table that consists of the species and the number of cells per species
CellsPerBrainArea.table <- table(seurat.obj$location)
CellsPerBrainArea.table <- as.data.frame(CellsPerBrainArea.table)

#calculate the percentage of cells in each species
CellsPerBrainArea.table$percentage <- round(100*(CellsPerBrainArea.table$Freq/sum(CellsPerBrainArea.table$Freq)),2)

#Transform the name of the species in characters
CellsPerBrainArea.table$Var1 <- as.character(CellsPerBrainArea.table$Var1)


#Plot the stacked bar graph % cells per species
CellsPerBrainArea.StackedBarPlot <- ggplot(CellsPerBrainArea.table, aes(x = " ", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_text(aes(label = sprintf("%0.2f", percentage)), position = position_fill(vjust=0.5), size = 4) +
  ggtitle("Number of Cells per Brain Area (%)") +
  xlab("Species") +
  ylab("Proportion") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1.5))

#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_CellsPerBrainArea_StackedBarPlot.png"), width = 800, height = 700)
print(CellsPerBrainArea.StackedBarPlot)
dev.off()




#-----------------------------------------------------------------------------------------------
#### NUMBER OF CELLS PER AGE ####

#Create a table that consists of the species and the number of cells per species
CellsPerAge.table <- as.data.frame(table(seurat.obj$Age))


#calculate the percentage of cells in each species
CellsPerAge.table$percentage <- round(100*(CellsPerAge.table$Freq/sum(CellsPerAge.table$Freq)),2)

#Transform the name of the species in characters
CellsPerAge.table$Var1 <- as.character(CellsPerAge.table$Var1)


#Plot the stacked bar graph % cells per species
CellsPerAge.StackedBarPlot <- ggplot(CellsPerAge.table, aes(x = " ", y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_text(aes(label = sprintf("%0.2f", percentage)), position = position_fill(vjust=0.5), size = 4) +
  ggtitle("Number of Cells per Age (%)") +
  xlab("Species") +
  ylab("Proportion") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1.5))

#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_CellsPerBrainArea_StackedBarPlot.png"), width = 800, height = 700)
print(CellsPerBrainArea.StackedBarPlot)
dev.off()


# Pie chart
ggplot(CellsPerBrainArea.table, aes(x="", y=Freq, fill=Var1)) +
  theme_bw(base_size = 15) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())



###############################################################################################
#### Export information to EXCEL ####


## Produce the tables with the % and # of cell types for the whole dataset, the human dataset and the rat dataset
# Subset for the human dataset
# use regex to exclude rat clusters
Striatum.CellTypePerBrainArea.table <- subset(CellTypePerBrainArea.table, grepl('^[^Nigra].+', Var1) & Var2 == 'Striatum', select = c("Var1", "Freq", "percentage"))
# Re-calculate the percentage now that we have just the human dataset
Striatum.CellTypePerBrainArea.table$percentage <- round(100*(Striatum.CellTypePerBrainArea.table$Freq/sum(Striatum.CellTypePerBrainArea.table$Freq)),2)
# Rename columns
colnames(Striatum.CellTypePerBrainArea.table) <- c("Cell Type", "Number of Cells", "Percentage of Cells")



# Subset for the rat dataset
# use regex to exclude human clusters
Nigra.CellTypePerBrainArea.table <- subset(CellTypePerBrainArea.table, grepl('^[^Striatum].+', Var1) & Var2 == 'Nigra', select = c("Var1", "Freq", "percentage"))
# Re-calculate the percentage now that we have just the human dataset
Nigra.CellTypePerBrainArea.table$percentage <- round(100*(Nigra.CellTypePerBrainArea.table$Freq/sum(Nigra.CellTypePerBrainArea.table$Freq)),2)
# Rename columns
colnames(Nigra.CellTypePerBrainArea.table) <- c("Cell Type", "Number of Cells", "Percentage of Cells")



# Subset for the whole dataset
Merged.CellTypePerBrainArea.table <- subset(CellTypePerBrainArea.table, select = c("Var1", "Freq", "percentage"))
# Sum the # and % of identical cell types (but different species)
Merged.CellTypePerBrainArea.table <- Merged.CellTypePerBrainArea.table %>%
  group_by(Var1) %>%
  summarise_all(sum)
# Rename columns
colnames(Merged.CellTypePerBrainArea.table) <- c("Cell Type", "Number of Cells", "Percentage of Cells")



#-----------------------------------------------------------------------------------------------
## Write tables to an excel workbook


#Create an excel workbook to store the stats
stats.wb = createWorkbook()

#Create a worksheet in the excel workbook for the dataframes
striatum.cells = addWorksheet(stats.wb, "Striatum Cells per Cell Type")
nigra.cells = addWorksheet(stats.wb, "Nigra Cells per Cell Type")
all.cells = addWorksheet(stats.wb, "Cells per Cell Type")

#Write the data in the worksheet
writeData(stats.wb, x= Merged.CellTypePerBrainArea.table, sheet=all.cells, rowNames=F)
writeData(stats.wb, x= Striatum.CellTypePerBrainArea.table, sheet=striatum.cells, rowNames=F)
writeData(stats.wb, x= Nigra.CellTypePerBrainArea.table, sheet=nigra.cells, rowNames=F)
  


#Save the workbook
saveWorkbook(stats.wb, file = paste0(plot.dir, str_to_title(species), "_CellsPerCellType.xlsx"), overwrite = T)


