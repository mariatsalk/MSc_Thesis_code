## THESIS FINAL PLOTS FOR REPORT

#-------------------------------------------------------------------------------------------
## Define the neseccary directories, files and objects
source(".\\Directories&Packages.R")
source("Markers.R")
source("SignatureMarkers.R")
source("snRNA_Functions.R")

obj.dir = "SN&STR_Filtered&RemovedLowQualitySamples&HumanGenes_HARMONY_normalized&scaled_annotated_obj.rds" #The object to be loaded in the script
dataset = "Rat" #Name of the dataset to be loaded in the script (in this case: Human, Rat, Merged)
species_prefix ="premRNARnor6--"


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <- c("#0072B2", "#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#E69F00",  "#D55E00")
safePalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
cbPalette10 <- c("#E69F00", "#56B4E9", "#117733", "#F0E442", "#0072B2", "#D55E00", "#CC6677", "#AA4499", "#44AA99", "#999933")
cbPalette10A <- c("#56B4E9", "#F0E442", "#CC6677", "#AA4499", "#E69F00", "#D55E00", "#0072B2", "#117733", "#44AA99", "#999933")
cbPalette10B <- c("#56B4E9", "#F0E442", "#CC6677", "#AA4499", "#117733", "#E69F00", "#0072B2", "#D55E00", "#44AA99", "#999933")
cbPalette10C <- c("#56B4E9", "cadetblue2", "#CC6677", "#AA4499", "#117733", "#E69F00", "#0072B2", "red3", "#44AA99", "#F0E442")
cbPalette10D <- c("#117733", "#E69F00", "#0072B2", "red3", "#56B4E9", "cadetblue2", "#CC6677", "#AA4499", "#44AA99", "#F0E442")

microglia.palette = "RdYIGn"
astroglia.palette = "Blues" #c("Blue3", "deepskyblue3", "cadetblue3", "deepskyblue")
oligo.palette = "Oranges"
opc.palette = "Reds"
glia.palette = data.frame(microglia.palette, astroglia.palette, oligo.palette, opc.palette)

#-------------------------------------------------------------------------------------------
## Load the data
seurat.obj <- readRDS(file = paste0(seurat.dir, obj.dir))

#Change the order of the seurat clusters
celltype.order <- c("MSN-D1", "MSN-D2", "Interneurons", "Glutaminergic", "Oligodendrocyte", "OPC", "Astrocyte", "Immune cells", "Endothelia", "Pericyte")
celltype.order2 <- c("Oligodendrocyte", "OPC", "Astrocyte", "Immune cells", "MSN-D1", "MSN-D2", "Interneurons", "Glutaminergic", "Endothelia", "Pericyte")
names(celltype.order) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj, celltype.order)

Idents(seurat.obj) <- "seurat_clusters"

annotated.seurat.obj <- seurat.obj
annotated.seurat.obj <- RenameIdents(seurat.obj, celltype.order2)
annotated.seurat.obj$celltype <- Idents(seurat.obj)
annotated.seurat.obj$celltype <- factor(annotated.seurat.obj$celltype,
                                        levels = c("Oligodendrocyte", "OPC", "Astrocyte", "Immune cells", "MSN-D1", "MSN-D2", "Interneurons", "Glutaminergic", "Endothelia", "Pericyte"))

names(celltype.order2) <- levels(annotated.seurat.obj)
annotated.seurat.obj <- RenameIdents(annotated.seurat.obj, celltype.order2)


# Remove stupid prefix
counts <- GetAssayData(annotated.seurat.obj, slot = "counts")
rownames(counts) <- gsub(species_prefix, "", rownames(counts))
# Extract the meta data information and the active.ident information from the initial object and store it in a variable
new.meta.data <- annotated.seurat.obj@meta.data
new.active.ident <- annotated.seurat.obj@active.ident
s.obj <- CreateSeuratObject(counts, meta.data = new.meta.data, assay = "RNA")
Idents(s.obj) <- new.active.ident
rownames(s.obj)

DimPlot(annotated.seurat.obj, group.by = "celltype", cols = cbPalette10D, label = T)

#-------------------------------------------------------------------------------------------
## UMAPs and barplots


DimPlot(seurat.obj, cols = cbPalette10C) +
  theme_void() +
  NoLegend()

DimPlot(seurat.obj, group.by = "location", shuffle = T) +
  theme_void() +
  NoLegend()



#-------------------------------------------------------------------------------------------
## VIOLIN AND FEATURE PLOTS FOR CELL TYPE IDENTIFICATION


## ASTROCYTES
astro_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("AQP4", "GFAP", "Gja1", "Agt", "Slc1a2", "Slc1a3"), prefix, dataset), pt.size = 0, cols = cbPalette10C)
astro_violinplot

astro_featplot <- FeaturePlot(seurat.obj, features = convert_markers(c("AQP4", "GFAP", "Gja1", "Agt", "Slc1a2", "Slc1a3"), prefix, dataset))
astro_featplot



## OLIGODENDROCYTE + OPC
oligoOPC_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("MOG", "MBP", "Opalin", "Pdgfra", "Plp1", "Sox10"), prefix, dataset), pt.size = 0, cols = cbPalette10C)
oligoOPC_violinplot




## OLIGODENDROCYTE
oligo_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("MOG", "MBP", "Opalin"), prefix, dataset), pt.size = 0, cols = cbPalette10C)
oligo_violinplot


oligo_violinplot <- VlnPlot(seurat.obj, features = convert_markers(Oligo.M, prefix, dataset), pt.size = 0, cols = cbPalette10A)
oligo_violinplot

oligo_featplot <- FeaturePlot(seurat.obj, features = convert_markers(Oligo.M, prefix, dataset))
oligo_featplot


## OPC
opcs_violinplot <- VlnPlot(seurat.obj, features = convert_markers(OPC.M, prefix, dataset), pt.size = 0, cols=cbPalette10C)
opcs_violinplot

opcs_featplot <- FeaturePlot(seurat.obj, features = convert_markers(c("S100b", "Mbp", "Sox10", "Pdgfra", "Plp1", "Cspg4"), prefix, dataset))
opcs_featplot



## MICROGLIA
microglia_violinplot <- VlnPlot(seurat.obj, features = convert_markers(Microglia.M, prefix, dataset), pt.size = 0)
microglia_violinplot

microglia_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("Hexb", "Sall1", "Cx3cr1", "Tgfbr1", "C1qa", "Cd84"), prefix, dataset), pt.size = 0, cols = cbPalette10C)
microglia_violinplot

microglia_featplot <- FeaturePlot(seurat.obj, features = convert_markers(c("Tgfbr1", "Cd84", "Trem2", "Aif1", "Tyrobp", "Cldn5"), prefix, dataset))
microglia_featplot



## MSN-D1
msn1_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("Tac1", "Drd1"), prefix, dataset), pt.size = 0)
msn1_violinplot
msn_featplot <- FeaturePlot(seurat.obj, features = convert_markers(rat.msn.markers, prefix, dataset))
msn_featplot



## MSN-D2
msn2_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("Drd2", "Adora2a", "Penk", "Gpr6"), prefix, dataset), pt.size = 0)
msn_violinplot
msn_featplot <- FeaturePlot(seurat.obj, features = convert_markers(rat.msn.markers, prefix, dataset))
msn_featplot


## INTERNEURONS
interneurons_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("Pvalb", "Npy", "Lhx6", "Satb1", "Sox6", "Kcnc2", "Kcnc1", "NKX2-1", "Kcna1", "Kcna2", "Kcnab1", "Caln1"), prefix, dataset), pt.size = 0)

GABAergic_interneurons_violinplot <-  VlnPlot(seurat.obj, features = convert_markers(c("Pthlh","Cox6a2", "Opn3", "Chrna3","Gfra2"), prefix, dataset), pt.size = 0)

## To see the difference between culster 0 and 1
a <- FindMarkers(seurat.obj, ident.1 = "Interneurons")
a <- a %>% arrange(desc(avg_log2FC))
b <- a
b$spec <- b$pct.1 -b$pct.2
b <- b %>%filter(spec > 0.5)
b
c <- rownames(b) 
c


## GLUTAMINERGIC



## ENDOTHELIA
endothelia_violinplot <- VlnPlot(seurat.obj, features = c(convert_markers(PericytesEndotelial.M, prefix, dataset), "premRNARnor6--Cdh5", "premRNARnor6--Pecam1", "premRNARnor6--Ocln", "premRNARnor6--Esam", "premRNARnor6--Ng2", "premRNARnor6--Anpep", "premRNARnor6--Ccl4", "premRNARnor6--Cd53", "premRNARnor6--Cd258"), pt.size = 0)
endothelia_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("Cldn5", "Slc2a1", "Flt1", "Ocln", "Esam"), prefix, dataset), pt.size = 0)
endothelia_violinplot

## PERICYTE
pericyte_violinplot <- VlnPlot(seurat.obj, features = convert_markers(c("Rgs5", "Esam"), prefix, dataset), pt.size = 0)




#-------------------------------------------------------------------------------------------
## Dotplot with marker expression per cluster

DimPlot(seurat.obj, cols = cbPalette10A)

feat.dotplot <- c("Tac1", "Drd1", "Drd2", "Adora2a", "Penk", "Gpr6", "Cldn5", "Kcnc1", "Kcnab1", "Sox6", "Lhx6", "Slc2a1", "Flt1", "Esam", "Rgs5", "MOG", "MBP", "S100b", "Sox10", "Pdgfra", "Plp1", "Cspg4", "AQP4", "GFAP", "Gja1", "Agt", "Slc1a2", "Slc1a3", "Tgfbr1", "Cd84", "Trem2", "Aif1", "Tyrobp")

feat.dotplot <- c("MOG", "MBP", "Plp1", "Opalin", "Sox10", "Pdgfra", "Cspg4", "Neu4", "Tns3", "Fyn", "Bmp4", "S100b", "AQP4", "GFAP", "Gja1", "Agt", "Slc1a2", "Slc1a3", "Tgfbr1", "Cd84", "Hexb", "Trem2", "Aif1", "Tyrobp", "Tac1", "Drd1", "Gpr6", "SLC1A1", "SLC17A6", "Drd2", "Adora2a", "Penk", "Kcnab1", "Kcnc1", "Sox6", "Lhx6", "Cldn5", "Slc2a1", "KDR", "Flt1", "VIM", "LAMC1", "Esam", "Rgs5")

DotPlot(annotated.seurat.obj, features = convert_markers(feat.dotplot, prefix, dataset), cols = cbPalette10D, dot.min = 0.1)+
  RotatedAxis()


#### WORKS!!!!!
DotPlot(s.obj, features = str_to_title(feat.dotplot), group.by = "celltype", split.by = "celltype", cols = cbPalette10D, dot.min = 0.1)+
  RotatedAxis()

DotPlot(annotated.seurat.obj, features = convert_markers(feat.dotplot, prefix, dataset), group.by = "celltype", split.by = "celltype", cols = cbPalette10C, dot.min = 0.1)+
  RotatedAxis()

DotPlot(annotated.seurat.obj, features = convert_markers(feat.dotplot, prefix, dataset), split.by = factor("celltype", levels = celltype.order2), cols = cbPalette10D, dot.min = 0.1)+
  RotatedAxis()


feat.glia.dotplot <- c("MOG", "MBP", "S100b", "Olig1", "Olig2", "Nkx2-2", "Sox10", "Pdgfra", "Plp1", "Cspg4", "AQP4", "GFAP", "Gja1", "Agt", "Slc1a2", "Slc1a3", "Tgfbr1", "Cd84", "Trem2", "Hexb", "Entpd1", "Cx3cr1", "Olfml3")

DotPlot(seurat.obj, features = convert_markers(feat.glia.dotplot, prefix, dataset), split.by = "location", cols = c("#F8766D", "#00BFC4"), idents = c("Astrocyte", "Microglia", "OPC", "Oligodendrocyte"), dot.min = 0.1)+
  RotatedAxis()
  
DotPlot(seurat.obj, features = convert_markers(feat.glia.dotplot, prefix, dataset), cols = cbPalette10C, idents = c("Astrocyte", "Microglia", "OPC", "Oligodendrocyte"), dot.min = 0.1)+
  RotatedAxis()




#-----------------------------------------------------------------------------------------------

## CELL TYPE PER SAMPLE STACKED BAR PLOT

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

celltype.order <- c("MSN-D1", "MSN-D2", "Interneurons", "Glutaminergic", "Oligodendrocyte", "OPC", "Astrocyte", "Microglia", "Endothelia", "Pericyte")


#Plot the stacked bar graph % cell type per sample
CellTypePerSample.StackedBarPlot <- ggplot(CellTypePerSample.table, aes(x = Var2, y = percentage)) +
  facet_grid(~ Var3, scales = "free_x", space = "free_x") +
  geom_col(aes(fill = factor(Var1, levels = celltype.order2)), position="fill") +
  scale_fill_manual(values = cbPalette10D) +
  theme_void() +
  NoLegend()


#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_POSTER_CellTypePerCellType_StackedBarPlot.png"), width = 1200, height = 700)
print(CellTypePerAge.StackedBarPlot)
dev.off()



#---------------------------------------------------------------------------------------------
## CELL TYPE PER AGE

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
CellTypePerAge.StackedBarPlot <- ggplot(CellTypePerAge.table, aes(x = Var2, y = Freq, fill = factor(Var1, levels = celltype.order2))) +
  facet_grid(~ Var3, scales = "free_x", space = "free_x") +
  geom_col(width = 0.5, position = "fill") +
  scale_fill_manual(values = cbPalette10D) +
  theme_void() +
  NoLegend()


#save the plot in a png image
png(filename=paste0(plot.dir, dataset, "_POSTER_CellTypePerAge_StackedBarPlot.png"), width = 1200, height = 700)
print(CellTypePerAge.StackedBarPlot)
dev.off()



