library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
library(stringr)

set.seed(1) 

saveext = "~/Desktop/Data/Endometrial/InVitro/Matteo/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#Get the dataset
D <- readRDS("/Users/christopherpenfold/Downloads/CREST_anotated.rds")

#This object contains more fine anotations, lets set to the broad anotations for Figure 6
Idents(D,cells=WhichCells(D,idents=c("Stroma_decid","Stroma_prolif"))) <- "Stroma"
Idents(D,cells=WhichCells(D,idents=c("Epithelia_gland","Epithelia_lum","Epithelia_SOX"))) <- "Unciliated epithelia"
Idents(D,cells=WhichCells(D,idents=c("Epithelia_cil"))) <- "Ciliated epithelia"
Idents(D,cells=WhichCells(D,idents=c("EmbDisc"))) <- "Embryonic Disc"
Idents(D,cells=WhichCells(D,idents=c("EVT"))) <- "Extravillous trophoblast"
Idents(D,cells=WhichCells(D,idents=c("STB","eSTB"))) <- "Syncitiotrophoblast"
Idents(D,cells=WhichCells(D,idents=c("CTB"))) <- "Cytotrophoblast"
Idents(D,cells=WhichCells(D,idents=c("ExMes"))) <- "Mesenchyme"

#Colours for plotting
cols = c('Unciliated epithelia'='#c4807f',
'Ciliated epithelia'='#7f3c3c',
'Stroma'='#ad71bc',
'Embryonic Disc'='#2d7ee0',
'Amnion'='#3cdddd',
'Hypoblast'='#ff8c64',
'Extravillous trophoblast'='#bf3737',
'Syncitiotrophoblast'='#cb8f46',
'Cytotrophoblast'='#e2c120',
'YolkSac'='#2e892e',
'Mesenchyme'='#9fd13f')

#Plot the UMAP.
#This is the UMAP aligned to the composite reference (Wang et al, and Mole et al) used in Figure 1
p<-DimPlot(D, cols=cols,  pt.size = 4, reduction = "oldUMAP",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Figure6A.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

#This is an alternative UMAP based on harmony aligned co-ordinates (using just our data)
p<-DimPlot(D, cols=cols, pt.size = 4, reduction = "harmumap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/HarmonisedUMAP_split.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

#This is a UMAP based on the PC from a seurat integrated alignment using just our data
p<-DimPlot(D, cols=cols, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SeuratAlignedUMAP_split.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


#Panel of markers for visualisation
Markers1 <- c("SOX15","CD9","PDGFA","SERPINE2","SFRP1","GAL","DNAJC15","L1TD1","PODXL","CDH11","GREB1L","PATJ","FUT8","DMD","VTCN1",
              "APOC3","RSPO3","PDGFRA","FRZB","VCAN","GATA6","HGF","APOB","LRP2","GPC3","COL18A1","TFF3","GATA4","BST2","APOE","CER1","LEFTY1","LEFTY2",
              "DPP4","APELA","OPHN1","TEAD3","GATA2","ERVV-1","CGB8","CGB5","CGA","SDC1","TBX3","DIO2","NOTUM","ASCL2","HLA-G","SNAI1")

#Calculate aveage expression. Alternative approaches now exist for this.
Av <- AverageExpression(D)
Z1 <- Av$RNA[Markers1,c("Embryonic Disc","Amnion","YolkSac","Mesenchyme","Hypoblast","Hypoblast","Cytotrophoblast","Syncitiotrophoblast","Extravillous trophoblast")]
Z1 <- t(scale(t(Z1) ))
Z1 <- Z1[rowSums(is.na(Z1)) == 0, ]

mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))

pheatmap(Z1, breaks=mat_breaks, color =  redblue1(20), border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Figure6B.pdf",sep=""),width=10,height=22)
