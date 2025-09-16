library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("Matrix")
library(pheatmap)
#library(destiny)
set.seed(1)

#Save folder
saveext = "./FinalAlignB/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))

mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
mammal.combined2$Cells <- Idents(mammal.combined2)
D1 <- subset(mammal.combined2,idents=c("STB_d14","CTB_d14","EVT_d14","STB","EVT","CTB"))

Idents(D1) <- paste(D1$Dataset,D1$Cl05,sep="_")


DefaultAssay(D1) <- "RNA"

M1 <- FindMarkers(D1, ident.1=c("10X Ours_20","10X Ours_2","10X Ours_23"), ident.2=c("10X Ours_14","10X Ours_11","10X Ours_21"), verbose = FALSE, test.use = "MAST", only.pos = TRUE) #, logfc.threshold = log(0))

write.table(as.data.frame(M1),file="DEExpressionTb.csv")
