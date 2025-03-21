library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("Matrix")
library(pheatmap)
set.seed(1)

#Save folder
saveext = "./FinalAlignB/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))


mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
Idents(mammal.combined2) <- mammal.combined2$Dataset
mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")
Idents(mammal.combined2) <- paste(mammal.combined2$ID3,mammal.combined2$Cl05,sep="_")

AE1 <- AverageExpression(mammal.combined2)
Ae1 <- AE1$RNA


GL <- rownames(GetAssayData(mammal.combined2,assay="integrated"))

D1 <- cor(log2(Ae1[GL,]+1),log2(Ae1[GL,]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/CrossCorb_Selfcor.pdf",sep=""),width=20,height=20)

mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))

mammal.combined2$Cells <- Idents(mammal.combined2)
Idents(mammal.combined2) <- mammal.combined2$Dataset
mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")
Idents(mammal.combined2) <- paste(mammal.combined2$ID3,mammal.combined2$Cells,sep="_")

AE1 <- AverageExpression(mammal.combined2)
Ae1 <- AE1$RNA


GL <- rownames(GetAssayData(mammal.combined2,assay="integrated"))

D1 <- cor(log2(Ae1[GL,]+1),log2(Ae1[GL,]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/CrossCorb_Selfcor2.pdf",sep=""),width=20,height=20)


sadadasd

D1 <- cor(log2(Ae1[GL,list1]+1),log2(Ae3[GL,list3]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, display_numbers = round(D1, digits = 2), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/CrossCorb_stromaB.pdf",sep=""),width=20,height=5)



GL2 <- intersect(rownames(GetAssayData(mammal.combined2,assay="RNA")),feats)
D1 <- cor(log2(Ae1[GL2,list1]+1),log2(Ae2[GL2,list2]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, display_numbers = round(D1, digits = 2), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/CrossCorb_stroma2.pdf",sep=""),width=20,height=5)

GL3 <- intersect(rownames(GetAssayData(mammal.combined2,assay="RNA")),feats2)
D1 <- cor(log2(Ae1[GL3,list1]+1),log2(Ae2[GL3,list2]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, display_numbers = round(D1, digits = 2), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/CrossCorb_stroma3.pdf",sep=""),width=20,height=5)

GL4 <- intersect(rownames(GetAssayData(mammal.combined2,assay="RNA")),feats3)
D1 <- cor(log2(Ae1[GL4,list1]+1),log2(Ae2[GL4,list2]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, display_numbers = round(D1, digits = 2), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/CrossCorb_stroma4.pdf",sep=""),width=20,height=5)

GL5 <- intersect(rownames(GetAssayData(mammal.combined2,assay="RNA")),feats4)
D1 <- cor(log2(Ae1[GL5,list1]+1),log2(Ae2[GL5,list2]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, display_numbers = round(D1, digits = 2), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/CrossCorb_stroma5.pdf",sep=""),width=20,height=5)



sdsadsadada
