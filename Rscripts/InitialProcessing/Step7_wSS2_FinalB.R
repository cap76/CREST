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
Idents(mammal.combined2) <- paste(mammal.combined2$Dataset,Idents(mammal.combined2),sep="_")

AE1 <- AverageExpression(mammal.combined2)
Ae1 <- AE1$RNA


other_endometrium_data <- readRDS(file="../../Endometrial/OtherEndo/other_endometrium_data.rds")
other_endometrium_data <- subset(other_endometrium_data,idents=c("eS","dS","SOX9","Ciliated","Lumenal","Fibroblast C7","Glandular"))
Idents(other_endometrium_data) <- other_endometrium_data$ID3
other_endometrium_data$Dataset <- "Ref"
other_endometrium_data$Dataset2 <- other_endometrium_data$ID3
other_endometrium_data$ID3 <- other_endometrium_data$ID3

Idents(other_endometrium_data) <- paste(other_endometrium_data$ID3,other_endometrium_data$ID5,sep="_")


other_endometrium_data <- FindVariableFeatures(other_endometrium_data,nfeatures = 2000)
feats <- other_endometrium_data@assays$RNA@var.features

other_endometrium_data <- FindVariableFeatures(other_endometrium_data,nfeatures = 5000)
feats2 <- other_endometrium_data@assays$RNA@var.features

other_endometrium_data <- FindVariableFeatures(other_endometrium_data,nfeatures = 10000)
feats3 <- other_endometrium_data@assays$RNA@var.features

other_endometrium_data <- FindVariableFeatures(other_endometrium_data,nfeatures = 15000)
feats4<- other_endometrium_data@assays$RNA@var.features

AE2 <- AverageExpression(other_endometrium_data)
Ae2 <- AE2$RNA


other_endometrium_data1 <- readRDS(file="../../Endometrial/OtherEndo/other_endometrium_data.rds")
uID <- as.character(Idents(other_endometrium_data1))
uID[which(uID%in%c("Fibroblast C7","eS","dS"))] <- "Stroma"
uID[which(uID%in%c("SOX9","Lumenal","Glandular"))] <- "Eptihelial"
uID[which(uID%in%c("Ciliated"))] <- "Ciliated"
Idents(other_endometrium_data1) <- paste(other_endometrium_data1$ID3,uID,sep="_")

AE3 <- AverageExpression(other_endometrium_data1)
Ae3 <- AE3$RNA

list1 <- c("10X Ours_Stromal fibroblasts","10X Ours_Unciliated epithelia","10X Ours_Ciliated")

list2 <- c("proliferative_Fibroblast C7",
"early-secretory_Fibroblast C7",
"early-mid-secretory_Fibroblast C7",
"mid-secretory_Fibroblast C7",
"late-secretory_Fibroblast C7",
"proliferative_eS",
"early-secretory_eS",
"early-mid-secretory_eS",
"mid-secretory_eS",
"late-secretory_eS",
"proliferative_dS",
"early-secretory_dS",
"early-mid-secretory_dS",
"mid-secretory_dS",
"late-secretory_dS",
"proliferative_SOX9",
"early-secretory_SOX9",
"early-mid-secretory_SOX9",
"mid-secretory_SOX9",
"late-secretory_SOX9",
"proliferative_Lumenal",
"early-secretory_Lumenal",
"early-mid-secretory_Lumenal",
"mid-secretory_Lumenal",
"late-secretory_Lumenal",
"proliferative_Glandular",
"early-secretory_Glandular",
"early-mid-secretory_Glandular",
"mid-secretory_Glandular",
"late-secretory_Glandular",
"proliferative_Ciliated",
"early-secretory_Ciliated",
"early-mid-secretory_Ciliated",
"mid-secretory_Ciliated",
"late-secretory_Ciliated")

list3 <- c("proliferative_Stroma",
"early-secretory_Stroma",
"early-mid-secretory_Stroma",
"mid-secretory_Stroma",
"late-secretory_Stroma",
"proliferative_Eptihelial",
"early-secretory_Eptihelial",
"early-mid-secretory_Eptihelial",
"mid-secretory_Eptihelial",
"late-secretory_Eptihelial", 
"proliferative_Ciliated",
"early-secretory_Ciliated",
"early-mid-secretory_Ciliated",
"mid-secretory_Ciliated",
"late-secretory_Ciliated")

GL <- intersect(rownames(GetAssayData(mammal.combined2,assay="integrated")),
rownames(GetAssayData(other_endometrium_data,assay="RNA")))

D1 <- cor(log2(Ae1[GL,list1]+1),log2(Ae2[GL,list2]+1))
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D1),color =  redblue1(20),  border_color = NA, display_numbers = round(D1, digits = 2), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/CrossCorb_stroma.pdf",sep=""),width=20,height=5)


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
