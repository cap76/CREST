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
Idents(mammal.combined2) <- mammal.combined2$Dataset
mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")

DefaultAssay(mammal.combined2) <- "RNA"
FeaturePlot(mammal.combined2,feature="ISL1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_ISL1_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="DKK1",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DKK1_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="SFRP1",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SFPR1_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SFRP4",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SFRP4_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DKK2",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DKK2_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="IGFBP1",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_IGFBP1_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="IGFBP3",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_IGFBP3_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="IGFBP4",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_IGFBP4_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="IGFBP2",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_IGFBP2_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
#FeaturePlot(mammal.combined2,feature="IGFBP4",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
#ggsave(filename=paste(saveext,"FinalPCA_IGFBP4_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

#FeaturePlot(mammal.combined2,feature="IGFBP7",split.by = "Dataset",reduction = "umap", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
#ggsave(filename=paste(saveext,"FinalPCA_IGFBP7_justEmbryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


sdsadsadsad

#X1 <- GetAssayData(D1)
#saveRDS(X1,file="Tb_expession.rds")
#saveRDS(Idents(D1),file="Tb_expession_label.rds")
#saveRDS(D1$Dataset,file="Tb_expession_Dataset.rds")

mammal.combined2 <- subset(mammal.combined2,idents=c("Unciliated epithelia","Ciliated"))
mammal.combined2$Cells <- Idents(mammal.combined2)
Idents(mammal.combined2) <- mammal.combined2$Dataset
mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")

Idents(mammal.combined2) <- mammal.combined2$Cells
DefaultAssay(mammal.combined2) <- "RNA"
Type1<-WhichCells(mammal.combined2,expression=SPP1>0.250961574,slot="data")
Type2<-WhichCells(mammal.combined2,expression=SCGB2A2>0.250961574,slot="data")
Gland <- intersect(Type1,Type2)
Idents(mammal.combined2,cells=Gland) <- "Glandular"

saveRDS(Gland,file="GlandIDs.rds")
DefaultAssay(mammal.combined2) <- "RNA"

FeatureScatter(mammal.combined2,feature1="PAX2",feature2="SCGB2A2")
ggsave(filename=paste(saveext,"FinalFS_PAX2_SCGB2A2",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2,feature1="SCGB2A2",feature2="SPP1")
ggsave(filename=paste(saveext,"FinalFS_PAX2_SPP1",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2,feature1="SCGB2A2",feature2="DPP4")
ggsave(filename=paste(saveext,"FinalFS_PAX2_DPP4",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2,feature1="SCGB2A2",feature2="GPX3")
ggsave(filename=paste(saveext,"FinalFS_PAX2_GPX3",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeatureScatter(mammal.combined2,feature1="SOX9",feature2="LGR5")
ggsave(filename=paste(saveext,"FinalFS_SOX2_LGR5",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)




mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
mammal.combined2$Cells <- Idents(mammal.combined2)

mammal.combined2 <- subset(mammal.combined2,idents=c("EmDisc1_d14","EmDisc2_d14","Am_d14","Am/EmDisc_d14"))
mammal.combined2$Cells <- Idents(mammal.combined2)

DefaultAssay(mammal.combined2) <- "RNA"

FeatureScatter(mammal.combined2,feature1="PDGFA",feature2="SOX15")
ggsave(filename=paste(saveext,"FinalFS_PDGFA_SOX152",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2,feature1="PDGFA",feature2="POU5F1")
ggsave(filename=paste(saveext,"FinalFS_PDGFA_POU5F1",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2,feature1="WNT6",feature2="VTCN1")
ggsave(filename=paste(saveext,"FinalFS_WNT6_VTCN1",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)




mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
#mammal.combined2$Cells <- Idents(mammal.combined2)
Idents(mammal.combined2) <- paste(mammal.combined2$Dataset,Idents(mammal.combined2),sep="_")

DefaultAssay(mammal.combined2) <- "RNA"
DE4 <- FindMarkers(mammal.combined2,ident.1 = "10X Ours_ExMes_d14", ident.2 = "10X Ours_Stromal fibroblasts",only.pos=TRUE,test.use="MAST")

write.table(as.data.frame(DE4),file="ExMes_DE_vs_Stroma.csv")

#DE1 <- 
#FeatureScatter(mammal.combined2,feature1="",feature2="SCGB2A2")
#ggsave(filename=paste(saveext,"FinalFS_PAX2_SCGB2A2",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

#FeatureScatter(mammal.combined2,feature1="PAX2",feature2="SCGB2A2")
#ggsave(filename=paste(saveext,"FinalFS_PAX2_SCGB2A2",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)




