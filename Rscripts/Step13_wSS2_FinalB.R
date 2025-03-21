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

mammal.combined3 <- readRDS(file="JustEmbryonic.rds")
l1<-colnames(mammal.combined3)[which(mammal.combined3$newCl=="16")]
l2<-colnames(mammal.combined3)[which(mammal.combined3$newCl=="4")]
l3<-colnames(mammal.combined3)[which(mammal.combined3$newCl%in%c("18","9") )]
l3A<-colnames(mammal.combined3)[which(Idents(mammal.combined3)%in%c("ExMes_d14") )]
l4<-colnames(mammal.combined3)[which(mammal.combined3$newCl%in%c("1","11") )]


#X1 <- GetAssayData(D1)
#saveRDS(X1,file="Tb_expession.rds")
#saveRDS(Idents(D1),file="Tb_expession_label.rds")
#saveRDS(D1$Dataset,file="Tb_expession_Dataset.rds")

#mammal.combined2 <- subset(mammal.combined2,idents=c("Unciliated epithelia","Ciliated"))
#mammal.combined2$Cells <- Idents(mammal.combined2)
#Idents(mammal.combined2) <- mammal.combined2$Dataset
#mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")

Idents(mammal.combined2) <- mammal.combined2$Cells
DefaultAssay(mammal.combined2) <- "RNA"
Type0<-WhichCells(mammal.combined2,idents=c("Unciliated epithelia","Ciliated"))

Type1<-WhichCells(mammal.combined2,expression=SPP1>0.250961574,slot="data")
Type2<-WhichCells(mammal.combined2,expression=SCGB2A2>0.250961574,slot="data")
Gland <- intersect(intersect(Type1,Type2),Type0)
Idents(mammal.combined2,cells=Gland) <- "Glandular"

#saveRDS()

Idents(mammal.combined2) <- mammal.combined2$Cl05
Type3<-WhichCells(mammal.combined2,idents=c(1,5))
Type3B<-WhichCells(mammal.combined2,idents=c(0))
Type3C<-WhichCells(mammal.combined2,idents=c(22,13,12,8,4,3))
Type3D<-WhichCells(mammal.combined2,idents=c(2))
Type3E<-WhichCells(mammal.combined2,idents=c(22))
Type3F<-WhichCells(mammal.combined2,idents=c(4,3,19,12))
Type3G<-WhichCells(mammal.combined2,idents=c(13))


Idents(mammal.combined2) <- mammal.combined2$Dataset
Type4<-WhichCells(mammal.combined2,idents=c("10X Ours"))

Idents(mammal.combined2) <- mammal.combined2$Cells
#Idents(mammal.combined2,cells=intersect(Type3,Type4)) <- "Glandular"
#Idents(mammal.combined2,cells=Gland) <- "Glandular"

#Idents(mammal.combined2) <- mammal.combined2$Cells
Idents(mammal.combined2,cells=intersect(Type3C,Type4)) <- "Prolif"
Idents(mammal.combined2,cells=intersect(Type3B,Type4)) <- "Prolif"
Idents(mammal.combined2,cells=intersect(Type3F,Type4)) <- "Prolif"

Idents(mammal.combined2,cells=intersect(Type3D,Type4)) <- "Lumenal"

Type5<-WhichCells(mammal.combined2,expression=SOX9>0.250961574,slot="data")
Type6<-WhichCells(mammal.combined2,expression=LGR5>0.250961574,slot="data")

Idents(mammal.combined2,cells=intersect(Type4,intersect(Type5,Type5))) <- "SOX9P"
Idents(mammal.combined2,cells=intersect(Type4,intersect(Type5,Type6))) <- "SOX9LRG5"

#Idents(mammal.combined2,cells=intersect(Type4,intersect(Type5,Type6))) <- "Ciliated"

#Idents(mammal.combined2,cells=intersect(Type3,Type4)) <- "Glandular"
Idents(mammal.combined2,cells=Gland) <- "Glandular"

notEVT <- readRDS("notEVT.rds")

EVT <- readRDS("newEVT.rds")
Idents(mammal.combined2,cells=EVT) <- "EVT_d14"
Idents(mammal.combined2,cells=notEVT) <- "STB_d14"

Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="EmDisc1_d14")) <- "EmDisc_d14"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents="EmDisc2_d14")) <- "EmDisc_d14"


#Idents(mammal.combined2,cells=l1) <- "EmDisc_d14"
#Idents(mammal.combined2,cells=l2) <- "EmD4_d14"
Idents(mammal.combined2,cells=intersect(l3,l3A)) <- "ExMes_d14"
#Idents(mammal.combined2,cells=intersect(l4,l3A)) <- "Mes2_d14"


saveRDS(mammal.combined2,file=paste(saveext,"FinalAnotatedMergedSeuratObj.rds",sep=""))
saveRDS(Idents(mammal.combined2),file=paste(saveext,"FinalAnotatedMergedSeuratObjMeta.rds",sep=""))


#Which()


mammal.combined2EVT <- subset(mammal.combined2,idents=c("EVT_d14"))
Idents(mammal.combined2EVT) <- colnames(mammal.combined2EVT)
p<-DimPlot(mammal.combined2EVT, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"PCA_Seurat-AllAnoEVT",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined2EVT, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"UMAP_Seurat-AllAnoEVT",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE,p)



p<-DimPlot(mammal.combined2, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-AllAno",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined2, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-AllAno",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE,p)


mammal.combined3 <- mammal.combined2
Idents(mammal.combined3) <- mammal.combined3$Cl05
p<-DimPlot(mammal.combined3, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-AllCl",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined3, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-AllCl",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE,p)



#Av <- AverageExpression(mammal.combined2)
#Av <- Av$RNA


#DE1 <- FindMarkers(mammal.combined2,ident.1 = "10X Ours_EVT_d14", ident.2 = "10X Ours_Stromal fibroblasts",only.pos=TRUE,test.use="MAST")
#DE2 <- FindMarkers(mammal.combined2,ident.1 = "10X Ours_EVT_d14", ident.2 = "10X Ours_Stromal fibroblasts",only.pos=TRUE,test.use="MAST")
#DE3 <- FindMarkers(mammal.combined2,ident.1 = "10X Ours_STB_d14", ident.2 = "10X Ours_Stromal fibroblasts",only.pos=TRUE,test.use="MAST")
#DE4 <- FindMarkers(mammal.combined2,ident.1 = "10X Ours_STB_d14", ident.2 = "10X Ours_Stromal fibroblasts",only.pos=TRUE,test.use="MAST")
#DE5 <- FindMarkers(mammal.combined2,ident.1 = "10X Ours_ExMes_d14", ident.2 = "10X Ours_Stromal fibroblasts",only.pos=TRUE,test.use="MAST")
#DE6 <- FindMarkers(mammal.combined2,ident.1 = "10X Ours_ExMes_d14", ident.2 = "10X Ours_Stromal fibroblasts",only.pos=TRUE,test.use="MAST")


#asdadada
#Idents(mammal.combined2) <- mammal.combined2$Cl05
#p<-DimPlot(mammal.combined2, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_Seurat-AllAnoCl",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE,p)


DefaultAssay(mammal.combined2) <- "RNA"

Markers <- c("JAM3",
"GATA2",
"GATA3",
"TFAP2A",
"TFAP2C",
"CGB8",
"CGB5",
"CGA",
"ERVW-1",
"PRDM6",
"TBX3",
"DIO2",
"NOTUM",
"HLA-G",
"ASCL2",
"SNAI1",
"HGF",
"SNAI2",
"HAND1",
"HAND2",
"PDGFRA",
"BST2",
"GATA6",
"GATA4",
"CER1",
"NODAL",
"LEFTY1",
"LEFTY2",
"APOA1",
"BAMBI",
"WNT6",
"VTCN1",
"AKAP12",
"PODXL",
"SFRP1",
"SOX15",
"PDGFA",
"POU5F1")

Av <- AverageExpression(mammal.combined2)
Av <- Av$RNA

saveRDS(Av,file="Avt.rds")
saveRDS(setdiff(c("CTB_d14","STB_d14","EVT_d14","ExMes_d14","Hyp_d14","Am_d14","Am/EmDisc_d14","EmDisc_d14"),colnames(Av)),file="Setdiff.rds")
D <- Av[Markers,c("CTB_d14","STB_d14","EVT_d14","ExMes_d14","Hyp_d14","Am_d14","Am/EmDisc_d14","EmDisc_d14")]

mat_breaks <- seq(-3, 3, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D),color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers.pdf",sep=""),scale="row",width=10,height=12)
sadada
Idents(mammal.combined2) <- paste(mammal.combined2$Dataset,Idents(mammal.combined2),sep="_")
Av <- AverageExpression(mammal.combined2)
Av <- Av$RNA

DE2<-readRDS(file="DE_Cil_Lum1.rds")
DE2A<-readRDS(file="DE_Cil_Lum2.rds")

DE3<-readRDS(file="DE_Lum_Gla1.rds")
DE3A<-readRDS(file="DE_Lum_Gla2.rds")

DE4<-readRDS(file="DE_Cil_Gla1.rds")
DE4A<-readRDS(file="DE_Cil_Gla2.rds")


TF<-read.table("../../Thorsten/Other/TF.txt",header = F)
TF <- TF$V1
SIGNAL<-read.table("../../Thorsten/Other/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]

Markers <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
list2 <- intersect(c(TF,SIGNAL1,SIGNAL2,SIGNAL3,Markers),rownames(Av))
#write.table(as.data.frame(Cl),file="ClustersEmbryonicOnly.csv",quote = FALSE, sep = ",")
AvExp3 <- as.data.frame(Av)
AvExp3$Sum <- AvExp3[,"10X Ours_Prolif"]+AvExp3[,"10X Ours_Glandular"]+AvExp3[,"10X Ours_Lumenal"]
AvExp3$Delta <- log(apply(AvExp3[,c("10X Ours_Prolif","10X Ours_Glandular","10X Ours_Lumenal")], 1, max)+1)-log(apply(AvExp3[,c("10X Ours_Prolif","10X Ours_Glandular","10X Ours_Lumenal")], 1, min) +1)
AvExp3$Delta <- 3*(AvExp3$Delta/max(AvExp3$Delta))
AvExp3$MAX <- 100*(apply(AvExp3[,c("10X Ours_Prolif","10X Ours_Glandular","10X Ours_Lumenal")], 1, max))
AvExp3$MAX[which(AvExp3$MAX<10)] <- 0.1
AvExp3$MAX[which(AvExp3$MAX>10)] <- 5
#Markers <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
#"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
#"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
#list2 <- rownames(DE4)

#FeatureScatter(mammal.combined2,feature1="SOX9",feature2="LGR5")
#ggsave(filename=paste(saveext,"FinalFS_SOX9_LRG5",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



AVE2 <-  as.data.frame(AvExp3[,c("10X Ours_Prolif","10X Ours_Glandular","10X Ours_Lumenal")])
colnames(AVE2) <- c("Prolif","Gland","Lumenal")
AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
AVE2<- AVE2[which(AVE2$MAX>1),]
AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
library(ggtern)
p1 <- ggtern(data=AVE2, aes(x =Prolif,z=Gland,y=Lumenal)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
p1 <- p1+ theme_showgrid()
p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Prolif)/(AVE2$Prolif+AVE2$Gland+AVE2$Lumenal),
z= (AVE2$Gland)/(AVE2$Prolif+AVE2$Gland+AVE2$Lumenal),
y= (AVE2$Lumenal)/(AVE2$Prolif+AVE2$Gland+AVE2$Lumenal),label = rownames(AVE2),color = c("black"))
ggsave(filename=paste(saveext,"GGTern_Markers_ProlifGlandLumenal.pdf",sep=""),width = 20, height = 20, plot = p1)


#sadddsadsadad

AvExp3 <- as.data.frame(Av)
AvExp3$Sum <- AvExp3[,"10X Ours_Unciliated epithelia"]+AvExp3[,"10X Ours_Glandular"]+AvExp3[,"10X Ours_Ciliated"]
AvExp3$Delta <- log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max)+1)-log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, min) +1)
AvExp3$Delta <- 3*(AvExp3$Delta/max(AvExp3$Delta))
AvExp3$MAX <- 100*(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max))
AvExp3$MAX[which(AvExp3$MAX<10)] <- 0.1
AvExp3$MAX[which(AvExp3$MAX>10)] <- 5
list2 <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
list2 <- rownames(DE2)
AVE2 <-  as.data.frame(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")])
colnames(AVE2) <- c("Epit","Gland","Cil")
AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
AVE2<- AVE2[which(AVE2$MAX>1),]
AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
library(ggtern)
p1 <- ggtern(data=AVE2, aes(x =Epit,z=Gland,y=Cil)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
p1 <- p1+ theme_showgrid()
p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Epit)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
z= (AVE2$Gland)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
y= (AVE2$Cil)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),label = rownames(AVE2),color = c("black"))
ggsave(filename=paste(saveext,"GGTern_Markers_Cil_v_Lum.pdf",sep=""),width = 20, height = 20, plot = p1)

AvExp3 <- as.data.frame(Av)
AvExp3$Sum <- AvExp3[,"10X Ours_Unciliated epithelia"]+AvExp3[,"10X Ours_Glandular"]+AvExp3[,"10X Ours_Ciliated"]
AvExp3$Delta <- log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max)+1)-log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, min) +1)
AvExp3$Delta <- 3*(AvExp3$Delta/max(AvExp3$Delta))
AvExp3$MAX <- 100*(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max))
AvExp3$MAX[which(AvExp3$MAX<10)] <- 0.1
AvExp3$MAX[which(AvExp3$MAX>10)] <- 5
list2 <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
list2 <- rownames(DE3)
AVE2 <-  as.data.frame(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")])
colnames(AVE2) <- c("Epit","Gland","Cil")
AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
AVE2<- AVE2[which(AVE2$MAX>1),]
AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
library(ggtern)
p1 <- ggtern(data=AVE2, aes(x =Epit,z=Gland,y=Cil)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
p1 <- p1+ theme_showgrid()
p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Epit)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
z= (AVE2$Gland)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
y= (AVE2$Cil)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),label = rownames(AVE2),color = c("black"))
ggsave(filename=paste(saveext,"GGTern_Markers_Lum_v_Gland.pdf",sep=""),width = 20, height = 20, plot = p1)


AvExp3 <- as.data.frame(Av)
AvExp3$Sum <- AvExp3[,"10X Ours_Unciliated epithelia"]+AvExp3[,"10X Ours_Glandular"]+AvExp3[,"10X Ours_Ciliated"]
AvExp3$Delta <- log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max)+1)-log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, min) +1)
AvExp3$Delta <- 3*(AvExp3$Delta/max(AvExp3$Delta))
AvExp3$MAX <- 100*(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max))
AvExp3$MAX[which(AvExp3$MAX<10)] <- 0.1
AvExp3$MAX[which(AvExp3$MAX>10)] <- 5
list2 <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
list2 <- rownames(DE4A)
AVE2 <-  as.data.frame(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")])
colnames(AVE2) <- c("Epit","Gland","Cil")
AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
AVE2<- AVE2[which(AVE2$MAX>1),]
AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
library(ggtern)
p1 <- ggtern(data=AVE2, aes(x =Epit,z=Gland,y=Cil)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
p1 <- p1+ theme_showgrid()
p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Epit)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
z= (AVE2$Gland)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
y= (AVE2$Cil)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),label = rownames(AVE2),color = c("black"))
ggsave(filename=paste(saveext,"GGTern_Markers_Gland_v_Cil.pdf",sep=""),width = 20, height = 20, plot = p1)

AvExp3 <- as.data.frame(Av)
AvExp3$Sum <- AvExp3[,"10X Ours_Unciliated epithelia"]+AvExp3[,"10X Ours_Glandular"]+AvExp3[,"10X Ours_Ciliated"]
AvExp3$Delta <- log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max)+1)-log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, min) +1)
AvExp3$Delta <- 3*(AvExp3$Delta/max(AvExp3$Delta))
AvExp3$MAX <- 100*(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max))
AvExp3$MAX[which(AvExp3$MAX<10)] <- 0.1
AvExp3$MAX[which(AvExp3$MAX>10)] <- 5
list2 <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
list2 <- rownames(DE2A)
AVE2 <-  as.data.frame(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")])
colnames(AVE2) <- c("Epit","Gland","Cil")
AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
AVE2<- AVE2[which(AVE2$MAX>1),]
AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
library(ggtern)
p1 <- ggtern(data=AVE2, aes(x =Epit,z=Gland,y=Cil)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
p1 <- p1+ theme_showgrid()
p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Epit)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
z= (AVE2$Gland)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
y= (AVE2$Cil)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),label = rownames(AVE2),color = c("black"))
ggsave(filename=paste(saveext,"GGTern_Markers_Lum_v_Cil.pdf",sep=""),width = 20, height = 20, plot = p1)

AvExp3 <- as.data.frame(Av)
AvExp3$Sum <- AvExp3[,"10X Ours_Unciliated epithelia"]+AvExp3[,"10X Ours_Glandular"]+AvExp3[,"10X Ours_Ciliated"]
AvExp3$Delta <- log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max)+1)-log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, min) +1)
AvExp3$Delta <- 3*(AvExp3$Delta/max(AvExp3$Delta))
AvExp3$MAX <- 100*(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max))
AvExp3$MAX[which(AvExp3$MAX<10)] <- 0.1
AvExp3$MAX[which(AvExp3$MAX>10)] <- 5
list2 <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
list2 <- rownames(DE3A)
AVE2 <-  as.data.frame(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")])
colnames(AVE2) <- c("Epit","Gland","Cil")
AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
AVE2<- AVE2[which(AVE2$MAX>1),]
AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
library(ggtern)
p1 <- ggtern(data=AVE2, aes(x =Epit,z=Gland,y=Cil)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
p1 <- p1+ theme_showgrid()
p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Epit)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
z= (AVE2$Gland)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
y= (AVE2$Cil)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),label = rownames(AVE2),color = c("black"))
ggsave(filename=paste(saveext,"GGTern_Markers_Gland_v_Lum.pdf",sep=""),width = 20, height = 20, plot = p1)



mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
mammal.combined2$Cells <- Idents(mammal.combined2)
Idents(mammal.combined2) <- mammal.combined2$Cells
mammal.combined2 <- subset(mammal.combined2,idents=c("Unciliated epithelia","Ciliated"))

DefaultAssay(mammal.combined2) <- "RNA"

Idents(mammal.combined2) <- paste(mammal.combined2$Dataset,mammal.combined2$Cl05,sep="_")

Markers <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PAX2","VTCN1","SLC26A7","MSLN")

Av <- AverageExpression(mammal.combined2)
Av <- Av$RNA
D <- Av[Markers,c("10X Ours_0","10X Ours_19","10X Ours_22","10X Ours_12","10X Ours_3","10X Ours_8","10X Ours_5","10X Ours_1","10X Ours_13","10X Ours_4","10X Ours_2","10X Ours_18")]


mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D),color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_matselectmarkers.pdf",sep=""),scale="row",width=10,height=12)

Markers2 <- c("ACTA2","C7","IGF1","PCOLCE","MMP11","ECM1","FOXO1","IL15","CFD","CEBPB","PDGFA")



mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
mammal.combined2$Cells <- Idents(mammal.combined2)
Idents(mammal.combined2) <- mammal.combined2$Cells
mammal.combined2 <- subset(mammal.combined2,idents=c("Stromal fibroblasts"))
Idents(mammal.combined2) <- paste(mammal.combined2$Dataset,mammal.combined2$Cl05,sep="_")

DefaultAssay(mammal.combined2) <- "RNA"


Av <- AverageExpression(mammal.combined2)
Av <- Av$RNA
D <- Av[Markers2,c("10X Ours_16","10X Ours_8","10X Ours_17","10X Ours_7","10X Ours_10","10X Ours_6")]


mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D),color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_matselectmarkers2.pdf",sep=""),scale="row",width=10,height=12)



mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
mammal.combined2$Cells <- Idents(mammal.combined2)
Idents(mammal.combined2) <- mammal.combined2$Dataset
mammal.combined2 <- subset(mammal.combined2,idents=c("10X Ours"))
Idents(mammal.combined2) <- mammal.combined2$Cells
mammal.combined2 <- subset(mammal.combined2,idents=c("Unciliated epithelia","Ciliated","Stromal fibroblasts"))
DefaultAssay(mammal.combined2) <- "RNA"
Idents(mammal.combined2) <- paste(mammal.combined2$Dataset,mammal.combined2$Cl05,sep="_")
PregMark <- c("CST6",
"FOS",
"VLDLR",
"ISG15",
"IFI6",
"MX2",
"C15H11ORF34",
"EIF3M",
"TINAGL1",
"PXDN",
"PRSS22",
"TUBD1",
"TMPRSS2",
"C1R",
"FLVCR2",
"SLC2A1",
"R3HDM1",
"MX1",
"BCAM",
"ISG15",
"IFI6",
"PENK",
"PRSS22",
"MS4A8",
"CLDN4",
"C15H11ORF34",
"MRS2",
"TINAGL1",
"R3HDM1",
"MX1",
"GPT2",
"OAS1Y",
"LRWD1",
"MX2",
"TRIM34",
"IRF9",
"NDRG2",
"SLC2A1",
"TRANK1")
DefaultAssay(mammal.combined2) <- "RNA"

FeaturePlot(mammal.combined2, features=intersect(PregMark,rownames(mammal.combined2)) )
ggsave(filename=paste(saveext,"PregMarker",".pdf",sep=""),width = 100, height = 100, limitsize = FALSE)

#width = 100, height = 100, limitsize = FALSE 

#Glands1 <- intersect(Type3,Type4)
##Glands2 <- Gland
#3Prolif1 <- intersect(Type3B,Type4)
#Lumenal1 <- intersect(Type3D,Type4)
#SOX9P <- intersect(Type3F,Type4)
#LRG5P <- intersect(Type3G,Type4)

#I#dents(mammal.combined2,cells=Glands1) <- "Glads"


#Type3<-WhichCells(mammal.combined2,idents=c(1,5))
#Type3B<-WhichCells(mammal.combined2,idents=c(0))
#Type3C<-WhichCells(mammal.combined2,idents=c(22,13,12,8,4,3))
#Type3D<-WhichCells(mammal.combined2,idents=c(2))
#Type3E<-WhichCells(mammal.combined2,idents=c(22))
#Type3F<-WhichCells(mammal.combined2,idents=c(4,3,19,12))
#oType3G<-WhichCells(mammal.combined2,idents=c(13))

