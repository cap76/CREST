library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library(stringr)
library("readxl")
set.seed(1) 

saveext = "~/Desktop/Data/Endometrial/InVitro/Matteo/MaternalAnalysis/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

PrimedNaiveForm <- read_excel("/Users/christopherpenfold/Downloads/mmc3.xls")
PrimedNaive <- read_excel("/Users/christopherpenfold/Downloads/mmc2.xls")
PrimedNaiveForm <- read_excel("/Users/christopherpenfold/Downloads/mmc3.xls")
AmnionMarkers <- read_excel("/Users/christopherpenfold/Downloads/AmnionMarkers.xlsx")

TomMarkers <- c("PLA2G2A","DIO2","SCARA5","CXCL8","CXCL14","TIMP3","IL15","ECM1","AXL","WNT7A","PTGS1","HEY1","DNAI1","FOXJ1","SPP1","PAEP","DPP4","LIF","AREG","EREG","NEAT1","KCNQ1OT1","WNT5A")

#Load the good samples
GoodSamples <- readRDS(file=paste("~/Desktop/Data/Endometrial/InVitro/Matteo//AnotatedMatteoDataset.rds",sep=""))
Idents(GoodSamples) <- GoodSamples$HarmAno


#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
#Also the assembloids data
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________

raw_counts1<-read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Assembloids/GSE168405_assembloid_counts.csv",sep=",", row.names=1, header=T)
BS<-read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Assembloids/GSE168405_assembloid_metadata.csv",sep=",",header = T, row.names=1)
assemb  <- CreateSeuratObject(counts = raw_counts1, assay = "RNA",min.cells = 0, min.features = 0)
Idents(assemb) <- BS$Cluster
assemb <- NormalizeData(assemb, verbose = FALSE)
assemb <- FindVariableFeatures(assemb, selection.method = "vst", nfeatures = 20000)
assemb$ID1 <- "Assemb"
assemb$ID2 <- BS$Cluster
assemb$ID3 <- BS$SampleName
assemb$ID4 <- BS$Patient
assemb$ID5 <- BS$Treatment
Idents(assemb) <- BS$Cluster
Assemb1 <- FindMarkers(assemb,ident.1 = "SS1", ident.2 = c("SS2","SS3","SS4","SS5"),test="MAST", only.pos=TRUE)
Assemb2 <- FindMarkers(assemb,ident.1 = "SS2", ident.2 = c("SS1","SS3","SS4","SS5"),test="MAST", only.pos=TRUE)
Assemb3 <- FindMarkers(assemb,ident.1 = "SS3", ident.2 = c("SS2","SS1","SS4","SS5"),test="MAST", only.pos=TRUE)
Assemb4 <- FindMarkers(assemb,ident.1 = "SS4", ident.2 = c("SS2","SS3","SS1","SS5"),test="MAST", only.pos=TRUE)
Assemb5 <- FindMarkers(assemb,ident.1 = "SS5", ident.2 = c("SS2","SS3","SS4","SS1"),test="MAST", only.pos=TRUE)



newdata1 <- Assemb1[order(Assemb1$avg_log2FC,decreasing = TRUE),]
newdata2 <- Assemb1[order(Assemb2$avg_log2FC,decreasing = TRUE),]
newdata3 <- Assemb1[order(Assemb3$avg_log2FC,decreasing = TRUE),]
newdata4 <- Assemb1[order(Assemb4$avg_log2FC,decreasing = TRUE),]
newdata5 <- Assemb1[order(Assemb5$avg_log2FC,decreasing = TRUE),]

Stroma <- AddModuleScore(Stroma, features = list( rownames(newdata1)[1:200] ), name = "TomSS1")
Stroma <- AddModuleScore(Stroma, features = list( rownames(newdata2)[1:200] ), name = "TomSS2")
Stroma <- AddModuleScore(Stroma, features = list( rownames(newdata3)[1:200] ), name = "TomSS3")
Stroma <- AddModuleScore(Stroma, features = list( rownames(newdata4)[1:200] ), name = "TomSS4")
Stroma <- AddModuleScore(Stroma, features = list( rownames(newdata5)[1:200] ), name = "TomSS5")


Mat <- as.data.frame(t(rbind(Stroma$TomSS11,Stroma$TomSS21,Stroma$TomSS31,Stroma$TomSS41,Stroma$TomSS51)))
test <- apply(Mat, 1, function(x) which(x == max(x), arr.ind = TRUE)[1])
IDs <- c("Dividing","E2","PreDec","EmDecid","SenDec")
newIDs <- as.character(IDs[test])
Stroma$TomAno <- newIDs
saveRDS(Stroma,file=paste(saveext,"/Stroma_wTomAno.rds",sep=""))


#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
#Load in the placenta data
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
Endo0 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Placenta/Placentaharm.rds")
Idents(Endo0) <- Endo0$ID2
Endo0 <- subset(Endo0,idents=c("dS1","dS2","dS3","Epi1","Epi2"))

#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
#And endometrial dataset
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
Endo2 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo2) <- Endo2$ID3
Endo2 <- subset(Endo2,idents=c("early-secretory","proliferative","late-secretory") )
Idents(Endo2) <- Endo2$ID5
Endo2 <- subset(Endo2,idents=c("Glandular","SOX9","Lumenal","Ciliated","dS","eS") ) #"dS","Fibroblast C7","eS"

#Load in the endometrial reference dataset
Endo3 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo3) <- Endo3$ID3
Endo3 <- subset(Endo3,idents=c("early-secretory","proliferative","late-secretory") )
Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents=c("Glandular","SOX9","Lumenal","Ciliated") ) #"dS","Fibroblast C7","eS"

#subset our dataset by epithelial and stromal cells
GoodSamples0 <- subset(GoodSamples,idents = c("EpithelialOCT4","Epithelial","Ciliated","Stroma1"))
#subset our dataset by epithelial cells
GoodSamples1 <- subset(GoodSamples,idents = c("EpithelialOCT4","Epithelial","Ciliated"))
DefaultAssay(GoodSamples1) <- "RNA"
#SOX <- WhichCells(GoodSamples,expression=SOX9> log(1+1) )
#LRG <- WhichCells(GoodSamples,expression=LGR5> log(1+1) )
#SPP <- WhichCells(GoodSamples,expression=SPP1> log(1+1) )
#SCGB <- WhichCells(GoodSamples,expression=SCGB2A2> log(1+1) )

#Do some DE
Idents(Endo2) <- Endo2$ID5
Idents(Endo2) <- paste(Endo2$ID3,Idents(Endo2),sep="")
#All general markers
ListA1 <- FindMarkers(Endo2,ident.1 = c("late-secretorydS","late-secretoryeS"),ident.2 = c("late-secretoryCiliated","late-secretoryLumenal","late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
ListA2 <- FindMarkers(Endo2,ident.2 = c("late-secretorydS","late-secretoryeS"),ident.1 = c("late-secretoryCiliated","late-secretoryLumenal","late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
ListB1 <- FindMarkers(Endo2,ident.1 = c("proliferativedS","proliferativeeS"),ident.2 = c("proliferativeSOX9","proliferativeCiliated","proliferativeLumenal","proliferativeGlandular"), test.use = "MAST", only.pos = TRUE )
ListB2 <- FindMarkers(Endo2,ident.2 = c("proliferativedS","proliferativeeS"),ident.1 = c("proliferativeSOX9","proliferativeCiliated","proliferativeLumenal","proliferativeGlandular"), test.use = "MAST", only.pos = TRUE )
ListC1 <- FindMarkers(GoodSamples0,ident.1 = c("Stroma1"),ident.2 = c("Ciliated","EpithelialOCT4","Epithelial"), test.use = "MAST", only.pos = TRUE )
ListC2 <- FindMarkers(GoodSamples0,ident.2 = c("Stroma1"),ident.1 = c("Ciliated","EpithelialOCT4","Epithelial"), test.use = "MAST", only.pos = TRUE )

#ListA3 <- FindMarkers(Endo2,ident.1 = c("late-secretorydS","late-secretoryeS"),ident.2 = c("late-secretoryCiliated","late-secretoryLumenal","late-secretoryGlandular"), test.use = "MAST" )
#ListB3 <- FindMarkers(Endo2,ident.1 = c("proliferativedS","proliferativeeS"),ident.2 = c("proliferativeSOX9","proliferativeCiliated","proliferativeLumenal","proliferativeGlandular"), test.use = "MAST")

Idents(GoodSamples0,WhichCells(GoodSamples0,idents=c("Ciliated","EpithelialOCT4"))) <- "Epithelial"
ListC3 <- FindMarkers(GoodSamples0,ident.2 = c("Stroma1"),ident.1 = c("Epithelial"), test.use = "MAST")

AvExp <- AverageExpression(GoodSamples0)
#AvExp <- AverageExpression(mammal.combined3)
AvExp <- as.data.frame(AvExp$RNA)
AvExp$LFC <- NA
AvExp$pval <- NA
AvExp[rownames(ListC3),"LFC"] <- ListC3$avg_log2FC
AvExp[rownames(ListC3),"pval"] <- ListC3$p_val_adj
AvExp$AvExp <- log2((AvExp$Epithelial+AvExp$Stroma1)/2)

p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
genes.to.label <- genes.to.label <- unique(c(intersect(rownames(ListA1),rownames(ListC1)),intersect(rownames(ListA2),rownames(ListC2)) ) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Stroma_v_Epithelia_secMrks.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
genes.to.label <- genes.to.label <- unique(c(intersect(rownames(ListB1),rownames(ListC1)),intersect(rownames(ListB2),rownames(ListC2)) ) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Stroma_v_Epithelia_prolifMrks.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

saveRDS(ListA1,file=paste(saveext,"/ListA1.rds",sep=""))
saveRDS(ListA2,file=paste(saveext,"/ListA2.rds",sep=""))
saveRDS(ListB1,file=paste(saveext,"/ListB1.rds",sep=""))
saveRDS(ListB2,file=paste(saveext,"/ListB2.rds",sep=""))
saveRDS(ListC1,file=paste(saveext,"/ListC1.rds",sep=""))
saveRDS(ListC2,file=paste(saveext,"/ListC2.rds",sep=""))
saveRDS(ListC3,file=paste(saveext,"/ListC3.rds",sep=""))
#mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo2,GoodSamples0), dims = 1:30, anchor.features = 5000)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithelstrom.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_epithel.pdf",sep=""),width = 10, height = 8,p)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,3)) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_3.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,8)) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_8.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,4)) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_4.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,5)) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_5.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,6)) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_6.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,7)) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_7.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split_epithel.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_epithel.pdf",sep=""),width = 10, height = 8,p)


####___________________________________________________________________________
####___________________________________________________________________________
####___________________________________________________________________________
####___________________________________________________________________________
#Now we will concentrate on the epithelial only
####___________________________________________________________________________
####___________________________________________________________________________
####___________________________________________________________________________
####___________________________________________________________________________

Endo4 <- subset(Endo3,idents=c("SOX9","Ciliated","Lumenal","Glandular"))
Idents(Endo4) <- paste(Endo4$IDun3,Idents(Endo4),sep="")
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_test.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

List1 <- FindMarkers(Endo4,ident.2 = "proliferativeLumenal",ident.1 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List2 <- FindMarkers(Endo4,ident.1 = "proliferativeLumenal",ident.2 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List3 <- FindMarkers(Endo4,ident.1 = "proliferativeSOX9",ident.2 = c("proliferativeCiliated","proliferativeGlandular","proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
List4 <- FindMarkers(Endo4,ident.1 = "proliferativeCiliated",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
List5 <- FindMarkers(Endo4,ident.1 = "proliferativeGlandular",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
List6 <- FindMarkers(Endo4,ident.2 = "proliferativeGlandular",ident.1 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )

saveRDS(List1,file=paste(saveext,"/List1.rds",sep=""))
saveRDS(List2,file=paste(saveext,"/List2.rds",sep=""))
saveRDS(List3,file=paste(saveext,"/List3.rds",sep=""))
saveRDS(List4,file=paste(saveext,"/List4.rds",sep=""))
saveRDS(List5,file=paste(saveext,"/List5.rds",sep=""))
saveRDS(List6,file=paste(saveext,"/List6.rds",sep=""))

newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List2$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4 <- List4[order(List4$avg_log2FC,decreasing = TRUE),]
newdata5 <- List5[order(List5$avg_log2FC,decreasing = TRUE),]
newdata6 <- List6[order(List6$avg_log2FC,decreasing = TRUE),]

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo4,GoodSamples1), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_epithel.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_3.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,8)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_8.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_4.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,5)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_5.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,6)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_6.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,7)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_7.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split_epithel.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_epithel.pdf",sep=""),width = 10, height = 8,p)

DefaultAssay(mammal.combined) <- "RNA"
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata1)[1:100] ), name = "Gland")
p1<-FeaturePlot(mammal.combined, features = "Gland1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Gland_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata2)[1:100] ), name = "Luminal")
p1<-FeaturePlot(mammal.combined, features = "Luminal1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Lum_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata3)[1:100] ), name = "SOX9P_")
p1<-FeaturePlot(mammal.combined, features = "SOX9P_1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SOX9_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata4)[1:100] ), name = "Cil")
p1<-FeaturePlot(mammal.combined, features = "Cil1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Cil_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata5)[1:100] ), name = "eGland")
p1<-FeaturePlot(mammal.combined, features = "eGland1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/earlyGland_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata6)[1:100] ), name = "eLuminal")
p1<-FeaturePlot(mammal.combined, features = "eLuminal1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/earlyLumenal_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)

ID6 <- as.character(mammal.combined$ID3)
ID6[which(ID6%in%c("C1","C2","C3","C4","C5","C6","C7","C8"))] <- "Ours"
mammal.combined$ID6 <- ID6

Idents(mammal.combined) <- mammal.combined$ID6
mammal.combinedours <- subset(mammal.combined,idents="Ours")
p1<-FeaturePlot(mammal.combinedours, features = "Gland1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE, dims=c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/Gland_ours_1_4.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "Luminal1", pt.size = 4, reduction = "pca",label = TRUE, repel = TRUE, dims=c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/Lum_ours_1_4.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "SOX9P_1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE, dims=c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/SOX9_ours_1_4.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "Cil1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/Cil_ours_1_4.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "eGland1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE, dims=c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/earlyGland_ours_1_4.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "eGland1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE, dims=c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/earlyLumenal_ours_1_4.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined$stachID <- Idents(mammal.combined)

#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- FindClusters(mammal.combined, resolution = 0.1)
#p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,4)) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_4_clust.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

Idents(mammal.combined) <- mammal.combined$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C6","C7","C8"))
notOurs1 <- subset(mammal.combined,idents=c("proliferative","late-secretory"))
Idents(notOurs1) <- paste(notOurs1$ID3,notOurs1$ID5,sep="_")
p <- FeatureScatter(onlyOurs, feature1 = "Gland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_Gland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "Gland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_Gland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- FindClusters(mammal.combined, resolution = 0.1)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_4_clust.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
       
mammal.combined$Clusters <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$Clusters


Idents(mammal.combined) <- mammal.combined$Dataset
mammal.combined <- subset(mammal.combined,idents=c("Batch1","10X Ours","Batch3","Batch6","Batch7","Batch8"))
#mammal.combined <- D1
Idents(mammal.combined) <- mammal.combined$Clusters
Idents(mammal.combined,WhichCells(mammal.combined,idents=c("6")) ) <- "Ciliated"
Idents(mammal.combined,WhichCells(mammal.combined,idents=c("2")) ) <- "SOX9"
Idents(mammal.combined,WhichCells(mammal.combined,idents=c("1","3","7")) ) <- "Glandular"
Idents(mammal.combined,WhichCells(mammal.combined,idents=c("0","4","5","8")) ) <- "Luminal"
Idents(mammal.combined,colnames(mammal.combined)[which(mammal.combined$HarmAno=="Ciliated")]) <- "Ciliated"

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_4_clustano.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

mammal.combined$EpithelialAno <- Idents(mammal.combined)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID6", label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_4_clustano.pdf",sep=""),width = 40, height = 8,p, limitsize = FALSE)
p<-FeatureScatter(mammal.combinedours, feature1 = "SCGB2A2", feature2 = "SPP1", pt.size = 4) 
ggsave(filename=paste(saveext,"/DimRed/SCGB2A2_v_SPP.pdf",sep=""),width = 10, height = 10,p, limitsize = FALSE)

Idents(mammal.combined) <- mammal.combined$ID6
mammal.combinedours2 <- subset(mammal.combined,idents="Ours")
Idents(mammal.combinedours2) <- mammal.combinedours2$EpithelialAno

DE <- FindMarkers(mammal.combinedours2,ident.1 = "Glandular", ident.2 = c("Luminal"), test = "MAST")
p <- FeatureScatter(onlyOurs, feature1 = "Cil1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_Cil_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "Cil1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_Cil_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(onlyOurs, feature1 = "SOX9P_1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_SOX9_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "SOX9P_1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_SOX9_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(onlyOurs, feature1 = "eGland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_earlyland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "eGland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_earlyGland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(onlyOurs, feature1 = "eGland1", feature2 = "eLuminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_earlyland_vs_eLuminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "eGland1", feature2 = "eLuminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_earlyGland_vs_eLuminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

saveRDS(mammal.combined,file=paste(saveext,"/AnotatedEpithelialCells.rds",sep=""))

#Now DE for our markers
ListX1 <- FindMarkers(mammal.combinedours2,ident.2 = "Luminal",ident.1 = c("Glandular"), test.use = "MAST", only.pos = TRUE )
ListX2 <- FindMarkers(mammal.combinedours2,ident.1 = "Luminal",ident.2 = c("Glandular"), test.use = "MAST", only.pos = TRUE )
ListX3 <- FindMarkers(mammal.combinedours2,ident.1 = "Luminal",ident.2 = c("Glandular"), test.use = "MAST")

ListY1 <- FindMarkers(mammal.combinedours2,ident.2 = c("Luminal","SOX9","Glandular","Luminal"),ident.1 = c("Ciliated"), test.use = "MAST", only.pos = TRUE )
ListY2 <- FindMarkers(mammal.combinedours2,ident.1 = c("Luminal","SOX9","Glandular","Luminal"),ident.2 = c("Ciliated"), test.use = "MAST", only.pos = TRUE )
ListY3 <- FindMarkers(mammal.combinedours2,ident.1 = c("Luminal","SOX9","Glandular","Luminal"),ident.2 = c("Ciliated"), test.use = "MAST")


#List1 <- FindMarkers(Endo4,ident.2 = "proliferativeLumenal",ident.1 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
#List2 <- FindMarkers(Endo4,ident.1 = "proliferativeLumenal",ident.2 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
#List3 <- FindMarkers(Endo4,ident.1 = "proliferativeSOX9",ident.2 = c("proliferativeCiliated","proliferativeGlandular","proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
#List4 <- FindMarkers(Endo4,ident.1 = "proliferativeCiliated",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
#List5 <- FindMarkers(Endo4,ident.1 = "proliferativeGlandular",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
#List6 <- FindMarkers(Endo4,ident.2 = "proliferativeGlandular",ident.1 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )


saveRDS(ListX1,file=paste(saveext,"ListX1.rds",sep=""))
saveRDS(ListX2,file=paste(saveext,"ListX2.rds",sep=""))
saveRDS(ListX3,file=paste(saveext,"ListX3.rds",sep=""))
saveRDS(ListY1,file=paste(saveext,"ListY1.rds",sep=""))
saveRDS(ListY2,file=paste(saveext,"ListY2.rds",sep=""))
saveRDS(ListY3,file=paste(saveext,"ListY3.rds",sep=""))

#
ListCil1 <- FindMarkers(Endo4,ident.1 = "proliferativeCiliated",ident.2 = c("proliferativeGlandular","proliferativeLumenal","proliferativeSOX9"), test.use = "MAST", only.pos = TRUE )
ListCil2 <- FindMarkers(Endo4,ident.2 = "proliferativeCiliated",ident.1 = c("proliferativeGlandular","proliferativeLumenal","proliferativeSOX9"), test.use = "MAST", only.pos = TRUE )
ListCil3 <- FindMarkers(Endo4,ident.1 = "late-secretoryCiliated",ident.2 = c("late-secretoryGlandular","late-secretoryLumenal"), test.use = "MAST", only.pos = TRUE )
ListCil4 <- FindMarkers(Endo4,ident.2 = "late-secretoryCiliated",ident.1 = c("late-secretoryGlandular","late-secretoryLumenal"), test.use = "MAST", only.pos = TRUE )
saveRDS(ListCil1,file=paste(saveext,"ListCil1.rds",sep=""))
saveRDS(ListCil2,file=paste(saveext,"ListCil2.rds",sep=""))
saveRDS(ListCil3,file=paste(saveext,"ListCil3.rds",sep=""))
saveRDS(ListCil4,file=paste(saveext,"ListCil4.rds",sep=""))

#ListZ1 <- FindMarkers(Endo4,ident.2 = "proliferativeLumenal",ident.1 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
#ListZ2 <- FindMarkers(Endo4,ident.1 = "proliferativeLumenal",ident.2 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )

mammal.combinedours2$stashID <- Idents(mammal.combinedours2)
Idents(mammal.combinedours2) <- mammal.combinedours2$stashID
DefaultAssay(mammal.combinedours2) <- "RNA"
#Now we can do some internal comparisons
AvExp <- AverageExpression(mammal.combinedours2)
AvExp <- as.data.frame(AvExp$RNA)
AvExp$LFC <- NA
AvExp$pval <- NA
AvExp[rownames(ListX3),"LFC"] <- ListX3$avg_log2FC
AvExp[rownames(ListX3),"pval"] <- ListX3$p_val_adj
AvExp$AvExp <- log2((AvExp$Luminal+AvExp$Glandular)/2)

p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
genes.to.label <- genes.to.label <- unique(c(intersect(rownames(ListX1),rownames(List1)),intersect(rownames(ListX2),rownames(List2)) ) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_LuminalEpithelia_vs_Glandular_prolifmarkers.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)


p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
genes.to.label <- genes.to.label <- unique(c(intersect(rownames(ListX1),rownames(List5)),intersect(rownames(ListX2),rownames(List6)) ) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_LuminalEpithelia_vs_Glandular_secmarkers.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

Idents(mammal.combinedours2) <- mammal.combinedours2$stashID
#Idents(mammal.combinedours2) <- mammal.combinedours2$stashID
Idents(mammal.combinedours2,cells=WhichCells(mammal.combinedours2,idents=c("Luminal","SOX9","Glandular"))) <- "NotCil"
#Now we can do some internal comparisons
AvExp <- AverageExpression(mammal.combinedours2)
AvExp <- as.data.frame(AvExp$RNA)
AvExp$LFC <- NA
AvExp$pval <- NA
AvExp[rownames(ListY3),"LFC"] <- ListY3$avg_log2FC
AvExp[rownames(ListY3),"pval"] <- ListY3$p_val_adj
AvExp$AvExp <- log2((AvExp$NotCil+AvExp$Ciliated)/2)

p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
genes.to.label <- genes.to.label <- unique(c(intersect(rownames(ListY1),rownames(ListCil3)),intersect(rownames(ListY2),rownames(ListCil4)) ) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Cil_vs_nonCil_secmarkers.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)


p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
genes.to.label <- genes.to.label <- unique(c(intersect(rownames(ListY1),rownames(ListCil1)),intersect(rownames(ListY2),rownames(ListCil2)) ) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Cil_vs_nonCil_prolifmarkers.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

###_________________________________________________________________________________________________________
###_________________________________________________________________________________________________________
###_________________________________________________________________________________________________________
###_________________________________________________________________________________________________________
###Now for the stroma
###_________________________________________________________________________________________________________
###_________________________________________________________________________________________________________
###_________________________________________________________________________________________________________
###_________________________________________________________________________________________________________


GoodSamples <- readRDS(file=paste("~/Desktop/Data/Endometrial/InVitro/Matteo//AnotatedMatteoDataset.rds",sep=""))
Idents(GoodSamples) <- GoodSamples$HarmAno
GoodSamples2 <- subset(GoodSamples,idents=c("Stroma1"))

#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
#Also the assembloids data
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________

raw_counts1<-read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Assembloids/GSE168405_assembloid_counts.csv",sep=",", row.names=1, header=T)
BS<-read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Assembloids/GSE168405_assembloid_metadata.csv",sep=",",header = T, row.names=1)
assemb  <- CreateSeuratObject(counts = raw_counts1, assay = "RNA",min.cells = 0, min.features = 0)
Idents(assemb) <- BS$Cluster
assemb <- NormalizeData(assemb, verbose = FALSE)
assemb <- FindVariableFeatures(assemb, selection.method = "vst", nfeatures = 20000)
assemb$ID1 <- "Assemb"
assemb$ID2 <- BS$Cluster
assemb$ID3 <- BS$SampleName
assemb$ID4 <- BS$Patient
assemb$ID5 <- BS$Treatment
Idents(assemb) <- BS$Cluster
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
#Load in the placenta data
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
Endo0 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Placenta/Placentaharm.rds")
Idents(Endo0) <- Endo0$ID2
Endo0 <- subset(Endo0,idents=c("dS1","dS2","dS3","Epi1","Epi2"))

#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
#And endometrial dataset again, becuase why not.
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
Endo3 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo3) <- Endo3$ID3
Endo3 <- subset(Endo3,idents=c("early-secretory","proliferative","late-secretory") )
Idents(Endo3) <- Endo3$ID5
Endo4 <- subset(Endo3,idents=c("Fibroblast C7","dS","eS"))
#D_s <- subset(D,idents=c("Stromal fibroblasts"))

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo4,GoodSamples2), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_stroma.pdf",sep=""),width = 10, height = 8,p)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_3.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_4.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,5)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_5.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,6)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_6.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,7)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_7.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_stroma.pdf",sep=""),width = 10, height = 8,p)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split_stroma.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID6", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_dataset.pdf",sep=""),width = 40, height = 8,p, limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID6", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_7_dataset.pdf",sep=""),width = 40, height = 8,p, limitsize = FALSE)

Idents(Endo4) <- paste(Endo4$ID3,Idents(Endo4),sep="")

List1 <- FindMarkers(Endo4,ident.2 = "proliferativeeS",ident.1 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List2 <- FindMarkers(Endo4,ident.1 = "proliferativeeS",ident.2 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List3 <- FindMarkers(Endo4,ident.1 = "proliferativeFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
#List4 <- FindMarkers(Endo4,ident.1 = "late-secretoryFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
#List5 <- FindMarkers(Endo4,ident.2 = "proliferativeFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
#List6 <- FindMarkers(Endo4,ident.2 = "late-secretoryFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )

saveRDS(List2,file=paste(saveext,"ProlifStroma_vs_SecdStroma.rds",sep=""))
saveRDS(List1,file=paste(saveext,"SecdStroma_vs_ProlifStroma.rds",sep=""))
saveRDS(List3,file=paste(saveext,"ProlfiFib_vs_ProlifStroma.rds",sep=""))
saveRDS(List4,file=paste(saveext,"ProlifStroma_vs_ProlfiFib.rds",sep=""))
saveRDS(List5,file=paste(saveext,"SecFib_vs_SecStroma.rds",sep=""))
saveRDS(List6,file=paste(saveext,"SecStroma_vs_SecFib.rds",sep=""))

newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List2$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4 <- List4[order(List4$avg_log2FC,decreasing = TRUE),]
newdata5 <- List5[order(List5$avg_log2FC,decreasing = TRUE),]
newdata6 <- List6[order(List6$avg_log2FC,decreasing = TRUE),]

newdata13 <- List13[order(List13$avg_log2FC,decreasing = TRUE),]
newdata14 <- List14[order(List14$avg_log2FC,decreasing = TRUE),]

Idents(mammal.combined) <- mammal.combined$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
notOurs <- subset(mammal.combined,idents=c("proliferative","late-secretory"))

DefaultAssay(mammal.combined) <- "RNA"
Idents(mammal.combined) <- mammal.combined$ID3
mammal.combined3 <- subset(mammal.combined,idents=c("C1","C2","C3","C6","C7","C8"))
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata1)[1:200] ), name = "dS")
p1<-FeaturePlot(mammal.combined3, features = "dS1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_dS_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata2)[1:200] ), name = "prolifS")
p1<-FeaturePlot(mammal.combined3, features = "prolifS1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_prolifS_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
p1<-FeaturePlot(mammal.combined3, features = "prolifC71", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_prolifFibC7_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata4)[1:100] ), name = "lateC7")
p1<-FeaturePlot(mammal.combined3, features = "lateC71", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_lateFibC7_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata13)[1:200] ), name = "dS2")
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata14)[1:200] ), name = "prolifS2")

p <- FeatureScatter(mammal.combined3, feature1 = "dS1", feature2 = "prolifS1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_Stroma.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(mammal.combined3, feature1 = "dS21", feature2 = "prolifS21", )
ggsave(filename=paste(saveext,"/DimRed/Ours_Stroma2.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(mammal.combined, feature1 = "dS21", feature2 = "prolifS21", )
ggsave(filename=paste(saveext,"/DimRed/All_Stroma2.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

StromaAno <- as.character(mammal.combined3$ID1)
StromaAno[1:length(StromaAno)] <- "dS"
StromaAno[which(mammal.combined3$prolifS1>mammal.combined3$dS1)] <- "pS"

Idents(mammal.combined3) <- StromaAno

p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_newAno.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)
mammal.combined3$StromaID <- StromaAno

saveRDS(mammal.combined3,file=paste(saveext,"/AnotatedStromaCells.rds",sep=""))
Idents(mammal.combined3) <- mammal.combined3$Genotype
p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_GT.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)

mammal.combined5 <- mammal.combined #subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))
mammal.combined5 <- AddModuleScore(mammal.combined5, features = list( rownames(newdata1)[1:100] ), name = "dS")
mammal.combined5 <- AddModuleScore(mammal.combined5, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
mammal.combined5 <- AddModuleScore(mammal.combined5, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
mammal.combined5 <- AddModuleScore(mammal.combined5, features = list( rownames(newdata4)[1:100] ), name = "lateC7")


mammal.combined4 <- subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata1)[1:100] ), name = "dS")
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata4)[1:100] ), name = "lateC7")


#mammal.combined3 <- subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata1)[1:100] ), name = "dS")
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata4)[1:100] ), name = "lateC7")

#mammal.combined3 <- subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata1)[1:100] ), name = "dS")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata4)[1:100] ), name = "lateC7")


Idents(mammal.combined) <- mammal.combined3$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C7","C6","C8"))
uID <- as.character( Idents(onlyOurs) )
uID[1:length(uID)] <- "dS"
uID[which(onlyOurs$prolifS1 > onlyOurs$dS1)] <- "pS"
Idents(onlyOurs) <- uID
DefaultAssay(onlyOurs) <- "RNA"

uID <- as.character( Idents(mammal.combined) )
uID[1:length(uID)] <- "dS"
uID[which(mammal.combined$prolifS1 > mammal.combined$dS1)] <- "pS"
uID[which(mammal.combined$ID3%in%c("proliferative","early-secretory","late-secretory"))] <- "Reference"
Idents(mammal.combined) <- uID

DefaultAssay(mammal.combined) <- "RNA"

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_samescale.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)

mammal.combined$Delta <- mammal.combined$dS1 - mammal.combined$prolifS1
twosplit <- as.character(mammal.combined$ID3)
twosplit[which(twosplit%in%c("C1","C2","C3","C6","C7","C8"))] <- "ours"
twosplit[which(twosplit%in%c("proliferative","early-secretory","late-secretory"))] <- "Ref"
mammal.combined$twosplit <- twosplit
Idents(mammal.combined) <- mammal.combined$Delta

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "Delta", split.by = "twosplit", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdS.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "Delta", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdScb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

Idents(mammal.combined) <- mammal.combined$twosplit
mammal.combinedsss <- subset(mammal.combined,idents=c("ours"))
p<-FeaturePlot(mammal.combinedsss,  reduction = "pca", features = "Delta", split.by = "twosplit", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdS2.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combinedsss,  reduction = "pca", features = "Delta", cols =  c("blue", "red"), pt.size = 2 ) # + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdScb2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

uID <- as.character( Idents(mammal.combinedsss) )
uID[1:length(uID)] <- "dS"
uID[which(mammal.combinedsss$prolifS1 > mammal.combinedsss$dS1)] <- "pS"
Idents(mammal.combinedsss) <- uID
DefaultAssay(mammal.combinedsss) <- "RNA"

mammal.combinedsss$ID4 <- NULL
p<-VlnPlot(mammal.combinedsss,  features = rownames(newdata1)[1:40], cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_decidualised.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p)


p<-VlnPlot(mammal.combinedsss,  features = rownames(newdata2)[1:40], cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_notdecidualised.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p)


p <- FeatureScatter(mammal.combined4, feature1 = "dS1", feature2 = "prolifS1", )
ggsave(filename=paste(saveext,"/DimRed/Ref_Stroma.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

write.table(mammal.combined3$dS1,file=paste(saveext,"/Ours_Stroma_Dec.csv",sep=""))
write.table(mammal.combined3$prolifS1,file=paste(saveext,"/Ours_Stroma_Prolif.csv",sep=""))
write.table(colnames(mammal.combined3),file=paste(saveext,"/StromaID.csv",sep=""))
write.table(mammal.combined3$ID3,file=paste(saveext,"/StromaGenome.csv",sep=""))

DEOursdSpS1 <- FindMarkers(mammal.combined3,ident.1 = "dS",ident.2 = c("pS"), test.use = "MAST",only.pos=TRUE)
DEOursdSpS2 <- FindMarkers(mammal.combined3,ident.2 = "dS",ident.1 = c("pS"), test.use = "MAST",only.pos=TRUE)
DEOursdSpS3 <- FindMarkers(mammal.combined3,ident.1 = "dS",ident.2 = c("pS"), test.use = "MAST")

AvExp <- AverageExpression(mammal.combined3)
AvExp <- as.data.frame(AvExp$RNA)
AvExp$LFC <- NA
AvExp$pval <- NA
AvExp[rownames(DEOursdSpS3),"LFC"] <- DEOursdSpS3$avg_log2FC
AvExp[rownames(DEOursdSpS3),"pval"] <- DEOursdSpS3$p_val_adj
AvExp$AvExp <- log2((AvExp$pS+AvExp$dS)/2)


p1 <- ggplot(AvExp, aes(x=AvExp, y=LFC)) + geom_point() + theme_classic() + geom_hline(aes(yintercept = 0)) + xlim(c(-10,10))
genes.to.label <- genes.to.label <- unique(c( intersect(rownames(List1),rownames(DEOursdSpS1)), intersect(rownames(List2),rownames(DEOursdSpS2)) ) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_dStroma_v_pStroma.pdf",sep=""),width = 20, height = 20, plot = p1, useDingbats=FALSE)

#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
#Now do merge for Vln plots of key markers
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________

Stroma <- readRDS(file=paste(saveext,"/AnotatedStromaCells.rds",sep=""))
Idents(Stroma,WhichCells(Stroma,idents="Impl")) <- "dS"
Idents(Stroma,WhichCells(Stroma,idents="Decid")) <- "dS"
Idents(Stroma,WhichCells(Stroma,idents="Prolif")) <- "pS"
DefaultAssay(Stroma) <- "RNA"
VlnPlot(Stroma, features =c("DIO2","PLA2G2A"), split.by = "ID3")
ggsave(filename=paste(saveext,"VlnPLAG2A_DIO2.pdf",sep=""),width = 40, height = 20, plot = p1, useDingbats=FALSE)

epith <- readRDS(file=paste(saveext,"/AnotatedEpithelialCells.rds",sep=""))
onlyOurs <- subset(epith,idents=c("Ours"))
Idents(onlyOurs) <- onlyOurs$EpithelialAno

DefaultAssay(onlyOurs) <- "RNA"
AllData3 <- merge(onlyOurs,y=c(Stroma,Endo4,Endo0),project = "Merged")

Idents(AllData3,cells=WhichCells(AllData3,idents=c("Luminal"))) <- "11) Luminal"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeLumenal"))) <- "12) proliferativeLuminal"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryLumenal"))) <- "13) early-secretoryLuminal"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("Glandular"))) <- "21) Glandular"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeGlandular"))) <- "22) proliferativeGlandular"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryGlandular"))) <- "23) early-secretoryGlandular"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("late-secretoryGlandular"))) <- "24) late-secretoryGlandular"

Idents(AllData3,cells=WhichCells(AllData3,idents=c("Ciliated"))) <- "31) Ciliated"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeCiliated"))) <- "32) proliferativeCiliated"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryCiliated"))) <- "33) early-secretoryCiliated"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("late-secretoryCiliated"))) <- "34) late-secretoryCiliated"

Idents(AllData3,cells=WhichCells(AllData3,idents=c("Epi1"))) <- "41) Epi1"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("Epi2"))) <- "42) Epi2"

#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS1"))) <- "51) EpS1"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS2"))) <- "52) EpS2"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS3"))) <- "53) EpS3"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS4"))) <- "54) EpS4"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS5"))) <- "55) EpS5"


gL <- c("C20orf85","MORN2","OMG","AGR3","PAEP","KRT17","CP","EDN1","G0S2","SAMD4A")

AllData4 <- subset(AllData3,idents=c("11) Luminal","21) Glandular","31) Ciliated","12) proliferativeLuminal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLuminal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("11) Luminal","21) Glandular","31) Ciliated","12) proliferativeLuminal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLuminal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"))


Dats <- subset(AllData4,idents=c("21) Glandular","11) Luminal","31) Ciliated"))
Dats$ID0 <- Idents(Dats)
Idents(Dats) <- Dats$Dataset
Dats <- subset(Dats,idents="10X Ours")
Idents(Dats) <- Dats$ID0

List1B <- FindMarkers(AllData4,ident.1 = "21) Glandular", ident.2 = "11) Lumenal", test.use = "MAST", only.pos = TRUE) #, min.cells.group = 1, 
                      #min.cells.feature = 1,
                      #min.pct = 0,
                      #logfc.threshold = 0)
List2Bs <- FindMarkers(AllData4,ident.1 = "11) Lumenal", ident.2 = c("21) Glandular"), test.use = "MAST", only.pos = TRUE )#, min.cells.group = 1, 
                      #min.cells.feature = 1,
                      #min.pct = 0,
                      #logfc.threshold = 0)
List4Bs <- FindMarkers(AllData4,ident.1 = "31) Ciliated",ident.2 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE) #, min.cells.group = 1, 
                      #min.cells.feature = 1,
                      #min.pct = 0,
                      #logfc.threshold = 0)

List4Bo <- FindMarkers(AllData4,ident.2 = "31) Ciliated",ident.1 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE)#, min.cells.group = 1, 
                       #min.cells.feature = 1,
                       #min.pct = 0,
                       #logfc.threshold = 0 )


#List1B <- FindMarkers(Dats,ident.1 = "21) Glandular", ident.2 = "11) Lumenal", test.use = "MAST", only.pos = TRUE )
#List2B <- FindMarkers(Dats,ident.1 = "11) Lumenal", ident.2 = c("21) Glandular"), test.use = "MAST", only.pos = TRUE )
#List4B <- FindMarkers(Dats,ident.1 = "31) Ciliated",ident.2 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE )

#List4Bo <- FindMarkers(Dats,ident.2 = "31) Ciliated",ident.1 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE )

#C20orf85,MORN2,OMG,AGR3
#PAEP, KRT17,CP,EDN1,G0S2,SAMD4A

newdata1B <- List1B[order(List1B$avg_log2FC,decreasing = TRUE),]
newdata2B <- List2Bs[order(List2Bs$avg_log2FC,decreasing = TRUE),]
newdata4Bs <- List4Bs[order(List4Bs$avg_log2FC,decreasing = TRUE),]

cl1 <- intersect(rownames(List1B),rownames(List1))
cl2 <- intersect(rownames(List2Bs),rownames(List2))
cl4 <- intersect(rownames(List4Bs),rownames(List4))
cl4o <- intersect(rownames(List4Bo),rownames(List4o))

cl5 <- intersect(rownames(List1B),rownames(List5))
cl6 <- intersect(rownames(List2Bs),rownames(List6))

newdata1t <- List1[cl1,]
newdata1t <- newdata1t[order(newdata1t$avg_log2FC,decreasing = TRUE),]
newdata2t <- List2[cl2,]
newdata2t <- newdata2t[order(newdata2t$avg_log2FC,decreasing = TRUE),]
#newdata3t <- List3[cl2,]
#newdata3t <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4t <- List4[cl4,]
newdata4t <- newdata4t[order(newdata4t$avg_log2FC,decreasing = TRUE),]
newdata5t <- List5[cl5,]
newdata5t <- newdata5t[order(newdata5t$avg_log2FC,decreasing = TRUE),]
newdata6t <- List6[cl6,]
newdata6t <- newdata6t[order(newdata6t$avg_log2FC,decreasing = TRUE),]


newdata1t2 <- List1B[cl1,]
newdata1t2 <- newdata1t2[order(newdata1t2$avg_log2FC,decreasing = TRUE),]
newdata2t2 <- List2Bs[cl2,]
newdata2t2 <- newdata2t2[order(newdata2t2$avg_log2FC,decreasing = TRUE),]
newdata4t2 <- List4Bs[cl4,]
newdata4t2 <- newdata4t2[order(newdata4t2$avg_log2FC,decreasing = TRUE),]
newdata5t2 <- List1B[cl5,]
newdata5t2 <- newdata5t2[order(newdata5t2$avg_log2FC,decreasing = TRUE),]
newdata6t2 <- List2Bs[cl6,]
newdata6t2 <- newdata6t2[order(newdata6t2$avg_log2FC,decreasing = TRUE),]



DefaultAssay(AllData4) <- "RNA"
p1<-VlnPlot(AllData4,  idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata1t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata2t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata4t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata5t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set5.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata6t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set6.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData4,  idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata1t2)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset1_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata2t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset2_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata4t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set4_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata5t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set5_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata6t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set6_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

Dats <- subset(AllData4,idents=c("21) Glandular","11) Lumenal","31) Ciliated"))
DefaultAssay(Dats) <- "RNA"


#List1Ca <- FindMarkers(Dats,ident.2 = "21) Glandular", ident.1 = "11) Lumenal", test.use = "MAST", only.pos = TRUE)
List1C <- FindMarkers(AllData4,ident.1 = "21) Glandular", ident.2 = "11) Lumenal", test.use = "MAST", only.pos = FALSE)#, min.cells.group = 1, 
  #                    min.cells.feature = 1,
  #                    min.pct = 0,
  #                    logfc.threshold = 0)
List4C <- FindMarkers(AllData4,ident.1 = c("31) Ciliated"),ident.2 = c("11) Lumenal","21) Glandular"), test.use = "MAST", only.pos = FALSE)#, min.cells.group = 1, 
   #                   min.cells.feature = 1,
  #                    min.pct = 0,
  #                    logfc.threshold = 0)



#AllData5 <- AllData4
#Idents(AllData5)
Av <- AverageExpression(AllData4)
Cl1 <- List1C
Ae5 <- as.data.frame(Av$RNA)
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_log2FC
pospos1 <- which( abs(Ae5$FC1)>log(1.1) & Ae5$Pval1>-log2(0.01))
Ae5$ColInd <- 0

Ae5[ which( (Ae5$FC1)>log(1.1) & Ae5$Pval1> -log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.1) & Ae5$Pval1> -log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
Ae5$AvExp <- log(0.5*(Ae5$'21) Glandular' + Ae5$'11) Lumenal') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))

genes.to.label <- genes.to.label <- unique(c(cl1,cl2,cl5,cl6) )

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Gland_vs_Lumenal_newnew.pdf",sep=""),width = 20, height = 10, plot = p1, useDingbats=FALSE)
dev.off()


AllData5 <- AllData4
Idents(AllData5,cells=WhichCells(AllData5,idents=c("21) Glandular"))) <- "11) Lumenal"



#AllData5 <- AllData4
#Idents(AllData5)
Av <- AverageExpression(AllData5)
Cl1 <- List4C
Ae5 <- as.data.frame(Av$RNA)
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_log2FC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01))
Ae5$ColInd <- 0

Ae5[ which( (Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
Ae5$AvExp <- log(0.5*(Ae5$'31) Ciliated' + Ae5$'11) Lumenal') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))

genes.to.label <- genes.to.label <- unique(c(cl4,cl4o) )

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Ciliated_vs_Other_newnew.pdf",sep=""),width = 20, height = 10, plot = p1, useDingbats=FALSE)
dev.off()


DefaultAssay(AllData4) <- "RNA"
p1<-FeatureScatter(AllData4,feature1 = "KRT17",feature2 = "SCGB2A1", cells = WhichCells(AllData4,idents = c("21) Glandular","11) Lumenal")))
ggsave(filename=paste(saveext,"MA_KRT17_SCGB2A1_CiliatedandGland_vs_Luminal.pdf",sep=""),width = 10, height = 10, plot = p1, useDingbats=FALSE)


Idents(AllData4,cells=WhichCells(AllData4,idents=c("42) Epi2"))) <- "41) Epi1"


#
List_imp1 <- FindMarkers(AllData4,ident.1 = c("41) Epi1"),ident.2 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List_imp2 <- FindMarkers(AllData4,ident.2 = c("41) Epi1"),ident.1 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )

List_imp1 <- List_imp1[order(List_imp1$avg_log2FC,decreasing = TRUE),]
List_imp2 <- List_imp2[order(List_imp2$avg_log2FC,decreasing = TRUE),]

List_imp3 <- FindMarkers(AllData4,ident.1 = c("21) Glandular"),ident.2 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List_imp4 <- FindMarkers(AllData4,ident.2 = c("21) Glandular"),ident.1 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )

List_imp3 <- List_imp3[order(List_imp3$avg_log2FC,decreasing = TRUE),]
List_imp4 <- List_imp4[order(List_imp4$avg_log2FC,decreasing = TRUE),]

cl1 <- intersect(rownames(List_imp1),rownames(List_imp3))
cl2 <- intersect(rownames(List_imp2),rownames(List_imp4))

newdata1t <- List_imp1[cl1,]
newdata1t <- newdata1t[order(newdata1t$avg_log2FC,decreasing = TRUE),]
newdata2t <- List_imp2[cl2,]
newdata2t <- newdata2t[order(newdata2t$avg_log2FC,decreasing = TRUE),]

newdata3t <- List_imp3[cl1,]
newdata3t <- newdata3t[order(newdata3t$avg_log2FC,decreasing = TRUE),]
newdata4t <- List_imp4[cl2,]
newdata4t <- newdata4t[order(newdata4t$avg_log2FC,decreasing = TRUE),]

write.table(as.data.frame(newdata1t),file=paste(saveext,"GlandImplantationDE1.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata2t),file=paste(saveext,"GlandImplantationDE2.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata3t),file=paste(saveext,"GlandImplantationDE3.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata4t),file=paste(saveext,"GlandImplantationDE4.rds",sep=""),sep=",",quote = FALSE)

AllData5 <- AllData4
AllData5 <- subset(AllData5,idents=c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"))
Idents(AllData5) <- factor(x = Idents(AllData5), levels = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"))


DefaultAssay(AllData5) <- "RNA"
AllData5$ID4 <- NULL
AllData5$ID1 <- NULL
AllData5$ID2 <- NULL
AllData5$ID3 <- NULL


#______________________________________________________________________________
#______________________________________________________________________________
#Now do dotplots
#______________________________________________________________________________
#______________________________________________________________________________

Endo4 <- subset(Endo3,idents=c("SOX9","Ciliated","Lumenal","Glandular","eS","dS","Fibroblast C7"))
Idents(Endo4) <- paste(Endo4$ID3,Idents(Endo4),sep="")


Stroma <- readRDS(file=paste(saveext,"/AnotatedStromaCells.rds",sep=""))
Idents(Stroma,WhichCells(Stroma,idents="Impl")) <- "dS"
Idents(Stroma,WhichCells(Stroma,idents="Decid")) <- "dS"
Idents(Stroma,WhichCells(Stroma,idents="Prolif")) <- "pS"
DefaultAssay(Stroma) <- "RNA"

epith <- readRDS(file=paste(saveext,"/AnotatedEpithelialCells.rds",sep=""))
onlyOurs <- subset(epith,idents=c("Ours"))
Idents(onlyOurs) <- onlyOurs$EpithelialAno

DefaultAssay(onlyOurs) <- "RNA"
AllData3 <- merge(onlyOurs,y=c(Stroma,Endo4,Endo0),project = "Merged")
Idents(AllData3,WhichCells(AllData3,idents=c("dS1","dS2","dS3"))) <- "pdS"

list1 <- c("proliferativeFibroblast C7","early-secretoryFibroblast C7","late-secretoryFibroblast C7",
"proliferativeeS","early-secretoryeS","late-secretoryeS","pS",
"proliferativedS","early-secretorydS","late-secretorydS","dS","pdS",
"proliferativeSOX9","early-secretorySOX9","late-secretorySOX9","SOX9",
"proliferativeLumenal","early-secretoryLumenal","late-secretoryLumenal","Luminal",
"proliferativeGlandular","early-secretoryGlandular","late-secretoryGlandular","Glandular",
"proliferativeCiliated","early-secretoryCiliated","late-secretoryCiliated","Ciliated","Epi1","Epi2")

#Idents(AllData3,cells=WhichCells(AllData3,idents=c("Luminal"))) <- "11) Luminal"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeLumenal"))) <- "12) proliferativeLuminal"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryLumenal"))) <- "13) early-secretoryLuminal"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("Glandular"))) <- "21) Glandular"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeGlandular"))) <- "22) proliferativeGlandular"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryGlandular"))) <- "23) early-secretoryGlandular"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("late-secretoryGlandular"))) <- "24) late-secretoryGlandular"

#Idents(AllData3,cells=WhichCells(AllData3,idents=c("Ciliated"))) <- "31) Ciliated"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeCiliated"))) <- "32) proliferativeCiliated"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryCiliated"))) <- "33) early-secretoryCiliated"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("late-secretoryCiliated"))) <- "34) late-secretoryCiliated"

#Idents(AllData3,cells=WhichCells(AllData3,idents=c("Epi1"))) <- "41) Epi1"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("Epi2"))) <- "42) Epi2"

#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS1"))) <- "51) EpS1"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS2"))) <- "52) EpS2"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS3"))) <- "53) EpS3"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS4"))) <- "54) EpS4"
#Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS5"))) <- "55) EpS5"
gL <- c("C20orf85","MORN2","OMG","AGR3","PAEP","KRT17","CP","EDN1","G0S2","SAMD4A")

AllData4 <- subset(AllData3,idents=list1)
#  "11) Luminal","21) Glandular","31) Ciliated","12) proliferativeLuminal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLuminal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = list1)
#                             c("11) Luminal","21) Glandular","31) Ciliated","12) proliferativeLuminal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLuminal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"))




#list1 <- c("proliferative Ciliated","early-secretory Ciliated","late-secretory Ciliated",
#           "proliferative Glandular","Ours Ciliated",
#           "early-secretory Glandular","late-secretory Glandular","Ours Glandular","proliferative Lumenal","early-secretory Lumenal","late-secretory Lumenal",
#           "Ours Lumenal","proliferative SOX9","early-secretory SOX9","late-secretory SOX9","Ours SOX9","Ours LRG5","proliferative eS",
#           "early-secretory eS","late-secretory eS","Ours Prolf Strom","proliferative dS","early-secretory dS","late-secretory dS","Ours Stroma")


Markers <- c(
  "SOX9",
  "LGR5",
  "PGR",
  "ESR1",
  "MMP7",
  "CPM",
  "MKI67",
  "HMGB2",
  "PLAU",
  "IL32",
  "TNF",
  "WNT7A",
  "FOXJ1",
  "PIFO",
  "ABCG1",
  "SCGB2A2",
  "C2CD4A",
  "SLC18A2",
  "PAEP",
  "CXCL14",
  "SPP1",
  "PAX2",
  "VTCN1",
  "SLC26A7",
  "MSLN",
  "ACTA2",
  "PCOLCE",
  "MMP11",
  "ECM1",
  "FOXO1")

DefaultAssay(AllData4) <- "RNA"
p<-DotPlot(AllData4, features = Markers, cols = c("lightgrey", "red"), assay = "RNA",  dot.scale =4,scale =FALSE )
ggsave(filename=paste(saveext,"/DimRed/EpitheliaDotPlot.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)
p<-DotPlot(AllData4, features = Markers, cols = c("lightgrey", "red"), assay = "RNA",  scale =TRUE  )
ggsave(filename=paste(saveext,"/DimRed/EpitheliaDotPlotScale.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)

###Idents(AllData3) <- factor(Idents(AllData3),levels=list1)

DefaultAssay(AllData3) <- "RNA"
p<-DotPlot(AllData3, features = Markers, cols = c("lightgrey", "red"), assay = "RNA",  dot.scale =4,scale =FALSE )
ggsave(filename=paste(saveext,"/DimRed/EpitheliaDotPlot.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)
p<-DotPlot(AllData3, features = Markers, cols = c("lightgrey", "red"), assay = "RNA",  scale =TRUE  )
ggsave(filename=paste(saveext,"/DimRed/EpitheliaDotPlotScale.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)


p1<-VlnPlot(AllData4, features = gL, pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllIEptith.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)


#p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata1t)[1:40], pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

#p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata2t)[1:40], pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

#p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata3t)[1:40], pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

#p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata4t)[1:40], pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

########_______________________________________________________________________
########_______________________________________________________________________
#Now some volcano and other plots
########_______________________________________________________________________
########_______________________________________________________________________
Stroma <- readRDS(file=paste(saveext,"/AnotatedStromaCells.rds",sep=""))
Idents(Stroma,WhichCells(Stroma,idents="Impl")) <- "dS"
Idents(Stroma,WhichCells(Stroma,idents="Decid")) <- "dS"
Idents(Stroma,WhichCells(Stroma,idents="Prolif")) <- "pS"
DefaultAssay(Stroma) <- "RNA"

Idents(Endo4) <- paste(Endo4$ID3,Endo4$ID5,sep="")
AllData3 <- merge(Stroma,y=c(Endo4,Endo0,assemb),project = "Merged")
AllData3 <- JoinLayers(AllData3)
AllData4 <- subset(AllData3,idents=c("proliferativeeS","early-secretoryeS","late-secretoryeS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","pS","dS"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("proliferativeeS","early-secretoryeS","late-secretoryeS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","pS","dS"))

List_imp1 <- FindMarkers(AllData4,ident.1 = c("dS1","dS2","dS3"),ident.2 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List_imp2 <- FindMarkers(AllData4,ident.2 = c("dS1","dS2","dS3"),ident.1 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )

saveRDS(List_imp1,file=paste(saveext,"/Stroma_implantation_vs_Sectretory.rds",sep=""))
saveRDS(List_imp2,file=paste(saveext,"/Stroma_Sectretory_vs_implantation.rds",sep=""))

List_imp3 <- FindMarkers(AllData4,ident.1 = c("dS1","dS2","dS3"),ident.2 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE )
List_imp4 <- FindMarkers(AllData4,ident.2 = c("dS1","dS2","dS3"),ident.1 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE )

saveRDS(List_imp3,file=paste(saveext,"/Stroma_implantation_vs_prolif.rds",sep=""))
saveRDS(List_imp4,file=paste(saveext,"/Stroma_prolf_vs_implantation.rds",sep=""))

List_imp1 <- List_imp1[order(List_imp1$avg_log2FC,decreasing = TRUE),]
List_imp2 <- List_imp2[order(List_imp2$avg_log2FC,decreasing = TRUE),]


List_imp5 <- List_imp3[order(List_imp3$avg_log2FC,decreasing = TRUE),]
List_imp6 <- List_imp4[order(List_imp4$avg_log2FC,decreasing = TRUE),]


List_imp5 <- FindMarkers(AllData4,ident.1 = c("dS"),ident.2 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List_imp6 <- FindMarkers(AllData4,ident.2 = c("dS"),ident.1 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE)

List_imp7 <- FindMarkers(AllData4,ident.1 = c("dS"),ident.2 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE )
List_imp8 <- FindMarkers(AllData4,ident.2 = c("dS"),ident.1 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE)

List_imp5 <- List_imp5[order(List_imp5$avg_log2FC,decreasing = TRUE),]
List_imp6 <- List_imp6[order(List_imp6$avg_log2FC,decreasing = TRUE),]

List_imp7 <- List_imp7[order(List_imp7$avg_log2FC,decreasing = TRUE),]
List_imp8 <- List_imp8[order(List_imp8$avg_log2FC,decreasing = TRUE),]

List_imp9 <- FindMarkers(AllData4,ident.1 = c("late-secretorydS","late-secretoryeS"),ident.2 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE )
List_imp10 <- FindMarkers(AllData4,ident.2 = c("late-secretorydS","late-secretoryeS"),ident.1 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE )

List_imp9 <- List_imp9[order(List_imp9$avg_log2FC,decreasing = TRUE),]
List_imp10 <- List_imp10[order(List_imp10$avg_log2FC,decreasing = TRUE),]


saveRDS(List_imp5,file=paste(saveext,"/Stroma_Ourimplantation_vs_Sectretory.rds",sep=""))
saveRDS(List_imp6,file=paste(saveext,"/Stroma_Sectretory_vs_Ourimplantation.rds",sep=""))

saveRDS(List_imp7,file=paste(saveext,"/Stroma_Ourimplantation_vs_prolif.rds",sep=""))
saveRDS(List_imp8,file=paste(saveext,"/Stroma_prolif_vs_ourImplantation.rds",sep=""))


saveRDS(List_imp9,file=paste(saveext,"/Stroma_sec_vs_prolif.rds",sep=""))
saveRDS(List_imp10,file=paste(saveext,"/Stroma_prolif_vs_sec.rds",sep=""))

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(List_imp1)[1:200] ), name = "ImpS1")
p1<-FeaturePlot(mammal.combined3, features = "ImpS11", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_ImpS1_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(List_imp2)[1:200] ), name = "ImpsecvS1")
p1<-FeaturePlot(mammal.combined3, features = "ImpsecvS11", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_ImpS2_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(List_imp3)[1:200] ), name = "ImpS2")
p1<-FeaturePlot(mammal.combined3, features = "ImpS21", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_ImpS3_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)


mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(List_imp4)[1:200] ), name = "ImpprolifS2")
p1<-FeaturePlot(mammal.combined3, features = "ImpprolifS21", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_ImpS4_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

Mat <- as.data.frame(t(rbind(mammal.combined3$dS1,mammal.combined3$prolifS1,mammal.combined3$ImpS21)))
test <- apply(Mat, 1, function(x) which(x == max(x), arr.ind = TRUE))
IDs <- c("Dec","Prolf","Imp2")
newLabel <- IDs[unlist(test)]
newLabel[1:length(newLabel)] <- "Prolif"
newLabel[which(mammal.combined3$dS1>mammal.combined3$prolifS1 & mammal.combined3$dS1>mammal.combined3$ImpS11)] <- "Decid"
#Now refine by implantation type
newLabel[which(mammal.combined3$dS1>mammal.combined3$prolifS1 & mammal.combined3$ImpS11>mammal.combined3$ImpsecvS11)] <- "Impl"
newLabel[which(mammal.combined3$dS1<mammal.combined3$prolifS1 & mammal.combined3$ImpS21>mammal.combined3$ImpprolifS21)] <- "Impl"
Idents(mammal.combined3) <-newLabel
mammal.combined3$StromaType <- Idents(mammal.combined3)

p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/StromaType2.pdf",sep=""),width = 30, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims=c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/StromaType2_1_3.pdf",sep=""),width = 30, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims=c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/StromaType2_1_4.pdf",sep=""),width = 30, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims=c(1,5)) 
ggsave(filename=paste(saveext,"/DimRed/StromaType2_1_5.pdf",sep=""),width = 30, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims=c(1,6)) 
ggsave(filename=paste(saveext,"/DimRed/StromaType2_1_6.pdf",sep=""),width = 30, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined3,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/StromaType2_nosplit.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)

cl1 <- intersect(rownames(List_imp1),rownames(List_imp3))
cl2 <- intersect(rownames(List_imp2),rownames(List_imp4))

newdata1t <- List_imp1[cl1,]
newdata1t <- newdata1t[order(newdata1t$avg_log2FC,decreasing = TRUE),]
newdata2t <- List_imp2[cl2,]
newdata2t <- newdata2t[order(newdata2t$avg_log2FC,decreasing = TRUE),]

newdata3t <- List_imp3[cl1,]
newdata3t <- newdata3t[order(newdata3t$avg_log2FC,decreasing = TRUE),]
newdata4t <- List_imp4[cl2,]
newdata4t <- newdata4t[order(newdata4t$avg_log2FC,decreasing = TRUE),]

AllData5 <- AllData4
Idents(AllData5,cells=WhichCells(AllData5,idents=c("dS2","dS3"))) <- "dS1"
Idents(AllData5,cells=WhichCells(AllData5,idents=c("proliferativeeS","proliferativedS"))) <- "pStrom"
Idents(AllData5,cells=WhichCells(AllData5,idents=c("late-secretorydS","late-secretoryeS"))) <- "secStrom"

write.table(as.data.frame(newdata1t),file=paste(saveext,"StromaImplantationDE1.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata2t),file=paste(saveext,"StromaImplantationDE2.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata3t),file=paste(saveext,"StromaImplantationDE3.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata4t),file=paste(saveext,"StromaImplantationDE4.rds",sep=""),sep=",",quote = FALSE)

DefaultAssay(AllData5) <- "RNA"
AllData5$ID4 <- NULL
AllData5$ID1 <- NULL
AllData5$ID2 <- NULL
AllData5$ID3 <- NULL

IDss<- as.character(Idents(AllData5))
IDss[which(is.na(Idents(AllData5))==1 )] <- "Other"
Idents(AllData5) <- IDss
Idents(AllData5) <- factor(x = Idents(AllData5), levels = c("pS","dS","pStrom","secStrom","dS1","early-secretorydS","early-secretoryeS","Other"))

List_imp1A <- FindMarkers(AllData5,ident.1 = c("dS1"),ident.2 = c("secStrom"), test.use = "MAST")
List_imp3A <- FindMarkers(AllData5,ident.1 = c("dS1"),ident.2 = c("pStrom"), test.use = "MAST" )
List_imp5A <- FindMarkers(AllData5,ident.1 = c("dS"),ident.2 = c("secStrom"), test.use = "MAST" )
List_imp7A <- FindMarkers(AllData5,ident.1 = c("dS"),ident.2 = c("pStrom"), test.use = "MAST")

Idents(Endo0,WhichCells(Endo0,idents=c("dS2","dS3"))) <- "dS1"
M1 <- AverageExpression(Endo0)
Idents(Endo4,WhichCells(Endo4,idents=c("proliferativeeS","proliferativedS"))) <- "pStrom"
Idents(Endo4,WhichCells(Endo4,idents=c("late-secretorydS","late-secretoryeS"))) <- "secStrom"
M2 <- AverageExpression(Endo4)
M3 <- AverageExpression(onlyOurs)

All
AvExpS <- AverageExpression(AllData5, assays = "RNA")
AvExpS <- as.data.frame(AvExpS$RNA)

AvExpS$X1 <- NA
AvExpS$Y1 <- NA
AvExpS$X2 <- NA
AvExpS$Y2 <- NA
AvExpS$X3 <- NA
AvExpS$Y3 <- NA
AvExpS$X4 <- NA
AvExpS$Y4 <- NA
AvExpS[rownames(List_imp1A),"X1"] <- List_imp1A$avg_log2FC
AvExpS[rownames(List_imp1A),"Y1"] <- List_imp1A$p_val_adj
AvExpS[rownames(List_imp3A),"X2"] <- List_imp3A$avg_log2FC
AvExpS[rownames(List_imp3A),"Y2"] <- List_imp3A$p_val_adj
AvExpS[rownames(List_imp5A),"X3"] <- List_imp5A$avg_log2FC
AvExpS[rownames(List_imp5A),"Y3"] <- List_imp5A$p_val_adj
AvExpS[rownames(List_imp7A),"X4"] <- List_imp7A$avg_log2FC
AvExpS[rownames(List_imp7A),"Y4"] <- List_imp7A$p_val_adj

AvExpS <- AvExpS[intersect(rownames(Endo0),rownames(Endo4)),] 

p1 <- ggplot(AvExpS, aes(X1,X3)) + geom_point() + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) + xlim(c(-10, 10))
genes.to.label <- genes.to.label <- unique(intersect(rownames(List_imp1A),rownames(List_imp5A)) )
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"Scatter_Imp_vs_Sec.pdf",sep=""),width = 40, height = 40, plot = p1, useDingbats=FALSE)

p1 <- ggplot(AvExpS, aes(X2,X4)) + geom_point() + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) + xlim(c(-10, 10))
genes.to.label <- genes.to.label <- unique(intersect(rownames(List_imp3A),rownames(List_imp7A)))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"Scatter_Imp_vs_Prolif.pdf",sep=""),width = 40, height = 40, plot = p1, useDingbats=FALSE)


list2 <- c("AXL","TYRO3","MERTK","PROS1","GAS6")


p1<-VlnPlot(AllData5,  idents = c("pS","dS","pStrom","secStrom","dS1"),features = list2, pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_AXL.pdf",sep=""),width = 40, height = 20, limitsize = FALSE,p1)


p1<-VlnPlot(AllData4,features = list2, pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplEpi_AXL.pdf",sep=""),width = 40, height = 20, limitsize = FALSE,p1)


list <- c("NFKBIA","TUBA1B","SPOCK1","TNFAIP3","HMGN2","GLRX","SOD2","TUBB","ID1","RND3","NES","WNT4")

p1<-VlnPlot(AllData5,  idents = c("pS","dS","pStrom","secStrom","dS1"),features = list, pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("pStrom","secStrom","Stromal fibroblasts","dS1"),features = rownames(newdata2t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_DownReg2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("pStrom","secStrom","Stromal fibroblasts","dS1"),features = rownames(newdata3t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_DownReg3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("pStrom","secStrom","Stromal fibroblasts","dS1"),features = rownames(newdata4t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_DownReg4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)


#Not do stroma
List7 <- FindMarkers(Endo4,ident.2 = "proliferativeeS",ident.1 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List8 <- FindMarkers(Endo4,ident.1 = "proliferativeeS",ident.2 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List9 <- FindMarkers(Endo4,ident.1 = "proliferativeFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List10 <- FindMarkers(Endo4,ident.1 = "late-secretoryFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List11 <- FindMarkers(Endo4,ident.2 = "proliferativeFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List12 <- FindMarkers(Endo4,ident.2 = "late-secretoryFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )


List13 <- FindMarkers(Endo4,ident.2 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.1 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List14 <- FindMarkers(Endo4,ident.1 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.2 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )


newdata7 <- List7[order(List7$avg_log2FC,decreasing = TRUE),]
newdata8 <- List8[order(List8$avg_log2FC,decreasing = TRUE),]
newdata9 <- List9[order(List9$avg_log2FC,decreasing = TRUE),]
newdata10 <- List10[order(List10$avg_log2FC,decreasing = TRUE),]
newdata11 <- List11[order(List11$avg_log2FC,decreasing = TRUE),]
newdata12 <- List12[order(List12$avg_log2FC,decreasing = TRUE),]

newdata13 <- List13[order(List13$avg_log2FC,decreasing = TRUE),]
newdata14 <- List14[order(List14$avg_log2FC,decreasing = TRUE),]


#Now load our stroma

mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_withModuleScores.rds",sep=""))
Idents(mammal.combined) <- mammal.combined$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
uID <- as.character( Idents(onlyOurs) )
uID[1:length(uID)] <- "dS"
uID[which(onlyOurs$prolifS1 > onlyOurs$dS1)] <- "pS"
Idents(onlyOurs) <- uID
DefaultAssay(onlyOurs) <- "RNA"

uID <- as.character( Idents(onlyOurs) )
uID[1:length(uID)] <- "dS"
uID[which(onlyOurs$prolifS1 > onlyOurs$dS1)] <- "pS"
Idents(onlyOurs) <- uID
DefaultAssay(onlyOurs) <- "RNA"

AllData3 <- merge(onlyOurs,y=c(Endo4,Endo0,assemb),project = "Merged")

AllData4 <- subset(AllData3,idents=c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))


DefaultAssay(AllData4) <- "RNA"
AllData4$ID4 <- NULL
AllData4$ID1 <- NULL
AllData4$ID2 <- NULL
AllData4$ID3 <- NULL

p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata7)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata8)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata13)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata14)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata7)[41:80], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set1_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata8)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set2_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata13)[41:80], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata14)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



AllData3 <- merge(onlyOurs,y=c(Endo4,Endo0,assemb),project = "Merged")

AllData4 <- subset(AllData3,idents=c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))


Idents(AllData4,cells=WhichCells(AllData4,idents=c("proliferativeeS","early-secretoryeS","late-secretoryeS"))) <- "proliferativeeS"
Idents(AllData4,cells=WhichCells(AllData4,idents=c("proliferativedS","early-secretorydS","late-secretorydS"))) <- "late-secretorydS"
Idents(AllData4,cells=WhichCells(AllData4,idents=c("dS2","dS3"))) <- "dS2"

AllData4 <- subset(AllData4,idents=c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"))



Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("pS","dS","proliferativeeS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"))

DefaultAssay(AllData4) <- "RNA"
AllData4$ID4 <- NULL
AllData4$ID1 <- NULL
AllData4$ID2 <- NULL
AllData4$ID3 <- NULL

p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata7)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set1_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata8)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set2_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata13)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata14)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



#Not do stroma
List7 <- FindMarkers(Endo4,ident.2 = "proliferativeeS",ident.1 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List8 <- FindMarkers(Endo4,ident.1 = "proliferativeeS",ident.2 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List9 <- FindMarkers(Endo4,ident.1 = "proliferativeFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List10 <- FindMarkers(Endo4,ident.1 = "late-secretoryFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List11 <- FindMarkers(Endo4,ident.2 = "proliferativeFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List12 <- FindMarkers(Endo4,ident.2 = "late-secretoryFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )


List13 <- FindMarkers(Endo4,ident.2 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.1 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List14 <- FindMarkers(Endo4,ident.1 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.2 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )

List15 <- FindMarkers(onlyOurs,ident.2 = c("pS"),ident.1 = c("dS"), test.use = "MAST", only.pos = TRUE )
List16 <- FindMarkers(onlyOurs,ident.1 = c("pS"),ident.2 = c("dS"), test.use = "MAST", only.pos = TRUE )

int1 <- intersect(rownames(List13),rownames(List15))
int2 <- intersect(rownames(List14),rownames(List16))

newdata15 <- List13[int1,]
  
newdata15 <- newdata15[order(newdata15$avg_log2FC,decreasing = TRUE),]

newdata16 <- List14[int2,]
newdata16 <- newdata16[order(newdata16$avg_log2FC,decreasing = TRUE),]




p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata15)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3_mergecommon.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata16)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4_mergecommon.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



onlyOurs <- subset(mammal.combined,cells=ours)
notOurs <- subset(mammal.combined,cells=ours,invert=TRUE)

DefaultAssay(onlyOurs) <- "RNA"
DefaultAssay(notOurs) <- "RNA"
List1 <- FindMarkers(onlyOurs,ident.2 = "Lumenal",ident.1 = c("Glandular"), test.use = "MAST")
List2 <- FindMarkers(onlyOurs,ident.2 = "Lumenal",ident.1 = c("Ciliated"), test.use = "MAST")
List3 <- FindMarkers(onlyOurs,ident.2 = "Glandular",ident.1 = c("Ciliated"), test.use = "MAST")
List4 <- FindMarkers(notOurs,ident.2 = c("proliferative_Lumenal","early-secretory_Lumenal","late-secretory_Lumenal"),ident.1 = c("proliferative_Glandular","early-secretory_Glandular","late-secretory_Glandular"), test.use = "MAST", only.pos = TRUE )
List5 <- FindMarkers(notOurs,ident.2 = c("proliferative_Lumenal","early-secretory_Lumenal","late-secretory_Lumenal"),ident.1 = c("Ciliated"), test.use = "MAST", only.pos = TRUE )
List6 <- FindMarkers(notOurs,ident.2 = c("proliferative_Glandular","early-secretory_Glandular","late-secretory_Glandular"),ident.1 = c("Ciliated"), test.use = "MAST", only.pos = TRUE )


newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List2$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4 <- List4[order(List4$avg_log2FC,decreasing = TRUE),]
newdata5 <- List5[order(List5$avg_log2FC,decreasing = TRUE),]
newdata6 <- List6[order(List6$avg_log2FC,decreasing = TRUE),]


AvE <- AverageExpression(onlyOurs)
AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$Lumenal+AvExp$Glandular)/2
AvExp$'log2FC' <- NA

AvExp[rownames(List1),'log2FC'] <- List1$avg_log2FC
AvExp[rownames(List1),'pval'] <- List1$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')

genes.to.label1 = unique(c(intersect(rownames(List1),rownames(List4))))

p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"LumenGland_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)




AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$Lumenal+AvExp$Ciliated)/2
AvExp$'log2FC' <- NA

AvExp[rownames(List2),'log2FC'] <- List2$avg_log2FC
AvExp[rownames(List2),'pval'] <- List2$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')

genes.to.label1 = unique(c(intersect(rownames(List2),rownames(List5))))

p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"LumenCiliated_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)


AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$Glandular+AvExp$Ciliated)/2
AvExp$'log2FC' <- NA

AvExp[rownames(List3),'log2FC'] <- List3$avg_log2FC
AvExp[rownames(List3),'pval'] <- List3$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')

genes.to.label1 = unique(c(intersect(rownames(List3),rownames(List6))))

p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"GlandularCiliated_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)





