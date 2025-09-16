library("Seurat")
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
library(stringr)
library(gtools)

set.seed(1) 

saveext = "/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping5/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#oldMarmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/Cyno3D/data/Plots/Marmoset_spatialSS2_allshots.rds")

#Set a key for the colours
my_cols <- c('EmDisc_CS5'='#0c9cf5','EmDisc_CS6'='#0767DA','EmDisc_CS7'='#0233BF','Epi_CS3'='#00BFBF',
             'Am_CS5'='#877bd6','Am_CS6'='#5F54C7','Am_CS7'='#1A0873',
             'PGC_CS5'='#E6E600','PGC_CS6'='#BFBF04','PGC_CS7'='#999903',
             'VE_CS5'='#F04C04','VE_CS6'='#D74404','VE_CS7'='#BF3C04','Hyp_CS3'='#E6B500',
             'SYS_CS5'='#E68600','SYS_CS6'='#d17600','SYS_CS7'='#BF7104',
             'ExMes_CS5'='#e6c800','ExMes_CS6'='#c49a00','ExMes_CS7'='#967700','Stalk_CS6'='#754C24','Stalk_CS7'='#603813',
             'Tb_CS5'='#921FE6','Tb_CS6'='#8017c2','Tb_CS7'='#7108a6','Tb_CS3'='#BF0489','newTSP3'='#BF0489',
             'cmTSCPAVS'='#7108a6','3D post'='#e6c800','6D post'='#c49a00','3D post BMP4'='#967700',
             'Amspheroids'='#877bd6','LN no MEF'='#1A0873', 'Epi BMP4'='#5F54C7', 'Epispheroids'='#0c9cf5', 'LN MEF'='#0767DA', '24h post'='#0233BF','PGC'='#BFBF04',
             'PGCLC_ExMes'='#e6c800', 'PGCLC_Am'='#877bd6', 'PGCL'='#E6E600','PGC'='#E6E600','PLAXA'='#cfecff')


#marmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping3/Matteo_aligntoMarm_and_cyno_noPreimp.rds")


#/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/AligntoMarm_latestAnowTb.rds 
marmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping4/AlignMatteo_withMarmosetwCynoMarmHumanSS.rds")


Dataset <- as.character(marmoset$Dataset)
Dataset[which(Dataset%in%c("Batch1","10X Ours","Batch3","Batch6","Batch7","Batch8"))] <- "InVitro"
marmoset$Dataset <- Dataset

#
library(destiny)
#mammal.combined <- subset(mammal.combined,cells=C2outlier1,invert=TRUE)
#mammal.combined <- subset(mammal.combined,cells=C2outlier2,invert=TRUE)
DN1 <- as.data.frame(marmoset[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(marmoset)
marmoset[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(marmoset))
p <- DimPlot(marmoset, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

ID3 <- as.character(marmoset$ID3)
ID3[which(is.na(ID3)==1)] <- "Cyno"
marmoset$ID3 <- ID3
p <- DimPlot(marmoset, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_split",".pdf",sep=""),width = 80, height = 8, limitsize = FALSE,p)

marmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/AligntoMarm_latestAno.rds")
marmoset$stashID <- Idents(marmoset)
Idents(marmoset) <- marmoset$Dataset
mammal.combinedsubs1 <- subset(marmoset,idents=c("InVitro"))
mammal.combinedsubs2 <- subset(marmoset,idents=c("2) Marmoset in vivo"))
mammal.combinedsubs3 <- subset(marmoset,idents=c("Cyno"))

p <- DimPlot(mammal.combinedsubs1, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_Matteo",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(mammal.combinedsubs2) <- mammal.combinedsubs2$stashID
p <- DimPlot(mammal.combinedsubs2, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_Marm",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(mammal.combinedsubs3, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_Cyno",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

intersect(c(intersect(WhichCells(mammal.combinedsubs1,expression=NANOS3>0.1), WhichCells(mammal.combinedsubs1,expression=TFAP2C>0.1)),intersect(WhichCells(mammal.combinedsubs1,expression=NANOS3>0.1), WhichCells(mammal.combinedsubs1,expression=TFAP2A>0.1))),intersect(WhichCells(mammal.combinedsubs1,expression=NANOS3>0.1), WhichCells(mammal.combinedsubs1,expression=SOX17>0.1)))

p <- FeatureScatter(mammal.combinedsubs1,feature1 = "NANOS3" ,feature2 = "TFAP2C")
ggsave(filename=paste(saveext,"/DimRed/PGC",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

DefaultAssay(mammal.combinedsubs1) <- "RNA"
p <- FeatureScatter(mammal.combinedsubs1,feature1 = "NANOS3" ,feature2 = "SOX17")
ggsave(filename=paste(saveext,"/DimRed/PGC3",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


DefaultAssay(mammal.combinedsubs1) <- "RNA"
p <- FeatureScatter(mammal.combinedsubs1,feature1 = "NANOS3" ,feature2 = "POU5F1")
ggsave(filename=paste(saveext,"/DimRed/PGC4",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


p <- FeatureScatter(mammal.combinedsubs1,feature1 = "NANOS3" ,feature2 = "TFAP2A")
ggsave(filename=paste(saveext,"/DimRed/PGC2",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

mammal.combined <- FindClusters(mammal.combined,resolution = 1)
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_Marm",".pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)

DefaultAssay(mammal.combined) <- "RNA"
p <- FeaturePlot(mammal.combined, reduction = "dm", feature = "TFAP2C",pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_Marm_TFAP2C",".pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)

DefaultAssay(mammal.combinedsubs1) <- "integrated"
mammal.combinedsubs1 <- FindClusters(mammal.combinedsubs1,resolution = 1)
p <- DimPlot(mammal.combinedsubs1, reduction = "dm", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_Ours_Cluster",".pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)

p <- DimPlot(mammal.combinedsubs1, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM1_Ours_Cluster_splt",".pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p <- DimPlot(marmoset, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DM1",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

p <- DimPlot(marmoset, reduction = "umap", pt.size = 2, label.size = 2, split.by="Species",label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DM2",".pdf",sep=""),width = 40, height = 8, limitsize = FALSE,p)



#marmoset <- readRDS(file=paste("~/Desktop/Data/Endometrial/InVitro/Matteo/MatteoLatestAlign_to_marmosetwithMaternal.rds",sep=""))

#marmoset$cellshape <- as.factor(as.double(as.factor(marmoset$Dataset)))
#marmoset$cellshape <- as.factor(2*as.double(as.factor(marmoset$Dataset)))

#p<-DimPlot(marmoset,  cols = my_cols, shape.by = 'cellshape', pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/AlignedData_PCA_withPreImpl.pdf",sep=""),width = 20, height = 10,p)

#p<-DimPlot(marmoset, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/AlignedData_PCA_withPreImpl.pdf",sep=""),width = 20, height = 10,p)

#Dataset <- as.character(marmoset$orig.ident)
#Dataset2 <- as.character(marmoset$rough.ident)
#Dataset[which(Dataset%in%c("D10","D12","D6","D8","CER1","cRCER1","cRH9","Day5-4-pTEICM","Day5-8-pTEICM","hm","hv","Day6-2-pTEICM","Day6-3-pTEICM","Day6-4-pTEICM","Day6-5-pTEICM","Day7-2-pTEICM","Day7-3-pTEICM","Day7-6-pTEICM","H9","hE"))] <- "Other"
#Dataset[which(Dataset%in%c("marmoset"))] <- "MarmInVivo"
#Dataset[which(Dataset%in%c("Petro","Tyser","Wang","Xiang"))] <- "HumanRef"
#Dataset[which(Dataset2%in%c("CER1_primed","CER1_PXGL","cRCER1_ACL_CER1-","cRCER1_ACL_CER1-_AB","cRCER1_ACL_CER1+","cRCER1_ACL_CER1+_AB","cRH9_ACL_PD-","cRH9_ACL_PD-_AB","cRH9_ACL_PD+","cRH9_ACL_PD+_AB"))] <- "InVitro"


#marmoset$Dataset <- Dataset

#TempID <- as.character(marmoset$Lineage)
#TempID[which(is.na(TempID)==1)] <- "Other"
#TempID[which(Dataset2%in%c("CER1_primed","CER1_PXGL","cRCER1_ACL_CER1-","cRCER1_ACL_CER1-_AB","cRCER1_ACL_CER1+","cRCER1_ACL_CER1+_AB","cRH9_ACL_PD-","cRH9_ACL_PD-_AB","cRH9_ACL_PD+","cRH9_ACL_PD+_AB"))] <- Dataset2[which(Dataset2%in%c("CER1_primed","CER1_PXGL","cRCER1_ACL_CER1-","cRCER1_ACL_CER1-_AB","cRCER1_ACL_CER1+","cRCER1_ACL_CER1+_AB","cRH9_ACL_PD-","cRH9_ACL_PD-_AB","cRH9_ACL_PD+","cRH9_ACL_PD+_AB"))]
#Idents(marmoset) <- TempID


#p<-DimPlot(marmoset, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/AlignedData_PCA_withPreImplID.pdf",sep=""),width = 20, height = 10,p)

#p<-DimPlot(marmoset, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/AlignedData_UMAP_withPreImplID.pdf",sep=""),width = 20, height = 10,p)
#marmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/AligntoMarm_latestAnowTb.rds")

#marmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping2/AlignMatteo_withMarmosetwMaternal.rds")
marmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping2/AlignMatteo_withMarmosetwMaternal.rds")

Dataset <- as.character(marmoset$Dataset)
Dataset[which(Dataset%in%c("Batch1","10X Ours","Batch3","Batch6","Batch7","Batch8"))] <- "InVitro"
marmoset$Dataset <- Dataset

#Set number of neighbours to base projections on
k = 10

#Now generate the mappings. First generate KNN/SNN graphs for various resolutions
marmoset <- FindNeighbors(marmoset, reduction = "pca", dims = 1:20, k.param = 100)
KNN_K100 <- marmoset@graphs$integrated_nn
SNN_K100 <- marmoset@graphs$integrated_snn
saveRDS(KNN_K100,file=paste(saveext,"KNN_K100_Batch1_withPreImpl.rds",sep=""))
saveRDS(SNN_K100,file=paste(saveext,"SNN_K100_Batch1_withPreImpl.rds",sep=""))

marmoset <- FindNeighbors(marmoset, reduction = "pca", dims = 1:20, k.param = 50)
KNN_K50 <- marmoset@graphs$integrated_nn
SNN_K50 <- marmoset@graphs$integrated_snn
saveRDS(KNN_K50,file=paste(saveext,"KNN_K50_Batch1_withPreImpl.rds",sep=""))
saveRDS(SNN_K50,file=paste(saveext,"SNN_K50_Batch1_withPreImpl.rds",sep=""))

marmoset <- FindNeighbors(marmoset, reduction = "pca", dims = 1:20, k.param = 150)
KNN_K150 <- marmoset@graphs$integrated_nn
SNN_K150 <- marmoset@graphs$integrated_snn
saveRDS(KNN_K150,file=paste(saveext,"KNN_K150_Batch1_withPreImpl.rds",sep=""))
saveRDS(SNN_K150,file=paste(saveext,"SNN_K150_Batch1_withPreImpl.rds",sep=""))

marmoset <- FindNeighbors(marmoset, reduction = "pca", dims = 1:20, k.param = 30)
KNN_K30 <- marmoset@graphs$integrated_nn
SNN_K30 <- marmoset@graphs$integrated_snn
saveRDS(KNN_K30,file=paste(saveext,"KNN_K30_Batch1_withPreImpl.rds",sep=""))
saveRDS(SNN_K30,file=paste(saveext,"SNN_K30_Batch1_withPreImpl.rds",sep=""))

#Index for reference dataset vs mapping dataset
inds1 <- which(marmoset$Dataset%in%c("InVitro") ) 
inds <- which(marmoset$Dataset=="2) Marmoset in vivo")


#NewMeta <- data.frame(LOC=oldMarmoset$LOC)
#rownames(NewMeta) <- names(oldMarmoset$LOC)
#marmoset <- AddMetaData(marmoset, NewMeta)
#marmoset$LOC <- marmoset$LOC
#LOC <- as.character(marmoset$LOC)
#LOC[which(is.na(LOC)==1)] <- NA

#Empty arrays for storing
neigh <- matrix(, ncol = k, nrow = length(Idents(marmoset)))
C1 <- matrix(, ncol = k, nrow = length(Idents(marmoset)))
S1 <- matrix(, ncol = k, nrow = length(Idents(marmoset)))
C1p <- matrix(, ncol = k, nrow = length(Idents(marmoset)))
lC1 <- matrix(, ncol = k, nrow = length(Idents(marmoset)))
lC1p <- matrix(, ncol = k, nrow = length(Idents(marmoset)))

#Spatial loctions
Locations <- as.data.frame(marmoset$LOC)

#Numbers of the various cell types in the referene datasets
basecounts <- table(Idents(marmoset)[which(marmoset$Dataset=="2) Marmoset in vivo")])

#More empty array
neigh <- matrix(,  ncol = length(basecounts)  ,nrow = length(Idents(marmoset)) )

#Save current ID just in case
marmoset$oldID <- Idents(marmoset)

#Loop over all cells
for (i in 1:length(Idents(marmoset)) ) {
  #Find SNNs and get cell IDs
  a1 <- sort(SNN_K50[i,inds],decreasing = TRUE)[1:k]
  C1[i,] <- names(a1)
  #Count number of cell types in the neighbourhood
  neigh[i,] <- table(Idents(marmoset)[names(a1)])
  #Store values
  S1[i,] <- a1
  lC1[i,] <- as.character(Locations[names(a1),])
  
  #Permutate the SNN graph, count number of instances and store
  aperm <- SNN_K50[i,inds]
  names(aperm) <- gtools::permute(names(aperm))
  a2 <- sort(aperm,decreasing = TRUE)[1:k]
  C1p[i,] <- names(a2)
  lC1p[i,] <- as.character(Locations[names(a2),])
}

#Write out the overlaps and various other things
write.table(Idents(marmoset)[inds1],file=paste(saveext,'AllMapIdents_set1_all_withPreImpl_byCl.csv',sep=""),sep=",")
write.table(marmoset$LOC[inds1],file=paste(saveext,'AllMapIdents_set1_all_withPreImpl_byCl.csv',sep=""),sep=",")
write.table(Idents(marmoset)[inds],file=paste(saveext,'AllMapIdents_set2_all_withPreImpl_byCl.csv',sep=""),sep=",")
write.table(marmoset$LOC[inds],file=paste(saveext,'AllMapIdents_set2_all_withPreImpl_byCl.csv',sep=""),sep=",")
write.table(data.frame(as.data.frame(colnames(marmoset)[inds1]),as.data.frame(marmoset$LOC[inds1]),as.data.frame(Idents(marmoset)[inds1]), as.data.frame(C1[inds1,]), as.data.frame(lC1[inds1,]), 1 ),file=paste(saveext,'C1_set1_withPreImpl_byCl.csv',sep=""),sep=",",row.names = F)
write.table(data.frame(as.data.frame(colnames(marmoset)[inds1]),as.data.frame(marmoset$LOC[inds1]),as.data.frame(Idents(marmoset)[inds1]), as.data.frame(C1p[inds1,]),as.data.frame(lC1p[inds1,]),1 ),file=paste(saveext,'C1perm_set1_withPreImpl_byCl.csv',sep=""),sep=",",row.names = F)
write.table(S1[inds1,],file=paste(saveext,'S1_set1_all_withPreImpl_byCl.csv',sep=""),sep=",")
write.table(marmoset$LOC,file=paste(saveext,'AllLocKey_all_withPreImpl_byCl.csv',sep=""),sep=",")
write.table(Idents(marmoset),file=paste(saveext,'AllLocIDKey_all_withPreImpl_byCl.csv',sep=""),sep=",")

#Files for running the 'correlation based mapping' based on integrated dataset
DefaultAssay(marmoset) <- "integrated"
d0 <- GetAssayData(marmoset, assay = "integrated")
d1 <- d0[,which(marmoset$Dataset=="2) Marmoset in vivo")]
d2 <- d0[,which(marmoset$Dataset%in%c("InVitro") )]

Locs1 <- marmoset$LOC[which(marmoset$Dataset=="2) Marmoset in vivo")] 
Locs2 <- marmoset$LOC[which(marmoset$Dataset%in%c("InVitro") )] 

ID1 <- Idents(marmoset)[which(marmoset$Dataset=="2) Marmoset in vivo")] 
ID2 <- Idents(marmoset)[which(marmoset$Dataset%in%c("InVitro") )] 

write.csv(data.frame(Loc=Locs1,Type=ID1,z=1), file=paste(saveext,"/Correlation_int_locations1_wall_withPreImpl_byCl.csv",sep=""))
write.csv(data.frame(Loc=Locs2,Type=ID2,z=1), file=paste(saveext,"/Correlation_int_locations2_wall_withPreImpl_byCl.csv",sep=""))

#Finally calculate the correlations all by all and save
C <- cor(as.data.frame(d1),as.data.frame(d2),method="pearson")
write.csv(as.data.frame(C), file=paste(saveext,"/Correlation_int_set1_wall_withPreImpl_byCl.csv",sep=""))






#####
Epithel <- subset(Matteo,idents=c("EpithelialOCT4","Epithelial","Unknown","putEmbryonic"))
Idents(Epithel) <- Epithel$Dataset
Epithel1 <- subset(Epithel,idents="Batch1")

Idents(Epithel1) <- Epithel1$HarmAno
Epithel1 <- FindVariableFeatures(Epithel1, selection.method = "vst", nfeatures = 20000)
Epithel1 <- ScaleData(Epithel1, verbose = FALSE)
Epithel1 <- RunPCA(Epithel1, npcs = 20, verbose = FALSE)
Epithel1 <- RunUMAP(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- RunTSNE(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- FindNeighbors(Epithel1, reduction = "pca", dims = 1:20)
#p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
#ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
#p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
#ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
#Idents(Epithel1) <- Epithel1$Genotype
#p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
#ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
#p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
#ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
#ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAP_Markers",".pdf",sep=""),width = 80, height = 80, limitsize = FALSE,p)
#p<-FeaturePlot(Epithel1, features =c("MSMB",
#"CST6",
#"CRISP3")  )
#ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAP_Markers",".pdf",sep=""),width = 80, height = 80, limitsize = FALSE,p)
#p<-FeaturePlot(Epithel1, features = c("POU5F1","SOX17","CER1","KRT7") ,reduction = "umap" )
#ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1PCA_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)

Epithel1 <- FindClusters(Epithel1,resolution = 3)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAP_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)

#POU5F1>2 = top 0.04%

#Cluster 12 and 16???
Epithel1 <- FindVariableFeatures(Epithel1, selection.method = "vst", nfeatures = 40000)
Epithel1 <- ScaleData(Epithel1, verbose = FALSE)
Epithel1 <- RunPCA(Epithel1, npcs = 20, verbose = FALSE)
Epithel1 <- RunUMAP(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- RunTSNE(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- FindNeighbors(Epithel1, reduction = "pca", dims = 1:20)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Epithel1) <- Epithel1$Genotype
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Epithel1, features = c("POU5F1","LAMB1","FN1","CER1") , reduction = "pca" )
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)
Epithel1 <- FindClusters(Epithel1,resolution = 3)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1UMAP_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "tsne", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1TSNE_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch1PCA_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
#Cluster 17 is distinct in PCA.

B1List1 <- WhichCells(Epithel1,expression=POU5F1>2)
B1List1A <- WhichCells(Epithel1,expression=POU5F1>2)

B1List2 <- intersect(B1List1A,colnames(Epithel1)[which(Epithel1$Genotype%in%c("NAssigned"))])
B1List3 <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("Embryonic_G1"))] #
B1List4 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("17"))]
B1List5 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("11"))] #Embryonic

#table(Epithel1$Genotype[which(Idents(Epithel1)%in%c("17"))])


#####
Epithel <- subset(Matteo,idents=c("EpithelialOCT4","Epithelial","Unknown","putEmbryonic"))
Idents(Epithel) <- Epithel$Dataset
Epithel1 <- subset(Epithel,idents="10X Ours")

Epithel1 <- FindVariableFeatures(Epithel1, selection.method = "vst", nfeatures = 40000)
Epithel1 <- ScaleData(Epithel1, verbose = FALSE)
Epithel1 <- RunPCA(Epithel1, npcs = 20, verbose = FALSE)
Epithel1 <- RunUMAP(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- RunTSNE(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- FindNeighbors(Epithel1, reduction = "pca", dims = 1:20)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Epithel1) <- Epithel1$Genotype
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Epithel1, features = c("POU5F1","CLDN6","SOX15","CER1")  )
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)
Epithel1 <- FindClusters(Epithel1,resolution = 3)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2UMAP_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "tsne", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2TSNE_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch2PCA_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
#Cluster 17 is distinct in PCA.

#B2List1 <- WhichCells(Epithel1,expression=POU5F1>2)
#B1List2 <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("NAssigned"))]
#B1List3 <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("Embryonic_G1"))] #
#B1List4 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("17"))]
#B1List5 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("11"))] #Embryonic
#For batch 2 there are no high OCT4 clusters. No putative embryonic genoypes. Done.



Epithel1 <- subset(Epithel,idents="Batch3")
Idents(Epithel1) <- Epithel1$HarmAno
Epithel1 <- FindVariableFeatures(Epithel1, selection.method = "vst", nfeatures = 25000)
Epithel1 <- ScaleData(Epithel1, verbose = FALSE)
Epithel1 <- RunPCA(Epithel1, npcs = 20, verbose = FALSE)
Epithel1 <- RunUMAP(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- RunTSNE(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- FindNeighbors(Epithel1, reduction = "pca", dims = 1:20)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch3",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch3UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Epithel1) <- Epithel1$Genotype
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch3GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch3UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Epithel1, features = c("POU5F1","APOA1","FN1","CER1")  ,  reduction = "pca" )
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch3UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)
Epithel1 <- FindClusters(Epithel1,resolution = 3)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch3UMAP_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch3PCA_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)


B3List1 <- WhichCells(Epithel1,expression=POU5F1>2)
B3List1A <- WhichCells(Epithel1,expression=POU5F1>1)
B3List2 <- intersect(B3List1A,colnames(Epithel1)[which(Epithel1$Genotype%in%c("NAssigned"))])
B3List3 <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("Embryonic_G5","Embryonic_G4","Epithelial_G1"))]
B3List4 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("16","19"))] #Embryonic
#List5 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("11"))] #Embryonic




Epithel1 <- subset(Epithel,idents="Batch6")
Idents(Epithel1) <- Epithel1$HarmAno
Epithel1 <- FindVariableFeatures(Epithel1, selection.method = "vst", nfeatures = 25000)
Epithel1 <- ScaleData(Epithel1, verbose = FALSE)
Epithel1 <- RunPCA(Epithel1, npcs = 20, verbose = FALSE)
Epithel1 <- RunUMAP(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- RunTSNE(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- FindNeighbors(Epithel1, reduction = "pca", dims = 1:20)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch6",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch6UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Epithel1) <- Epithel1$Genotype
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch6GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch6UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Epithel1, features = c("POU5F1","APOA1","FN1","CER1")  ,  reduction = "pca" )
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch6UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)
Epithel1 <- FindClusters(Epithel1,resolution = 3)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch6UMAP_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch6PCA_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)



B6List1 <- WhichCells(Epithel1,expression=POU5F1>2)
B6List1A <- WhichCells(Epithel1,expression=POU5F1>1)
B6List2 <- intersect(B6List1A,colnames(Epithel1)[which(Epithel1$Genotype%in%c("NAssigned"))])
B6List3 <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("Embryonic_G4","Embryonic_G2"))] #
B6List4 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("1","11","13"))]
#B6List5 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("11"))] #Embryonic



Epithel1 <- subset(Epithel,idents="Batch7")
Idents(Epithel1) <- Epithel1$HarmAno
Epithel1 <- FindVariableFeatures(Epithel1, selection.method = "vst", nfeatures = 25000)
Epithel1 <- ScaleData(Epithel1, verbose = FALSE)
Epithel1 <- RunPCA(Epithel1, npcs = 20, verbose = FALSE)
Epithel1 <- RunUMAP(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- RunTSNE(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- FindNeighbors(Epithel1, reduction = "pca", dims = 1:20)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch7",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch7UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Epithel1) <- Epithel1$Genotype
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch7GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch7UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Epithel1, features = c("POU5F1","APOA1","FN1","CER1")  ,  reduction = "pca" )
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch7UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)
Epithel1 <- FindClusters(Epithel1,resolution = 3)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch7UMAP_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch7PCA_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)


B7List1 <- WhichCells(Epithel1,expression=POU5F1>2)
B7List1A <- WhichCells(Epithel1,expression=POU5F1>1)
B7List2 <- intersect(B7List1A,colnames(Epithel1)[which(Epithel1$Genotype%in%c("NAssigned"))])
B7List3 <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("Embryonic_G2"))] #
#B7List4 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("1","11","13"))]
#B6List5 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("11"))] #Embryonic





Epithel1 <- subset(Epithel,idents="Batch8")
Idents(Epithel1) <- Epithel1$HarmAno
Epithel1 <- FindVariableFeatures(Epithel1, selection.method = "vst", nfeatures = 25000)
Epithel1 <- ScaleData(Epithel1, verbose = FALSE)
Epithel1 <- RunPCA(Epithel1, npcs = 20, verbose = FALSE)
Epithel1 <- RunUMAP(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- RunTSNE(Epithel1, reduction = "pca", dims = 1:20)
Epithel1 <- FindNeighbors(Epithel1, reduction = "pca", dims = 1:20)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch8",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch8UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Epithel1) <- Epithel1$Genotype
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch8GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch8UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Epithel1, features = c("POU5F1","APOA1","FN1","CER1") ,  reduction = "pca" )
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch8UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)
Epithel1 <- FindClusters(Epithel1,resolution = 3)
p <- DimPlot(Epithel1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch8UMAP_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)
p <- DimPlot(Epithel1, reduction = "pca", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/EpithelialBatch8PCA_Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)

B8List1 <- WhichCells(Epithel1,expression=POU5F1>2)
B8List1A <- WhichCells(Epithel1,expression=POU5F1>1)
B8List2 <- intersect(B8List1A,colnames(Epithel1)[which(Epithel1$Genotype%in%c("NAssigned"))] )
B8List3 <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("Embryoic_G3"))] #
#B7List4 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("1","11","13"))]
#B6List5 <- colnames(Epithel1)[which(Idents(Epithel1)%in%c("11"))] #Embryonic


Idents(Matteo,cells=B1List5) <- "putEmbryonic"
Idents(Matteo,cells=WhichCells(Matteo,idents="Unknown")) <- "Epithelial"

#Epithelial cells with a embryonic GT
Idents(Matteo,cells=B1List2) <- "HighOCT4_noGT"
Idents(Matteo,cells=B3List2) <- "HighOCT4_noGT"
Idents(Matteo,cells=B6List2) <- "HighOCT4_noGT"
Idents(Matteo,cells=B7List2) <- "HighOCT4_noGT"
Idents(Matteo,cells=B8List2) <- "HighOCT4_noGT"

Idents(Matteo,cells=B1List3) <- "EpitEmbryonicGT"
Idents(Matteo,cells=B3List3) <- "EpitEmbryonicGT"
Idents(Matteo,cells=B6List3) <- "EpitEmbryonicGT"
Idents(Matteo,cells=B7List3) <- "EpitEmbryonicGT"
Idents(Matteo,cells=B8List3) <- "EpitEmbryonicGT"

Idents(Matteo,cells=B1List1) <- "HighOCT4"
Idents(Matteo,cells=B3List1) <- "HighOCT4"
Idents(Matteo,cells=B6List1) <- "HighOCT4"
Idents(Matteo,cells=B7List1) <- "HighOCT4"
Idents(Matteo,cells=B8List1) <- "HighOCT4"

p<-DimPlot(Matteo, pt.size = 4, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_extraAno",".pdf",sep=""),width = 60, height = 20, limitsize = FALSE,p)
saveRDS(Matteo,paste(saveext,"/MatteoAno_withOCT4ano.rds",sep=""))


Temp1 <- subset(Matteo,idents="Stroma1")
#Idents(Temp1) <- Temp1$Genotype
#p<-DimPlot(Temp1, pt.size = 4, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/Stroma",".pdf",sep=""),width = 60, height = 20, limitsize = FALSE,p)

#saveRDS(Matteo,paste(saveext,"/MatteoAno_Stroma.rds",sep=""))


Idents(Temp1) <- Temp1$Dataset
Stroma1 <- subset(Temp1,idents="Batch1")

Idents(Stroma1) <- Stroma1$HarmAno
Stroma1 <- FindVariableFeatures(Stroma1, selection.method = "vst", nfeatures = 20000)
Stroma1 <- ScaleData(Stroma1, verbose = FALSE)
Stroma1 <- RunPCA(Stroma1, npcs = 20, verbose = FALSE)
Stroma1 <- RunUMAP(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- RunTSNE(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- FindNeighbors(Stroma1, reduction = "pca", dims = 1:20)

p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch1UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Stroma1) <- Stroma1$Genotype
p <- DimPlot(Stroma1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch1GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch1UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Stroma1, features = c("TTR","HGF","PCOLCE","DCN")  ,  reduction = "umap" )
ggsave(filename=paste(saveext,"/DimRed/StromaBatch1UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)


B1ListA <- colnames(Stroma1)[which(Stroma1$Genotype%in%c("Embryonic_G1"))] #
B1ListB <- colnames(Stroma1)[which(Stroma1$Genotype%in%c("Epithelial_G4"))] #


Stroma1 <- subset(Temp1,idents="Batch3")

Idents(Stroma1) <- Stroma1$HarmAno
Stroma1 <- FindVariableFeatures(Stroma1, selection.method = "vst", nfeatures = 20000)
Stroma1 <- ScaleData(Stroma1, verbose = FALSE)
Stroma1 <- RunPCA(Stroma1, npcs = 20, verbose = FALSE)
Stroma1 <- RunUMAP(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- RunTSNE(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- FindNeighbors(Stroma1, reduction = "pca", dims = 1:20)

p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch3UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Stroma1) <- Stroma1$Genotype
p <- DimPlot(Stroma1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch3GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch3UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Stroma1, features = c("TTR","HGF","PCOLCE","DCN")  ,  reduction = "umap" )
ggsave(filename=paste(saveext,"/DimRed/StromaBatch3UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)

B3ListB <- colnames(Epithel1)[which(Epithel1$Genotype%in%c("Epithelial_G3","Epithelial_G4","Epithelial_G6","Epithelial_G2"))]



Stroma1 <- subset(Temp1,idents="Batch6")

Idents(Stroma1) <- Stroma1$HarmAno
Stroma1 <- FindVariableFeatures(Stroma1, selection.method = "vst", nfeatures = 20000)
Stroma1 <- ScaleData(Stroma1, verbose = FALSE)
Stroma1 <- RunPCA(Stroma1, npcs = 20, verbose = FALSE)
Stroma1 <- RunUMAP(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- RunTSNE(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- FindNeighbors(Stroma1, reduction = "pca", dims = 1:20)

p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch3UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Stroma1) <- Stroma1$Genotype
p <- DimPlot(Stroma1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch6GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch6UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Stroma1, features = c("TTR","HGF","PCOLCE","DCN")  ,  reduction = "umap" )
ggsave(filename=paste(saveext,"/DimRed/StromaBatch6UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)

Stroma1 <- FindClusters(Stroma1,resolution = 0.3)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch6UMAP_Cl",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE,p)


B6ListA <- colnames(Stroma1)[which(Stroma1$Genotype%in%c("Embryonic_G4","NAssigned"))]


Stroma1 <- subset(Temp1,idents="Batch7")

Idents(Stroma1) <- Stroma1$HarmAno
Stroma1 <- FindVariableFeatures(Stroma1, selection.method = "vst", nfeatures = 20000)
Stroma1 <- ScaleData(Stroma1, verbose = FALSE)
Stroma1 <- RunPCA(Stroma1, npcs = 20, verbose = FALSE)
Stroma1 <- RunUMAP(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- RunTSNE(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- FindNeighbors(Stroma1, reduction = "pca", dims = 1:20)

p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch7UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Stroma1) <- Stroma1$Genotype
p <- DimPlot(Stroma1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch7GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch7UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Stroma1, features = c("TTR","HGF","PCOLCE","DCN")  ,  reduction = "umap" )
ggsave(filename=paste(saveext,"/DimRed/StromaBatch7UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)

Stroma1 <- FindClusters(Stroma1,resolution = 0.3)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch7UMAP_Cl",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE,p)


B7ListA <- colnames(Stroma1)[which(Stroma1$Genotype%in%c("Embryonic_G2"))] #


Stroma1 <- subset(Temp1,idents="Batch8")

Idents(Stroma1) <- Stroma1$HarmAno
Stroma1 <- FindVariableFeatures(Stroma1, selection.method = "vst", nfeatures = 20000)
Stroma1 <- ScaleData(Stroma1, verbose = FALSE)
Stroma1 <- RunPCA(Stroma1, npcs = 20, verbose = FALSE)
Stroma1 <- RunUMAP(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- RunTSNE(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- FindNeighbors(Stroma1, reduction = "pca", dims = 1:20)

p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch8UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Stroma1) <- Stroma1$Genotype
p <- DimPlot(Stroma1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch8GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch8UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Stroma1, features = c("TTR","HGF","PCOLCE","DCN")  ,  reduction = "umap" )
ggsave(filename=paste(saveext,"/DimRed/StromaBatch8UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)

Stroma1 <- FindClusters(Stroma1,resolution = 0.3)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch8UMAP_Cl",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE,p)


B8ListA <- colnames(Stroma1)[which(Stroma1$Genotype%in%c("Embryoic_G3"))] #






Stroma1 <- subset(Temp1,idents="10X Ours")

Idents(Stroma1) <- Stroma1$HarmAno
Stroma1 <- FindVariableFeatures(Stroma1, selection.method = "vst", nfeatures = 20000)
Stroma1 <- ScaleData(Stroma1, verbose = FALSE)
Stroma1 <- RunPCA(Stroma1, npcs = 20, verbose = FALSE)
Stroma1 <- RunUMAP(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- RunTSNE(Stroma1, reduction = "pca", dims = 1:20)
Stroma1 <- FindNeighbors(Stroma1, reduction = "pca", dims = 1:20)

p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch2UMAP",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
Idents(Stroma1) <- Stroma1$Genotype
p <- DimPlot(Stroma1, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch2GT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) #+xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch2UMAPGT",".pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(Stroma1, features = c("TTR","HGF","PCOLCE","DCN")  ,  reduction = "umap" )
ggsave(filename=paste(saveext,"/DimRed/StromaBatch2UMAP_Markers",".pdf",sep=""),width = 20, height = 20, limitsize = FALSE,p)

Stroma1 <- FindClusters(Stroma1,resolution = 0.3)
p <- DimPlot(Stroma1, reduction = "umap", pt.size = 2, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/StromaBatch2UMAP_Cl",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE,p)


B2ListA <- colnames(Stroma1)[which(Stroma1$Genotype%in%c("Embryonic_G1"))] #




Idents(Matteo,cells=WhichCells(Matteo,idents="Unknown")) <- "Epithelial"

#Epithelial cells with a embryonic GT

Idents(Matteo,cells=B1ListA) <- "StromaEmbryonicGT"
Idents(Matteo,cells=B2ListA) <- "StromaEmbryonicGT"

#Idents(Matteo,cells=B3ListA) <- "StromaEmbryonicGT"
Idents(Matteo,cells=B6ListA) <- "StromaEmbryonicGT"
Idents(Matteo,cells=B7ListA) <- "StromaEmbryonicGT"
Idents(Matteo,cells=B8ListA) <- "StromaEmbryonicGT"



p<-DimPlot(Matteo, pt.size = 4, label.size = 2, split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_extraAno",".pdf",sep=""),width = 60, height = 20, limitsize = FALSE,p)
saveRDS(Matteo,paste(saveext,"/MatteoAno_withOCT4AndStromaano.rds",sep=""))


OldAno <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno3.rds")


