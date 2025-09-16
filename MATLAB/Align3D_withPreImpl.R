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

marmoset <- readRDS("/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping4/AlignMatteo_withMarmosetwCynoMarmHumanSS.rds")

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
