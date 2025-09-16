library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("Matrix")
library(e1071)
library(pracma)
set.seed(1)
library(pheatmap)
#Same folder as before

saveext = "./FinalAlignB/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))

#Cell anotation and colours
cols <- c("#0767DA",
"#BF9600",
"#E68600",
"#E6B500",
"#F04C04",
"#D74404",
"#D74404",
"#BF3C04",
"#00BFBF",
"#0c9cf5",
"#0767DA",
"#0233BF",
"#0767DA",
"#877BD6",
"#7368CF",
"#5F54C7",
"#5F54C7",
"#5F54C7",
"#BF9600",
"#E68600",
"#AA8600",
"#921FE6",
"#8717D1",
"#7C10BB",
"#7108A6",
"#7C10BB",
"#877BD6",
"#6355B5",
"#3E2E94",
"#1A0873",
"#3E2E94",
"#BF9600",
"#BF9600",
"#BF0489",
"#00E6E6",
"#BF9600",
"#BF9600",
"#BF9600",
"#BF9600",
"#BF9600",
"#BF9600",
"#E68600",
"#E68600",
"#E68600",
"#E68600",
"#E68600",
"#E68600",
"#0767DA",
"#0767DA",
"#0767DA",
"#0767DA",
"#0767DA",
"#0767DA",
"lightgrey",
"#00BFBF",
"red",
"#5F54C7")


cType <-c("Embryonic",
"Stromal",
"Epithelial",
"Hyp_d9",
"Hyp_d11",
"Hyp_d12",
"VE",
"VE_d14",
"EmDisc_d9",
"EmDisc_d11",
"EmDisc_d12",
"EmDisc_d14",
"EmDisc",
"Am_d9",
"Am_d11",
"Am_d12",
"Am_d14",
"Am",
"Stromal",
"Uncilliated epithelial",
"Cilliated epithelial",
"CTB_d9",
"CTB_d11",
"CTB_d12",
"CTB_d14",
"CTB",
"STB_d9",
"STB_d11",
"STB_d12",
"STB_d14",
"STB",
"EVT_d14",
"EVT",
"Tr",
"ICM",
"Stromal_G1",
"Stromal_G2",
"Stromal_G3",
"Stromal_G4",
"Stromal_G5",
"Stromal_G6",
"Epithelial_G1",
"Epithelial_G2",
"Epithelial_G3",
"Epithelial_G4",
"Epithelial_G5",
"Epithelial_G6",
"Embryonic_G1",
"Embryonic_G2",
"Embryonic_G3",
"Embryonic_G4",
"Embryonic_G5",
"Embryonic_G6",
"NGT",
"Epi",
"Cont",
"Am/EmDisc_d14") 



#Load in the previous run
mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))
mammal.combined$iid <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$Dataset
mammal.combined <- subset(mammal.combined,idents=c("10X Ours"))
mito.genes <- grep(pattern = "^MT-", x = rownames(mammal.combined@assays[["RNA"]]), value = TRUE)
percent.mito2 <- Matrix::colSums(mammal.combined@assays[["RNA"]][mito.genes, ])/Matrix::colSums(mammal.combined@assays[["RNA"]])
mammal.combined <- AddMetaData(object = mammal.combined, metadata = percent.mito2, col.name = "percent.mito2")

FeaturePlot(mammal.combined, reduction = "pca", features = "nFeature_RNA", pt.size = 2, cols = c("lightgrey", "black")) # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_nFeature_RNA_nosplit.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, reduction = "pca", features = "nCount_RNA", pt.size = 2, cols = c("lightgrey", "black")) # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_nCount_RNA_nosplit.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, reduction = "pca", features = "percent.mito2", pt.size = 2, cols = c("lightgrey", "black")) # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_mito_nosplit.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

VlnPlot(mammal.combined,features=c("nCount_RNA"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "_nCount_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)
VlnPlot(mammal.combined,features=c("nFeature_RNA"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "nFeature_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)

VlnPlot(mammal.combined,features=c("percent.mito2"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "mito_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)


list1<-c("GATA4",
"PDGFRA",
"CER1",
"SOX17",
"APOA1",
"APOB",
"VIM",
"BST2",
"HGF",
"HAND2",
"POU5F1",
"SOX15",
"SFRP1",
"DNMT3B",
"PDGFA",
"TFAP2A",
"TFAP2C",
"WNT6",
"VTCN1",
"GATA2",
"GATA3",
"CGA",
"HLA-G",
"JAM3",
"DCN",
"MSLN")


mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))
mammal.combined2$iid <- Idents(mammal.combined)

Idents(mammal.combined2) <- mammal.combined2$Dataset
ourcells <- WhichCells(mammal.combined2,idents="10X Ours")

dat1 <- subset(mammal.combined2,idents="10X Ours")
Idents(dat1) <- dat1$iid
dat1 <- subset(dat1,idents=c("ExMes","EmD/Hyp"))
Idents(dat1) <- dat1$iid

Idents(mammal.combined2) <- mammal.combined2$iid
ourcells1 <- WhichCells(mammal.combined2,idents=c("ExMes","EmD/Hyp") )
ourcells2 <- intersect(ourcells1, ourcells )

#dat1 <- mammal.combined2
D <- GetAssayData(dat1,assay="RNA")

saveRDS(D,file="subsetdata.rds")
D <- as.data.frame(D)[list1,]

mat_breaks <- seq(-2, 2, length.out = 20)

redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes_Tb.pdf",sep=""),scale="row",width=20,height=10)

#pheatmap(log2(D+1),color =  redblue1(20),kmeans_k=5, breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesKmea.pdf",sep=""),scale="row",width=10,height=10)

res5 <- pheatmap((D),color =  redblue1(20), cutree_cols = 25,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes_Tb.pdf",sep=""),scale="row",width=20,height=10)#
#clust <- cbind(log2(D+1), cluster = cutree(res$tree_col, k = 25))
#saveRDS(clust, file = "Cluster.rds")

res1 <- pheatmap((D),color =  redblue1(20), cutree_cols = 5,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes5_Tb.pdf",sep=""),scale="row",width=20,height=10)
res2 <- pheatmap((D),color =  redblue1(20), cutree_cols = 10,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes10_Tb.pdf",sep=""),scale="row",width=20,height=10)
res3 <- pheatmap((D),color =  redblue1(20), cutree_cols = 15,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes15_Tb.pdf",sep=""),scale="row",width=20,height=10)
res4 <- pheatmap((D),color =  redblue1(20), cutree_cols = 20,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes20_Tb.pdf",sep=""),scale="row",width=20,height=10)
res6 <- pheatmap((D),color =  redblue1(20), cutree_cols = 30,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes30_Tb.pdf",sep=""),scale="row",width=20,height=10)
res7 <- pheatmap((D),color =  redblue1(20), cutree_cols = 35,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes35_Tb.pdf",sep=""),scale="row",width=20,height=10)


#clust <- cbind(log2(D+1), cluster = cutree(res$tree_col, k = 25))
saveRDS(D, file = "ClusterDataTb.rds")
saveRDS(res1, file = "Cluster1Tb.rds")
saveRDS(res2, file = "Cluster2Tb.rds")
saveRDS(res3, file = "Cluster3Tb.rds")
saveRDS(res4, file = "Cluster4Tb.rds")
saveRDS(res5, file = "Cluster5Tb.rds")
saveRDS(res6, file = "Cluster6Tb.rds")
saveRDS(res7, file = "Cluster7Tb.rds")

ct1 <- sort(cutree(res1$tree_col, k=5))
ct2 <- sort(cutree(res2$tree_col, k=10))
ct3 <- sort(cutree(res3$tree_col, k=15))
ct4 <- sort(cutree(res4$tree_col, k=20))
ct5 <- sort(cutree(res5$tree_col, k=25))
ct6 <- sort(cutree(res6$tree_col, k=30))
ct7 <- sort(cutree(res7$tree_col, k=35))

write.table(ct1,file="Cluster1Tb.csv",quote = FALSE)
write.table(ct2,file="Cluster2Tb.csv",quote = FALSE)
write.table(ct3,file="Cluster3Tb.csv",quote = FALSE)
write.table(ct4,file="Cluster4Tb.csv",quote = FALSE)
write.table(ct5,file="Cluster5Tb.csv",quote = FALSE)
write.table(ct6,file="Cluster6Tb.csv",quote = FALSE)
write.table(ct7,file="Cluster7Tb.csv",quote = FALSE)


cl1 <- names(which(ct2==1))
cl2 <- names(which(ct2==2))
cl3 <- names(which(ct2==3))
cl4 <- names(which(ct2==4))
cl5 <- names(which(ct2==5))
cl6 <- names(which(ct2==6))
cl7 <- names(which(ct2==7))
cl8 <- names(which(ct2==8))
cl9 <- names(which(ct2==9))
cl10 <- names(which(ct2==10))

saveRDS(dat1,file="dat1Tb.rds")

mammal.combined4B <- dat1
Idents(mammal.combined4B,cells=cl1) <- "Clu1"
Idents(mammal.combined4B,cells=cl2) <- "Clu2"
Idents(mammal.combined4B,cells=cl3) <- "Clu3"
Idents(mammal.combined4B,cells=cl4) <- "Clu4"
Idents(mammal.combined4B,cells=cl5) <- "Clu5"
Idents(mammal.combined4B,cells=cl6) <- "Clu6"
Idents(mammal.combined4B,cells=cl7) <- "Clu7"
Idents(mammal.combined4B,cells=cl8) <- "Clu8"
Idents(mammal.combined4B,cells=cl9) <- "Clu9"
Idents(mammal.combined4B,cells=cl10) <- "Clu10"

p1 <- DimPlot(mammal.combined4B,  cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_fullano",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)




#mammal.combined4B <- subset(mammal.combined4B,idents=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
AvE <- AverageExpression(mammal.combined4B)
saveRDS(AvE,file="AvExp.rds")
AvE <- AvE$RNA
AvE <- AvE[list1,]

PB <- pheatmap(log2(AvE+1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesTbPB.pdf",sep=""),scale="row",width=10,height=10)



asdadsadad


cluster1cells <- cl1


res1 <- pheatmap((D[,cl1]),color =  redblue1(20), cutree_cols = 5,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes5_Tbsse.pdf",sep=""),scale="row",width=20,height=10)
res2 <- pheatmap((D[,cl1]),color =  redblue1(20), cutree_cols = 10,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes10_Tbsse.pdf",sep=""),scale="row",width=20,height=10)
res3 <- pheatmap((D[,cl1]),color =  redblue1(20), cutree_cols = 15,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes15_Tbsse.pdf",sep=""),scale="row",width=20,height=10)
res4 <- pheatmap((D[,cl1]),color =  redblue1(20), cutree_cols = 20,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes20_Tbsse.pdf",sep=""),scale="row",width=20,height=10)
res5 <- pheatmap((D[,cl1]),color =  redblue1(20), cutree_cols = 25,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes_Tbsse.pdf",sep=""),scale="row",width=20,height=10)#
res6 <- pheatmap((D[,cl1]),color =  redblue1(20), cutree_cols = 30,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes30_Tbsse.pdf",sep=""),scale="row",width=20,height=10)
res7 <- pheatmap((D[,cl1]),color =  redblue1(20), cutree_cols = 35,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes35_Tbsse.pdf",sep=""),scale="row",width=20,height=10)



ct6 <- sort(cutree(res6$tree_col, k=30))

cl1 <- names(which(ct6==1))
cl2 <- names(which(ct6==2))
cl3 <- names(which(ct6==3))
cl4 <- names(which(ct6==4))
cl5 <- names(which(ct6==5))
cl6 <- names(which(ct6==6))
cl7 <- names(which(ct6==7))
cl8 <- names(which(ct6==8))
cl9 <- names(which(ct6==9))
cl10 <- names(which(ct6==10))
cl11 <- names(which(ct6==11))
cl12 <- names(which(ct6==12))
cl13 <- names(which(ct6==13))
cl14 <- names(which(ct6==14))
cl15 <- names(which(ct6==15))
cl16 <- names(which(ct6==16))
cl17 <- names(which(ct6==17))
cl18 <- names(which(ct6==18))
cl19 <- names(which(ct6==19))
cl20 <- names(which(ct6==20))
cl21 <- names(which(ct6==21))
cl22 <- names(which(ct6==22))
cl23 <- names(which(ct6==23))
cl24 <- names(which(ct6==24))
cl25 <- names(which(ct6==25))
cl26 <- names(which(ct6==26))
cl27 <- names(which(ct6==27))
cl28 <- names(which(ct6==28))
cl29 <- names(which(ct6==29))
cl30 <- names(which(ct6==30))


#mammal.combined4B <- dat1
Idents(mammal.combined4B,cells=cl1) <- "Clu1B"
Idents(mammal.combined4B,cells=cl2) <- "Clu2B"
Idents(mammal.combined4B,cells=cl3) <- "Clu3B"
Idents(mammal.combined4B,cells=cl4) <- "Clu4B"
Idents(mammal.combined4B,cells=cl5) <- "Clu5B"
Idents(mammal.combined4B,cells=cl6) <- "Clu6B"
Idents(mammal.combined4B,cells=cl7) <- "Clu7B"
Idents(mammal.combined4B,cells=cl8) <- "Clu8B"
Idents(mammal.combined4B,cells=cl9) <- "Clu9B"
Idents(mammal.combined4B,cells=cl10) <- "Clu10B"
Idents(mammal.combined4B,cells=cl11) <- "Clu11B"
Idents(mammal.combined4B,cells=cl12) <- "Clu12B"
Idents(mammal.combined4B,cells=cl13) <- "Clu13B"
Idents(mammal.combined4B,cells=cl14) <- "Clu14B"
Idents(mammal.combined4B,cells=cl15) <- "Clu15B"
Idents(mammal.combined4B,cells=cl16) <- "Clu16B"
Idents(mammal.combined4B,cells=cl17) <- "Clu17B"
Idents(mammal.combined4B,cells=cl18) <- "Clu18B"
Idents(mammal.combined4B,cells=cl19) <- "Clu19B"
Idents(mammal.combined4B,cells=cl20) <- "Clu20B"
Idents(mammal.combined4B,cells=cl21) <- "Clu21B"
Idents(mammal.combined4B,cells=cl22) <- "Clu22B"
Idents(mammal.combined4B,cells=cl23) <- "Clu23B"
Idents(mammal.combined4B,cells=cl24) <- "Clu24B"
Idents(mammal.combined4B,cells=cl25) <- "Clu25B"
Idents(mammal.combined4B,cells=cl26) <- "Clu26B"
Idents(mammal.combined4B,cells=cl27) <- "Clu27B"
Idents(mammal.combined4B,cells=cl28) <- "Clu28B"
Idents(mammal.combined4B,cells=cl29) <- "Clu29B"
Idents(mammal.combined4B,cells=cl30) <- "Clu30B"


Idents(mammal.combined4B,cells=AEVT) <- "CluEVT"


AvE <- AverageExpression(mammal.combined4B)
#xsaveRDS(AvE,file="AvExp.rds")
AvE <- AvE$RNA
AvE <- AvE[list1,]

PB <- pheatmap(log2(AvE+1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesTbPBsse.pdf",sep=""),scale="row",width=15,height=10)




saveRDS(res1, file = "Cluster1Tbsse.rds")
saveRDS(res2, file = "Cluster2Tbsse.rds")
saveRDS(res3, file = "Cluster3Tbsse.rds")
saveRDS(res4, file = "Cluster4Tbsse.rds")
saveRDS(res5, file = "Cluster5Tbsse.rds")
saveRDS(res6, file = "Cluster6Tbsse.rds")
saveRDS(res7, file = "Cluster7Tbsse.rds")



sadadad
mammal.combined2test <- mammal.combined2

Idents(mammal.combined2test,cells=cl1) <- "VE_d14"
Idents(mammal.combined2test,cells=cl2) <- "Am/EmDisc_d14"
Idents(mammal.combined2test,cells=cl3) <- "Am/EmDisc_d14"
Idents(mammal.combined2test,cells=cl4) <- "Cont"
Idents(mammal.combined2test,cells=cl5) <- "VE_d14"
Idents(mammal.combined2test,cells=cl6) <- "VE_d14"
Idents(mammal.combined2test,cells=cl7) <- "Cont"
Idents(mammal.combined2test,cells=cl8) <- "Cont"
Idents(mammal.combined2test,cells=cl9) <- "Am_d14"
Idents(mammal.combined2test,cells=cl10) <- "Cont"
Idents(mammal.combined2test,cells=cl11) <- "EmDisc_d14"
Idents(mammal.combined2test,cells=cl12) <- "EmDisc_d14"
Idents(mammal.combined2test,cells=cl13) <- "Cont"
Idents(mammal.combined2test,cells=cl14) <- "EmDisc_d14"
Idents(mammal.combined2test,cells=cl15) <- "CTB_d14"
Idents(mammal.combined2test,cells=cl16) <- "CTB_d14"
Idents(mammal.combined2test,cells=cl17) <- "VE_d14"
Idents(mammal.combined2test,cells=cl18) <- "VE_d14"
Idents(mammal.combined2test,cells=cl19) <- "VE_d14"
Idents(mammal.combined2test,cells=cl20) <- "VE_d14"

colind <- integer( length( levels(Idents(mammal.combined2test)) )  )
for (i in 1:length( levels(Idents(mammal.combined2test)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2test))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined2test,  cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_fullano",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



#mammal.combined4B <- subset(mammal.combined4B,idents=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
AvE <- AverageExpression(mammal.combined4B)
saveRDS(AvE,file="AvExp.rds")
AvE <- AvE$RNA
AvE <- AvE[list1,]

PB <- pheatmap(log2(AvE+1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesPB.pdf",sep=""),scale="row",width=10,height=10)


mammal.combined4B <- subset(mammal.combined4B,idents=c("Clu1","Clu2","Clu3","Clu4","Clu5","Clu6","Clu7","Clu8","Clu9","Clu10","Clu11","Clu12","Clu13","Clu14","Clu15","Clu16","Clu17","Clu18","Clu19","Clu20"))



PB <- pheatmap(log2(AvE+1),color =  redblue1(20),  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesPB2.pdf",sep=""),scale="row",width=10,height=10)


mammal.combined3 <- mammal.combined2
Concells <- WhichCells(mammal.combined3,expression=MSLN>0.2)
Idents(mammal.combined3,cells=intersect(Concells,ourcells)) <- "Cont"
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_Cont1","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined3,idents="Cont",invert=TRUE)


mammal.combined3 <- mammal.combined2
Concells <- WhichCells(mammal.combined3,expression=DCN>0.1)
Idents(mammal.combined3,cells=intersect(Concells,ourcells)) <- "Cont"
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_Cont2","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

mammal.combined2 <- subset(mammal.combined3,idents="Cont",invert=TRUE)


mammal.combined3 <- mammal.combined2
VEcells3 <- WhichCells(mammal.combined3,expression=APOA1>0.2)
Idents(mammal.combined3,cells=intersect(VEcells3,ourcells2)) <- "Cont"

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_VE3","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)


mammal.combined3 <- mammal.combined2
VEcells <- WhichCells(mammal.combined3,expression=CER1>0.2)
Idents(mammal.combined3,cells=intersect(VEcells,ourcells2) ) <- "Cont"

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined2, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_VE1","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)


mammal.combined3 <- mammal.combined2
EmDcells1 <- WhichCells(mammal.combined3,expression=POU5F1>0.2)
Idents(mammal.combined3,cells=intersect(EmDcells1,ourcells2) ) <- "Cont"

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_EmD1","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



mammal.combined3 <- mammal.combined2
VEcells2 <- WhichCells(mammal.combined3,expression=GATA4>0.2)
Idents(mammal.combined3,cells=intersect(VEcells2,ourcells2)) <- "Cont"

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_VE2","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



mammal.combined3 <- mammal.combined2
VEcells3 <- WhichCells(mammal.combined3,expression=APOA1>0.2)
Idents(mammal.combined3,cells=intersect(VEcells3,ourcells2)) <- "Cont"

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_VE3","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)




Idents(mammal.combined2) <- mammal.combined2$Cl05

p1 <- DimPlot(mammal.combined2, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_Cl",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)





sadsdsadasdasd



p1 <- DimPlot(mammal.combined2,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

#sdsadsad

colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined2, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

Idents(mammal.combined2) <- mammal.combined2$Cl05

p1 <- DimPlot(mammal.combined2, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_Cl",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



saveRDS(mammal.combined2,file = paste(saveext,"Seurat_combined_filtered.rds",sep=""))

DefaultAssay(mammal.combined2) <- "RNA"


Idents(mammal.combined2) <- mammal.combined2$Dataset
mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")

FeaturePlot(mammal.combined2,feature="LGALS9",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_LGALS9_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="HAND2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_HAND2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CD44",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_CD44_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="SFRP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SFRP1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="RARRES2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_RARRES2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="GPR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GPR1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="C5AR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_C5AR1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="NRG1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_NRG1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="ERBB4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_ERBB4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="EGFR",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_EGFR_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="COPA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_COPA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


#FeaturePlot(mammal.combined2,feature="COPA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
#ggsave(filename=paste(saveext,"FinalPCA_COPA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
#

FeaturePlot(mammal.combined2,feature="ROR2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_ROR2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="FOXJ1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_FOXJ1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="TFAP2C",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_TFAP2C_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="PDGFA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_PDGFA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="MUC4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_MUC4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="MSLN",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_MSLN_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FGF2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_FGF2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="VTCN1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_VTCN1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="TFAP2A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_TFAP2A_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="ISL1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_ISL1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="POU5F1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_POU5F1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="OTX2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_OTX2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="FSTL3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_FSTL3_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="EOMES",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_EOMES_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="S100P",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_S100P_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SDC1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SDC1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="PGF",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_PGF_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DKK1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DKK1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="SFRP2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SFRP2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="TDGF1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_TDGF1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DPPA5",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DPPA5_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined2,feature="TIMP3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_TIMP2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="ID2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_ID2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="NRP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_NRP1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DMKN",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DMKN_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined2,feature="SOX15",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SOX15_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="PDGFRA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_PDGFRA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="BAMBI",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_BAMBI_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DIO2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DIO2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="HLA-G",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_HLA-G_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="JAM3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_JAM3_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="APOA1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_APOA1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="GATA4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GATA4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="CER1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_CER1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="CGA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_CGA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="GATA3",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GATA3_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="GATA3",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GATA3_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="DIO2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DIO2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="BAMBI",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_BAMBI_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FOXJ1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_FOXJ1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DCN",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DCN_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="WNT7A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_WNT7A_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="LGR5",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_LGR5_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="SOX9",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SOX9_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined2,feature="PAEP",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_PAEP_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SPP1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SPP1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="SCGB2A2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SCGB2A2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_WNT6_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "umap",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalUMAP_WNT6_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="GABRP",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GABRP_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="DNMT3B",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DNMT3B_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="SFRP2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SFRP2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="STAT3",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_STAT3_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="HAND2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_HAND2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="TBX4",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_TBX4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


#Idents(mammal.combined2) <- mammal.combined2$Cells 
mammal.combined2a <- mammal.combined2 #subset(mammal.combined2,idents="Embryonic")

FeatureScatter(mammal.combined2a,feature1="CGA",feature2="GATA4")
ggsave(filename=paste(saveext,"Scatter_CGA_GATA4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2a,feature1="GATA2",feature2="GATA4")
ggsave(filename=paste(saveext,"Scatter_GATA2_GATA4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2a,feature1="CGA",feature2="PDGFA")
ggsave(filename=paste(saveext,"Scatter_CGA_PDGFA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2a,feature1="GATA4",feature2="PDGFA")
ggsave(filename=paste(saveext,"Scatter_GATA4_PDGFA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeatureScatter(mammal.combined2a,feature1="POU5F1",feature2="PDGFA")
ggsave(filename=paste(saveext,"Scatter_POU5F1_PDGFA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2a,feature1="POU5F1",feature2="GATA4")
ggsave(filename=paste(saveext,"Scatter_POU5F1_GATA4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2a,feature1="CGA",feature2="POU5F1")
ggsave(filename=paste(saveext,"Scatter_CGA_POU5F1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2a,feature1="GATA2",feature2="POU5F1")
ggsave(filename=paste(saveext,"Scatter_POU5F1_GATA2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined2a,feature1="SOX15",feature2="POU5F1")
ggsave(filename=paste(saveext,"Scatter_POU5F1_SOX15_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)#
#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)



#require(tidyverse)
#X <- mammal.combined2$Dataset
#Y <- mammal.combined2$Cl15   
#Z <- as.character(Idents(mammal.combined2))

#Z2 <- mammal.combined2$Genotype

#Dat <- data.frame(x=Y[which(X=="SS2 Reference 1")],y=Z[which(X=="SS2 Reference 1")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table1","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

#Dat <- data.frame(x=Y[which(X=="10X Reference 1")],y=Z[which(X=="10X Reference 1")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


#Dat <- data.frame(x=Y[which(X=="10X Ours")],y=Z2[which(X=="10X Ours")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table3","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


saddsadsdasds


mammal.combined3 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))
Idents(mammal.combined3) <- mammal.combined3$Dataset
mammal.combined3 <- subset(mammal.combined3,idents=c("10X Ours"))


uID <- as.character(mammal.combined3$Genotype)
uID[which(is.na(uID)==1)] <- "NGT"
Idents(mammal.combined3) <- uID

Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("NAssigned"))) <- "NGT"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("G1"))) <- "Stromal_G1"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("G3"))) <- "Epithelial_G3"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("G4"))) <- "Epithelial_G4"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("Embryoni_G3"))) <- "Embryonic_G3"

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, split.by="ID3", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_Genotype_ano","boundary.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE)

mammal.combined3 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))

Idents(mammal.combined3) <- mammal.combined3$Cl05

cType <- c("4",
"8",
"12",
"14",
"16",
"19",
"20",
"21",
"0",
"2",
"3",
"5",
"6",
"7",
"10",
"13",
"15",
"18",
"1",
"9",
"11",
"17")

cols <- c("#00BFBF",
"#00ABBF",
"#0197BF",
"#0183BF",
"#016FBF",
"#015BBF",
"#0247BF",
"#0233BF",
"#FAA734",
"#F8A32E",
"#F6A028",
"#F39C23",
"#F1981D",
"#EF9517",
"#ED9111",
"#EA8D0C",
"#E88A06",
"#E68600",
"#D9AD0D",
"#C49C0A",
"#AE8A06",
"#997903")

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_Cl05_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE, limitsize = FALSE)






mammal.combined3 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))

uID <- as.character(Idents(mammal.combined3))
uID[which(uID%in%c("CTb_CS4","EPI_CS4","TB_CS3"))] <- "Tr"
uID[which(uID%in%c("CTb_CS5A/B","CTb_CS5C","CTb_CS6"))] <- "CTB"
uID[which(uID%in%c("EmDisc_CS5A/B","EmDisc_CS5C","EmDisc_CS6","EmDiscPS_CS6"))] <- "EmDisc"
uID[which(uID%in%c("EVTb_CS5C","EVTb_CS6"))] <- "EVT"
uID[which(uID%in%c("STb_CS5A/B","STb_CS5C","STb_CS6"))] <- "STB"
uID[which(uID%in%c("VE_CS4","VE_CS5","VE_CS5A/B","VE_CS5C","SYS_CS6"))] <- "VE"
uID[which(uID%in%c("ICM_CS3"))] <- "ICM"
uID[which(uID%in%c("Unciliated epithelia 1","Unciliated epithelia 2","Epith"))] <- "Unciliated epithelia"
uID[which(uID%in%c("Cil"))] <- "Ciliated epithelia"
uID[which(uID%in%c("Ciliated"))] <- "Ciliated epithelia"
uID[which(uID%in%c("Stroma"))] <- "Stroma fibroblasts"
uID[which(Idents(mammal.combined3)%in%c("STB_d9"))] <- "CTB_d9"
uID[which(Idents(mammal.combined3)%in%c("STB_d11"))] <- "CTB_d11"
uID[which(Idents(mammal.combined3)%in%c("STB_d12"))] <- "CTB_d12"
uID[which(Idents(mammal.combined3)%in%c("CTB_d9"))] <- "STB_d9"
uID[which(Idents(mammal.combined3)%in%c("CTB_d11"))] <- "STB_d11"
uID[which(Idents(mammal.combined3)%in%c("CTB_d12"))] <- "STB_d12"
#uID[which(uID%in%c("Hyp_d9","Hyp_d11","Hyp_d12"))] <- "VE"
#uID[which(uID%in%c("EmDisc_d9","EmDisc_d11","EmDisc_d12"))] <- "EmDisc"
uID[which(Idents(mammal.combined3)%in%c("CTB") & mammal.combined3$Dataset==c("10X Ours"))] <- "STB"
uID[which(Idents(mammal.combined3)%in%c("STB") & mammal.combined3$Dataset==c("10X Ours"))] <- "CTB"
uID[which(mammal.combined3$Cl15%in%c(8,18,19,23,4,12,15,20,22,28,29,31,32,37,38) & mammal.combined3$Genotype%in%c("Epithelial_G1","Epithelial_G2","Epithelial_G3","Epithelial_G4"))] <- "Tb_Epith_Mix"
uID[which(mammal.combined3$Cl15%in%c(8,18,19,23,4,12,15,20,22,28,29,31,32,37,38) & mammal.combined3$Genotype%in%c("Stromal_G1","Stromal_G2","Stromal_G3","Stromal_G4","Stromal_G5","Stromal_G6"))] <- "Tb_Stroma_Mix"
Idents(mammal.combined3) <- uID

Idents(mammal.combined3,cells=AmList) <- "Am"
Idents(mammal.combined3,cells=CTBList) <- "CTB"
Idents(mammal.combined3,cells=eEVTList) <- "CTB"
Idents(mammal.combined3,cells=eSTBList) <- "STB"
Idents(mammal.combined3,cells=EVTList) <- "EVT"
Idents(mammal.combined3,cells=ICMList) <- "ICM"
Idents(mammal.combined3,cells=PostEpiAME) <- "Am"
Idents(mammal.combined3,cells=PostEPIE1) <- "EmDisc"
Idents(mammal.combined3,cells=PostEPIE2) <- "EmDisc"

Idents(mammal.combined3,cells=PostEpiGast) <- "EmDisc"
Idents(mammal.combined3,cells=PrEn) <- "VE"
Idents(mammal.combined3,cells=Epi) <- "Epi"
Idents(mammal.combined3,cells=EpiPrEn) <- "Epi"
Idents(mammal.combined3,cells=STB) <- "STB"
Idents(mammal.combined3,cells=Tr) <- "Tr"

Idents(mammal.combined3,cells=pICM) <- "Unknown"
Idents(mammal.combined3,cells=pEVT) <- "Unknown"
Idents(mammal.combined3,cells=pSTB) <- "STB"
Idents(mammal.combined3,cells=pVE) <- "VE"
Idents(mammal.combined3,cells=pEpi) <- "Epi"
Idents(mammal.combined3,cells=pPS) <- "Unknown"
Idents(mammal.combined3,cells=pEmD) <- "Unknown"

Idents(mammal.combined3,cells=WhichCells(mammal.combined3,
idents=c("CTB_d9","CTB_d11","CTB_d12","Hyp_d9","Hyp_d11","Hyp_d12","EmD/Hyp","EmDisc","Unkown_Emb","EmDisc_d9","EmDisc_d11","EmDisc_d12","STB_d9","STB_d11","STB_d12","Tr","STB","VE","ICM","EVT","Epi","STB","EVT","CTB","Am") )) <- "Embryonic"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("Unciliated epithelia","Ciliated epithelia") )) <- "Epithelial"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("Stroma fibroblasts","Stromal fibroblasts") )) <- "Stromal"

mammal.combined3 <- subset(mammal.combined3,idents="Unknown",invert=TRUE)


require(tidyverse)
X <- mammal.combined2$Dataset
X2 <- mammal.combined2$ID3
Y <- factor(mammal.combined2$Cl05,levels=c(4,8,12,14,16,19,20,15,5,7,21,0,2,3,6,10,13,18,1,9,11,17))
Z <- as.character(Idents(mammal.combined2))

Z2 <- mammal.combined2$Genotype
Z2[which(Z2=="NAssigned")] <- "NGT"

Z2[which(Z2=="G1")] <- "Stromal_G1"
Z2[which(Z2=="G3")] <- "Epithelial_G3"
Z2[which(Z2=="G4")] <- "Epithelial_G4"
Z2[which(Z2=="Embryoni_G3")] <- "Embryonic_G3"


Z2[which(Z2=="Epithelial_G1")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G2")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G3")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G4")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G5")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G6")] <- "Epithelail"

Z2[which(Z2=="Stromal_G1")] <- "Stromal"
Z2[which(Z2=="Stromal_G2")] <- "Stromal"
Z2[which(Z2=="Stromal_G3")] <- "Stromal"
Z2[which(Z2=="Stromal_G4")] <- "Stromal"
Z2[which(Z2=="Stromal_G5")] <- "Stromal"
Z2[which(Z2=="Stromal_G6")] <- "Stromal"

Z2[which(Z2=="Embryonic_G1")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G2")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G3")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G4")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G5")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G6")] <- "Embryonic"





Dat <- data.frame(x=Y[which(X=="SS2 Reference 1")],y=Z[which(X=="SS2 Reference 1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"Table1","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


Dat <- data.frame(x=Y[which(X=="10X Reference 1")],y=Z[which(X=="10X Reference 1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"Table2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)



#c("#0767DA","#E68600","#BF9600")

#Dat <- data.frame(x=Y[which(X=="10X Reference 1")],y=Z[which(X=="10X Reference 1")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


Dat <- data.frame(x=Y[which(X2=="C1")],y=Z2[which(X2=="C1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC1","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

Dat <- data.frame(x=Y[which(X2=="C2")],y=Z2[which(X2=="C2")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
Dat <- data.frame(x=Y[which(X2=="C3")],y=Z2[which(X2=="C3")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC3","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
Dat <- data.frame(x=Y[which(X2=="C4")],y=Z2[which(X2=="C4")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC4","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

Dat <- data.frame(x=Y[which(X2=="C6")],y=Z2[which(X2=="C6")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC6","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
Dat <- data.frame(x=Y[which(X2=="C7")],y=Z2[which(X2=="C7")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC7","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

Dat <- data.frame(x=Y[which(X2=="C8")],y=Z2[which(X2=="C8")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC8","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)



sadsdadsadadsadsdas
#Idents(SSD) <- uID

#uID[which(uID=="EPI_CS4")] <- "TB_CS4"

uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(13))] <- "EmDisc_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(16))] <- "Hyp_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(3,5,12,17))] <- "CTB_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(0,1,2,4,6,7,8,9,10,11,14,15,18,19,20))] <- "STB_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl15%in%c(5))] <- "EVT_d14"


Idents(mammal.combined1) <- uID #mammal.combined1$Tissue

saveRDS(mammal.combined1,file = paste(saveext,"Seurat_combined_justEmbryonic_ano.rds",sep=""))


colind <- integer( length( levels(Idents(mammal.combined1)) )  )
for (i in 1:length( levels(Idents(mammal.combined1)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined1))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined1, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_2_Embryoniconly_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_3_Embryoniconly_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(2,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_2_3_Embryoniconly_ano","boundary.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE)

p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_1_2Embryoniconly_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, split.by="ID3", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_1_2Embryoniconly_ano","boundary_byID3.pdf",sep=""),width = 120, height = 10, useDingbats = FALSE, limitsize = FALSE)


mammal.combined1$ano <- Idents(mammal.combined1)
Idents(mammal.combined1) <- mammal.combined1$Dataset

mammal.combined3 <- subset(mammal.combined1,idents=c("10X Ours"))
mammal.combined2 <- subset(mammal.combined1,idents=c("SS2 Reference 1"))
mammal.combined1 <- subset(mammal.combined1,idents=c("10X Reference 1"))


#sadsadsad

uID <- as.character(mammal.combined1$ano)
uID[which(uID%in%c("EmDisc_d9","EmDisc_d11","EmDisc_d12"))] <- "1) EmDisc"
uID[which(uID%in%c("CTB_d9","CTB_d11","CTB_d12","STB_d9","STB_d11","STB_d12","Hyp_d9","Hyp_d11","Hyp_d12"))] <- "2) Not EmDisc"
Idents(mammal.combined1) <- uID

#Grid of points for later
make.grid = function(x, n = 1000) {
  grange = apply(x, 2, range)
  x1 = seq(from = grange[1,1], to = grange[2,1], length = n)
  x2 = seq(from = grange[1,2], to = grange[2,2], length = n)
  expand.grid(X1 = x1, X2 = x2)
}

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = .5, scale = TRUE)
xgrid = make.grid(as.matrix(dat[,1:2]))
ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined1) <- mammal.combined1$ano
colind <- integer( length( levels(Idents(mammal.combined1)) ) )  
for (i in 1:length( levels(Idents(mammal.combined1)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined1))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined1, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Embryonic","boundary_v_maternal_embryoniconly.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


Idents(mammal.combined3) <- mammal.combined3$ano
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C8_Embryonic_v_maternal","boundary_embryoniconly.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined3$ano) 
uID[which(mammal.combined3$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined3$Cl05%in%c(10))] <- "Novel1"
uID[which(mammal.combined3$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined3$Cl05%in%c(19))] <- "Novel2"
uID[which(mammal.combined3$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined3$Cl05%in%c(20))] <- "Novel3"
Idents(mammal.combined3) <- uID

DefaultAssay(mammal.combined3) <- "RNA"
#Cl <- FindAllMarkers(mammal.combined3, test.use = "MAST")

saveRDS(uID,file="IDs.rds")
saveRDS(GetAssayData(mammal.combined3,assay="RNA"),file="RNAData.rds")



Idents(mammal.combined3) <- mammal.combined3$Genotype
p1 <- DimPlot(mammal.combined3,  pt.size = 4, reduction = "umap", label = TRUE, split.by="ID3", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_1_2Embryoniconly_ano","genotype_byID3.pdf",sep=""),width = 120, height = 10, useDingbats = FALSE, limitsize = FALSE)

mammal.combined3B <- subset(mammal.combined3,idents=c("Epithelial_G2","Epithelial_G4","Stromal_G1","Stromal_G3","Epithelial_G1","Epithelial_G6","Stromal_G4","Epithelial_G3","Stromal_G2"),invert=TRUE)

Idents(mammal.combined3B) <- mammal.combined3B$ano

colind <- integer( length( levels(Idents(mammal.combined3B)) )  )
for (i in 1:length( levels(Idents(mammal.combined3B)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3B))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined3B, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C8_Embryonic_testfilter","boundary_embryoniconly.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

saveRDS(mammal.combined3B,file=paste(saveext,"Seurat_combined_justEmbryonic_filtered.rds",sep=""))



saveRDS(mammal.combined3B$ano,file="IDs.rds")
saveRDS(GetAssayData(mammal.combined3B,assay="RNA"),file="RNAData.rds")

#write.table(as.data.frame(Cl),file="ClustersEmbryonicOnly.csv",quote = FALSE, sep = ",")

DefaultAssay(mammal.combined3) <- "integrated"
genelist <- rownames(GetAssayData(mammal.combined3))

Idents(mammal.combined3) <- mammal.combined3$Cl05
DefaultAssay(mammal.combined3) <- "RNA"
Ae <- AverageExpression(mammal.combined3)
Ae <- Ae$RNA
#Ae$gene <- rownames(Ae)

library(pheatmap)

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
mat_breaks <- seq(0, 6, length.out = 20)

C1 <- cor(  log2(Ae +1),log2(Ae +1) )
C2 <- cor(  log2(Ae[genelist,] +1),log2(Ae[genelist,] +1) )

pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/LineageHeatmap_C1_RNAAll.pdf",sep=""))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C1, digits = 3), border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/LineageHeatmap_C1_RNAIntgenes.pdf",sep=""))



sadsadsadsad


mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined_justEmbryonic.rds",sep=""))

uID <- as.character(mammal.combined$Cells)
#uID[which(uID=="EPI_CS4")] <- "TB_CS4"

uID[which(mammal.combined$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined$Cl05%in%c(13))] <- "EmDisc_d14"
uID[which(mammal.combined$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined$Cl05%in%c(16))] <- "Hyp_d14"
uID[which(mammal.combined$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined$Cl05%in%c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,17,18,19,20))] <- "Tb_d14"

#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("CTB_d9"))] <- "STB_d9"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("CTB_d11"))] <- "STB_d11"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("CTB_d12"))] <- "STB_d12"

#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("STB_d9"))] <- "CTB_d9"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("STB_d11"))] <- "CTB_d11"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("STB_d12"))] <- "CTB_d12"




mammal.combined$anno <- uID
Idents(mammal.combined) <- mammal.combined$Dataset

mammal.combined1 <- subset(mammal.combined,idents=c("10X Reference 1"))
mammal.combined2 <- subset(mammal.combined,idents=c("SS2 Reference 1"))
mammal.combined3 <- subset(mammal.combined,idents=c("10X Ours"))


Idents(mammal.combined1) <- mammal.combined1$anno
colind <- integer( length( levels(Idents(mammal.combined1)) )  )
for (i in 1:length( levels(Idents(mammal.combined1)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined1))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA1_Embryoniconly","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


#Idents(mammal.combined2) <- mammal.combined2$anno 
#colind <- integer( length( levels(Idents(mammal.combined2)) )  )
#for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
#  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
#}
#coluse <- cols[colind]

#p1 <- DimPlot(mammal.combined2, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,3))
#ggsave(filename=paste(saveext,"/DimRed/PCA2_Embryoniconly","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


Idents(mammal.combined3) <- mammal.combined3$anno 
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA3_Embryoniconly","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

#p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_B1-B6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)



sdsadsaadadad


uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Tb"))] <- "1) Tb"
uID[which(uID%in%c("EmDisc_d12","Embryonic","Hyp_d12","Stromal","Epithelial"))] <- "2) Not Tb"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

xgrid = make.grid(as.matrix(dat[,1:2]))

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Tb","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Epithelial"))] <- "1) Epithelial"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12","Stromal"))] <- "2) Noot Epithelial"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Epithelial","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Stromal"))] <- "1) Stromal"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12","Epithelial"))] <- "2) Not Stromal"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Stromal","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)



uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Stromal","Epithelial"))] <- "1) Maternal"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12"))] <- "2) Not Maternal"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Mat","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)



Idents(mammal.combined3) <- mammal.combined3$Tissue2
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_2_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_3_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(2,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_2_3_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined3$Tissue2)
uID[which(uID%in%c("Tb","Stromal","Epithelial"))] <- "2) Not embryonic"
uID[which(uID%in%c("EmDisc_d12","Embryonic","Hyp_d12"))] <- "1) EmDisc/Hyp"
Idents(mammal.combined3) <- uID

#Grid of points for later
make.grid = function(x, n = 1000) {
  grange = apply(x, 2, range)
  x1 = seq(from = grange[1,1], to = grange[2,1], length = n)
  x2 = seq(from = grange[1,2], to = grange[2,2], length = n)
  expand.grid(X1 = x1, X2 = x2)
}

dat <- as.data.frame(mammal.combined3[["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined3)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = .5, scale = TRUE)

xgrid = make.grid(as.matrix(dat[,1:2]))
ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined3) <- mammal.combined3$Tissue2
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C8_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


uID <- as.character(mammal.combined3$Tissue2)
uID[which(uID%in%c("Tb"))] <- "1) Tb"
uID[which(uID%in%c("EmDisc_d12","Embryonic","Hyp_d12","Stromal","Epithelial"))] <- "2) Not Tb"
Idents(mammal.combined3) <- uID

dat <- as.data.frame(mammal.combined3[["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined3)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

xgrid = make.grid(as.matrix(dat[,1:2]))

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined3) <- mammal.combined3$Tissue2
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C6_Tb","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)




uID <- as.character(mammal.combined3$Tissue2)
uID[which(uID%in%c("Stromal","Epithelial"))] <- "1) Maternal"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12"))] <- "2) Not Maternal"
Idents(mammal.combined3) <- uID

dat <- as.data.frame(mammal.combined3[["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined3)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined3) <- mammal.combined3$Tissue2
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C6_Mat","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


x


sdsadas

newID <- as.character(mammal.combined$Cells)


newID[which(mammal.combined$Cl05%in%c(0,10,12) & mammal.combined$species2=="New")] <- "Stromal fibroblasts"
newID[which(mammal.combined$Cl05%in%c(8) & mammal.combined$species2=="New")] <- "Ciliated"

newID[which(mammal.combined$Cl05%in%c(3,4,5,7) & mammal.combined$species2=="New")] <- "Unciliated epithelia 2"
newID[which(mammal.combined$Cl05%in%c(1) & mammal.combined$species2=="New")] <- "Unciliated epithelia 1"
newID[which(mammal.combined$Cl05%in%c(13) )] <- "Unciliated epithelia 3"

newID[which(mammal.combined$Cl05%in%c(17) & mammal.combined$species2=="New")] <- "Hyp_d14"

newID[which(mammal.combined$Cl05%in%c(2,15) & mammal.combined$species2=="New")] <- "CTB_d14"
newID[which(mammal.combined$Cl05%in%c(6,9,11,14) & mammal.combined$species2=="New")] <- "STB_d14"
newID[which(mammal.combined$Cl05%in%c(19) & mammal.combined$species2=="New")] <- "EmDisc_d14"

#newID[which(mammal.combined$Cl05%in%c(16) & mammal.combined$species2=="New")] <- "Unciliated epithelia 4"
newID[which(mammal.combined$Cl05%in%c(21) & mammal.combined$species2=="New")] <- "EVT"

Idents(mammal.combined) <- newID

#colind <- integer( length( levels(Idents(mammal.combined)) )  )
#for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
#  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
#}
#coluse <- cols[colind]

DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitIICol",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_SeuratIICol",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE)

sadsadsa

DimPlot(mammal.combined, cols = coluse,  reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitIICol",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, cols = coluse, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_SeuratIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)



DimPlot(mammal.combined,cols = coluse,  reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitIICol",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, cols = coluse, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_SeuratIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_ano1.rds",sep=""))

